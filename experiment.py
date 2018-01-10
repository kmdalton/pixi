import re,xds
from os import chdir, getcwd, devnull, remove
from os.path import exists
from subprocess import call,STDOUT
from glob import glob

NULL = open(devnull, 'w')


class experiment(list):
    """
    What's the point of this class? It's a container for other classees that compose the experiment. It will have the root methods for batch integration. Experiments contain crystals. Crystals contain images. It will be used to keep track of what imseries correspond to what obervable. For instance, use this class to keep track of replicates of a given pump probe delay. 
    """
    def integrate(self, reference, nxdsinp=None):
        for crystal in self:
            crystal[reference].integrate(nxdsinp)
            crystal.sync_integration_parameters(reference)

    def ratio(self, numerator, denominator):
        """
        Later we will replace this an arbitrary intensity algebra. But for now, we will just calculate the ratio of image series.

        Parameters
        ----------
        numerator : iterable
            An iterable containing strings which are keys in the crystal objects. These will be pooled to estimate the numerator. 
        denominator : iterable
            An iterable containing strings which are keys in the crystal objects. These will be pooled to estimate the denominator. 
        """
        pass

class crystal(dict):
    """
    A crystal is a dictionary of measurements at the same phi step.
    """
#Premature optimization is ___________
    def sync_integration_parameters(self, reference):
        for im_series in self.values():
            for image in im_series:
                for ref in self[reference]:
                    if image.imagenumber == ref.imagenumber:
                        image.nxdsin = ref.nxdsin
                        image.xparm  = ref.xparm

class image_series(list):
    """
    A collection of image objects.
    Parameters
    ----------
    imageFN : str (optional)
        The first image filename in a series. Defaults to none to create an empty object.

    Attributes
    ----------
    imageFN : str
        The first image file in the series.
    pattern : str
        A glob pattern corresponding to the images
    dirname : str
        The directory in which the images reside
    """
    def __init__(self, imageFN=None):
        self.imageFN = imageFN #The first image path
        self.pattern= None
        self.dirname= None
        self._populate()

    def _populate(self):
        suffix = re.search(r"[0-9]+\..*?$", self.imageFN)
        ext    = re.search(r"\.(([0-9]|[A-z])*?)$", self.imageFN)
        self.pattern = self.imageFN[:suffix.start()] + "[0-9]"*(ext.start() - suffix.start()) + ext.group()
        [self.append(image(i)) for i in sorted(glob(self.pattern))]
        self.dirname= self[0].dirname

    def __str__(self):
        return ' '.join([i.filename for i in self])

    def scale(self, scaledict):
        """
        Updates image.scale for each image object.

        Parameters
        ----------
        scaledict : dict
            A dictionary like {imagefilename: float}. 
        """
        for image in self:
            if image.filename in scaledict:
                image.scale = scaledict[image.filename]

    def generate_xdsin(self, **kw):
        """
        Use Kay Diederich's generate_XDS.INP bash script to generate an XDS input file from the images in the dataset. 

        Returns
        -------
        xdsin : xds.xdsinp
        """
        verbose = kw.get('verbose', False)
        stdout  = NULL
        stderr  = STDOUT
        if verbose:
            stdout = None
            stderr = None

        if exists("XDS.INP"):
            with open("XDS.INP", "r") as f:
                backup = f.read()
        else:
            backup = None
        call(["generate_XDS.INP", self.pattern], stdout=stdout, stderr=stderr)
        xdsin = xds.xdsinp("XDS.INP")
        if backup is not None:
            with open("XDS.INP", "w") as f:
                f.write(backup)
        return xdsin

    def generate_nxdsin(self, **kw):
        """
        Use Kay Diederich's generate_XDS.INP bash script to generate an nXDS input file from the images in the dataset. This is slightly suspect but mostly okay. 

        Returns
        -------
        nxdsin : xds.nxdsinp
        """
        xdsin = self.generate_xdsin(**kw)
        nxdsin = xds.nxdsinp()
        nxdsin.update(xdsin)
        return nxdsin

    def write_imlist(self, filename):
        """
        Write a list of image filenames for nXDS
        """
        with open(filename, 'w') as out:
            out.write('\n'.join([i.filename for i in self]))

    def integrate(self, nxdsin=None, **kw):
        """
        Use nXDS to integrate this series in batch. Copy integration parameters to images. 
        Parameters
        ----------
        nxdsin : xds.nxdsinp (optional)
            Optional parameters for integration. Defaults to self.generate_nxdsin output.

        Kwargs
        ------
        verbose : bool (optional)
            If True, do not redirect nXDS output os.dev_null. Default is False.
        """
        verbose = kw.get('verbose', False)
        stdout  = NULL
        stderr  = STDOUT
        if verbose:
            stdout = None
            stderr = None

        if nxdsin is None:
            nxdsin = self.generate_nxdsin()
        self.write_imlist("LISTIM")
        nxdsin['JOB='] = 'XYCORR INIT COLSPOT POWDER IDXREF INTEGRATE'
        nxdsin['IMAGE_LIST='] = "LISTIM"
        nxdsin['IMAGE_DIRECTORY='] = self.dirname
        nxdsin.write()
        call(['nxds_par'], stdout=stdout, stderr=stderr)

        hkl = xds.uncorrectedhkl("INTEGRATE.HKL")
        nxdsin['JOB='] = " INTEGRATE"
        xparm = xds.xparm("XPARM.nXDS")
        for i,imagedatum in hkl.imagedata.iterrows():
            filename = imagedatum.file_name
            imnxdsin = xds.nxdsinp()
            imnxdsin.update(nxdsin)
            imnxdsin['BEAM_DIVERGENCE='] = imagedatum.beam_divergence
            imnxdsin['BEAM_DIVERGENCE_E.S.D.=']  = "{} {}".format(imagedatum.sigma1, imagedatum.sigma2)
            imnxdsin['REFLECTING_RANGE_E.S.D.='] = "{} {}".format(imagedatum.reflecting_range_esd_1, imagedatum.reflecting_range_esd_2)
            self[filename].nxdsin = nxdsin
#TODO: There is a better way to do this but requires modifying xds.xparm
            self[filename].xparm = xparm.copy()
            for k in self[filename].xparm:
                if k != filename:
                    del self[filename].xparm[k]

class image():
    """
    Representation of an image file. Right now, this class has no methods, but the idea is to implement per image methods in the future (ie image.integrate, image.index, etc). 

    Parameters
    ----------
    image_path : str
        Full path of the image file

    Attributes
    ----------
    path : str
        Full path of image file
    filename : str
        Name of the image file
    dirname : str
        Directory which contains the image file
    imagenumber : int
        Index of the image in the rotation series. Determined from the numbering in image_path
    xparm : xds.xparm
        xparm instance giving the A matrix for this image. 
    scale : float
        Optional scale parameter for image. Defaults to 1.
    """
    def __init__(self, image_path):
        self.path     = image_path
        self.nxdsin   = xds.nxdsinp()
        self.xparm    = None
        self.hkl      = None
        self.scale    = 1.
        self.filename = re.search(r"[^\/]*?$", image_path).group() 
        if "/" in self.path:
            self.dirname = re.match(r".*\/", self.path).group()
        else:
            self.dirname = "./"
        self.imagenumber = int(re.search(r"[0-9]+\..*?$", image_path).group().split('.')[0])

    def __str__(self):
        return """xds.image object
            path    : {}
            filename: {}
            dirname : {}""".format(self.path, self.filename, self.dirname)

    def index(self, nxdsin=None, **kw):
        """
        Index the image using nXDS. Updates self.xparm accordingly.

        Parameters
        ----------
        xdsin : xds.nxdsinp (optional)
            Specify nXDS parameter file object for indexing run. 
        """
        verbose = kw.get('verbose', False)
        stdout  = NULL
        stderr  = STDOUT
        if verbose:
            stdout,stderr = None,None

        nxdsin = kw.get('nxdsin', xds.nxdsinp())
        nxdsin.update(self.nxdsin)
        xdsinp['JOB='] = " XYCORR INIT COLSPOT POWDER IDXREF"
        xdsinp['IMAGE_LIST='] = "LISTIM"
        xdsinp['IMAGE_DIRECTORY='] = self.dirname
        nxdsin.write('nXDS.INP')
        with open('LISTIM', 'w') as out:
            out.write(self.filename)

    def integrate(self, nxdsin=None, **kw):
        """
        Index the image using nXDS. Updates self.xparm accordingly.

        Parameters
        ----------
        xdsin : xds.nxdsinp (optional)
            Specify nXDS parameter file object for indexing run. 
        """
#TODO: check that self.xparm is not None
        verbose = kw.get('verbose', False)
        stdout  = NULL
        stderr  = STDOUT
        if verbose:
            stdout,stderr = None,None

        nxdsin = kw.get('nxdsin', xds.nxdsinp())
        nxdsin.update(self.nxdsin)
        xdsinp['JOB='] = " INTEGRATE"
        xdsinp['IMAGE_LIST='] = "LISTIM"
        xdsinp['IMAGE_DIRECTORY='] = self.dirname
        nxdsin.write('nXDS.INP')
        with open('LISTIM', 'w') as out:
            out.write(self.filename)
        self.xparm.write()

        remove('INTEGRATE.HKL')
        call(['nxds_par'], stdout=stdout, stderr=stderr)
        if exists('INTEGRATE.HKL'):
            self.hkl = xds.uncorrectedhkl('INTEGRATE.HKL')
