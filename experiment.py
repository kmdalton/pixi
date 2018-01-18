import re,xds
from copy import copy,deepcopy
import numpy as np
from matplotlib import pyplot as plt
from os import chdir, getcwd, devnull, remove
from os.path import exists
from subprocess import call,STDOUT
from glob import glob
import pandas as pd

NULL = open(devnull, 'w')


class experiment(list):
    """
    What's the point of this class? It's a container for other classees that compose the experiment. It will have the root methods for batch integration. Experiments contain crystals. Crystals contain images. It will be used to keep track of what imseries correspond to what obervable. For instance, use this class to keep track of replicates of a given pump probe delay. 
    """
    def integrate(self, reference, nxdsinp=None):
        for crystal in self:
            crystal.integrate(reference, nxdsinp)

    def doeke_scaling(self, numerator, denominator):
        """
        Later we will replace this an arbitrary intensity algebra. But for now, we will just calculate the ratio of image series.

        Parameters
        ----------
        numerator : iterable
            An iterable containing strings which are keys in the crystal objects. These will be pooled to estimate the numerator. 
        denominator : iterable
            An iterable containing strings which are keys in the crystal objects. These will be pooled to estimate the denominator. 
        """
        gamma = None
        for crystal in self:
            for image in crystal:
                num = [image[i] for i in numerator   if image[i].hkl is not None]
                den = [image[i] for i in denominator if image[i].hkl is not None]
                if len(num) > 0 and len(den) > 0:
                    omega = sum([i.scale for i in den]) / sum([i.scale for i in num])
                    g = sum([i.hkl for i in num]) * float(len(den)) / (sum([i.hkl for i in den]) * float(len(num)))
                    g = g * omega
                    gamma = pd.concat((gamma, g.data))
        return gamma

class crystal(dict):
    def __init__(self):
        self.indices = []
    """
    A crystal is a container for image_series which allows easier iteration over equivalent images. Dict-like assignment is supported. However, keys may only be strings and values may only be image_series.
    """
    def sync_integration_parameters(self, reference):
        for imdict in self:
            for v in imdict.values():
                v.nxdsin = imdict[reference].nxdsin
                v.xparm  = imdict[reference].xparm.copy()
                v.xparm.directory = v.dirname

    def integrate(self, reference, nxdsin):
        self[reference].integrate(nxdsin)
        self.sync_integration_parameters(reference)
        for imdict in self:
            for v in imdict.values():
                if v.xparm is not None:
                    v.integrate()

    def __setitem__(self, k, v):
        if isinstance(v, image_series) and isinstance(k, str):
            for im in v:
                if im.imagenumber not in self.indices:
                    self.indices.append(im.imagenumber)
            self.indices = sorted(self.indices)
            dict.__setitem__(self, k, v)
        elif not isinstance(k, str):
            raise TypeError("key has Type {}, but crystal values can only have type image_series")
        elif not isinstance(v, image_series):
            raise TypeError("Value has Type {}, but crystal values can only have type image_series")

    def __iter__(self):
        for index in self.indices:
            yield {k: v.get_image_by_index(index) for k,v in self.items()}

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
        self.indices = {}
        self._populate()

    def get_image_by_index(self, index):
        """
        Parameters
        ----------
        index : int
            The index of the image in the rotation series you want to extract.
        Returns
        -------
        im : image
            image object corresponding to index.
        """
        if index in self.indices:
            return self[self.indices[index]]
        else:
            return None

    def _populate(self):
        suffix = re.search(r"[0-9]+\..*?$", self.imageFN)
        ext    = re.search(r"\.(([0-9]|[A-z])*?)$", self.imageFN)
        self.pattern = self.imageFN[:suffix.start()] + "[0-9]"*(ext.start() - suffix.start()) + ext.group()
        [self.append(image(i)) for i in sorted(glob(self.pattern))]
        for i, im in enumerate(self):
            self.indices[im.imagenumber] = i
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
        return self

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
        align   = kw.get('align', True)
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
        if exists("INTEGRATE.HKL"):
            remove("INTEGRATE.HKL")
        call(['nxds_par'], stdout=stdout, stderr=stderr)

        imagedata = None
        if exists("INTEGRATE.HKL"):
            imagedata = xds.uncorrectedhkl("INTEGRATE.HKL").imagedata
        nxdsin['JOB='] = " INTEGRATE"
        xparm = xds.xparm("XPARM.nXDS")
        if align:
            xparm.align_parms()
        spotnxds = xds.spotnxds("SPOT.nXDS")
        for im in self:
            if im.filename in spotnxds:
                im.spot = spotnxds.copy()
                im.spot.clear()
                im.spot[im.filename] = spotnxds[im.filename]
            if imagedata is not None:
                integration_params = imagedata.loc[imagedata['file_name']==im.filename]
                im.nxdsin = xds.nxdsinp()
                im.nxdsin.update(nxdsin)
                if len(integration_params) == 1:
                    im.nxdsin['BEAM_DIVERGENCE='] = float(integration_params.beam_divergence)
                    im.nxdsin['BEAM_DIVERGENCE_E.S.D.=']  = "{} {}".format(float(integration_params.sigma1), float(integration_params.sigma2))
                    im.nxdsin['REFLECTING_RANGE_E.S.D.='] = "{} {}".format(float(integration_params.reflecting_range_esd_1), float(integration_params.reflecting_range_esd_2))
                    im.xparm = xparm.copy()
                    im.xparm.clear()
                    im.xparm[im.filename] = xparm[im.filename]
        return self

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
    spot : xds.spot
        Spotlist from nxds COLSPOT
    """
    def __init__(self, image_path=None):
        self.path     = image_path
        self.nxdsin   = xds.nxdsinp()
        self.xparm    = None
        self.hkl      = None
        self.scale    = 1.
        self.spot = None #xds.spotlist
        self.filename = None
        self.dirname = None
        self.imagenumber = None
        if image_path is not None:
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

        n = xds.nxdsinp()
        if self.nxdsin is not None:
            n.update(self.nxdsin)
        n.update(kw.get('nxdsin', xds.nxdsinp()))
        n['JOB='] = " COLSPOT IDXREF"
        n['IMAGE_LIST='] = "LISTIM"
        n['IMAGE_DIRECTORY='] = self.dirname
        n.write('nXDS.INP')
        if exists('SPOT.nXDS'):
            remove('SPOT.nXDS')
        with open('LISTIM', 'w') as out:
            out.write(self.filename)
        call(['nxds_par'], stdout=stdout, stderr=stderr)
        if exists('SPOT.nXDS'):
            self.spot = xds.spotnxds('SPOT.nXDS')

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

        n = xds.nxdsinp()
        if self.nxdsin is not None:
            n.update(self.nxdsin)
        n.update(kw.get('nxdsin', xds.nxdsinp()))
        n['JOB='] = " INTEGRATE"
        n['IMAGE_LIST='] = "LISTIM"
        n['IMAGE_DIRECTORY='] = self.dirname
        n.write('nXDS.INP')
        with open('LISTIM', 'w') as out:
            out.write(self.filename)
        self.xparm.write('XPARM.nXDS')

        if exists('INTEGRATE.HKL'):
            remove('INTEGRATE.HKL')
        call(['nxds_par'], stdout=stdout, stderr=stderr)
        if exists('INTEGRATE.HKL'):
            self.hkl = xds.uncorrectedhkl('INTEGRATE.HKL')

    def plot(self, **kw):
        #All pixels greater than this percentile get the same value
        percentile_cutoff = kw.get("percentile_cutoff", 0.999)
        cmap_name         = kw.get("cmap", "gray_r")
        data = plt.imread(self.dirname + self.filename)
        cutoff = np.sort(data.flatten())[int(percentile_cutoff*data.size)]
        data[data>=cutoff]=cutoff
        plt.imshow(data, cmap=cmap_name)
        X,Y = self.spot[self.filename].data.X,self.spot[self.filename].data.Y
        plt.scatter(X, Y, s=100, facecolors='None', edgecolors='w')
