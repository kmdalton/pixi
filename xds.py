import pandas as pd
from copy import copy
from xds_inp import xds_params, nxds_params
from subprocess import call
from os.path import exists,dirname,realpath
from os import devnull
from glob import glob
from StringIO import StringIO
import re
import numpy as np



class symops(dict):
    def __init__(self, libFN=None):
        if libFN is None:
            libFN = dirname(realpath(__file__)) + "/symop.lib"
        self._parse(libFN)
    def _parse(self, libFN):
        with open(libFN, 'r') as f:
            for match in re.findall(r"(?<=\n)[0-9].*?(?=\n[0-9])", '\n'+f.read(), re.DOTALL):
                k = int(match.split()[0])
                self[k] = symop(match)

class symop(dict):
    def __init__(self, text):
        self.number = int(text.split()[0])
        self.name = re.findall(r"(?<=').*?(?=')", text)[0]
        self._parse(text)
    def _parse(self, text):
        for line in text.split('\n')[1:]:
            self[line] = op(line)

class op():
    def __init__(self, text):
        self.rot_mat = np.zeros((3,3))
        ax  = { 
            'X':  np.array([1., 0., 0.]), 
            'Y':  np.array([0., 1., 0.]), 
            'Z':  np.array([0., 0., 1.]),
           }

        for i,t in enumerate(text.split(',')):
            for k,v in ax.items():
                if '-' + k in t:
                    self.rot_mat[:,i] -= v
                elif k in t:
                    self.rot_mat[:,i] += v

        self.trans = np.zeros(3)
        div = lambda x: float(x[0])/float(x[1])
        x,y,z = text.split(',')
        self.trans[0] = 0. if '/' not in x else div(re.sub(r"[^\/0-9]", "", x).split('/'))
        self.trans[1] = 0. if '/' not in y else div(re.sub(r"[^\/0-9]", "", y).split('/'))
        self.trans[2] = 0. if '/' not in z else div(re.sub(r"[^\/0-9]", "", z).split('/'))
    def __call__(self, vector):
        return np.matmul(self.rot_mat, vector)

    def translate(self, vector):
        """
        There is a decent chance this is garbage. Not necessary now, but TODO: fix this
        """
        return vector + self.trans*vector

def get_file(textfile):
    """Helper for parser classes. 
    Parameters
    ----------
    textfile : str or file

    Returns
    -------
    f : file object or None
    """
    if isinstance(textfile, str):
        return open(textfile, 'r')
    elif hasattr(textfile, 'read'):
        return textfile
    else:
        return None


class xparm(dict):
    """
    A class to parse, modify and write XPARM.nXDS files output by
    the IDXREF program of nXDS (http://nxds.mpimf-heidelberg.mpg.de/)
    
    Parameters
    ----------
    input_file : str or file
        Name of 'XPARM.nXDS' file or file object.
    
    Attributes
    ----------
    pixel_size_x : float
        The size of a pixel in the detector x-axis (mm). 
    pixel_size_y : float
        The size of a pixel in the detector y-axis (mm).
    directory : str
        The directory which contains the indexed images. 
    space_group : int
        Space group number corresponding to the images. 
    """
    def __init__(self, input_file=None):
        f = get_file("XPARM.nXDS") if input_file is None else get_file(input_file)
        lines = f.readlines()
        f.close()
        self.directory = lines[0].strip()
        length,self.n_detector_segs,self.nx,self.ny,self.pixel_size_x,self.pixel_size_y = lines[1].split()
        self.space_group_number = int(re.sub(r"[^0-9]", '', lines[2]))
        length,self.n_detector_segs,self.nx,self.ny = map(int, (length,self.n_detector_segs,self.nx,self.ny))
        self.pixel_size_x,self.pixel_size_y = map(float, (self.pixel_size_x,self.pixel_size_y))

        while len(lines) > 5:
            p = [lines.pop(5) for i in range(9)]
            self[p[0].strip()] = parm(''.join(p))
        self.header = lines[2:]

    def copy(self):
        return copy(self)

    def align_parms(self, **kw):
        """
        Flip unit cell axes to align the closest axis with the +Y direction of the detector. 
        
        See Als
        --------
        parm.flip_axes
        """
        deltaphi = kw.get('deltaphi', None)
        #k = self.keys()[0]
        #if self[k].sign == '-':
        #    self[k].flip_axes()
        ref = self[iter(self).next()]
        for k,v in self.items():
            if deltaphi is not None:
                ref = copy(self.values().next())
                suffix = re.search(r"[0-9]+?\.([0-9]|[A-z])*?$", v.lines[0]).group() 
                image_number = float(suffix.split('.')[0])
                #print -deltaphi*image_number
                phi0 = ref.phi
                ref.roty(-deltaphi*image_number)
                #print ref.phi - phi0
            if self.space_group_number in IDXAMBOPS:
                self[k] = v.align(ref, SYMOPS[self.space_group_number], IDXAMBOPS[self.space_group_number])
            else:
                self[k] = v.align(ref, SYMOPS[self.space_group_number])
            ref = self[k]
        return self

    def write(self, outFN, *images):
        """
        Save current entries to XPARM.nXDS format file.
        
        Parameters
        ----------
        outFN : str
            Name of the output file to write.
        *images : str
            Name(s) of images to write to outFN. Must match xparm.keys().
        """
        if len(images) == 0:
            images = sorted(self)
        with open(outFN, 'w') as out:
            out.write("{}\n".format(self.directory))
            out.write("{: 10d}{: 10d}{: 10d}{: 10d}{: 12.6f}{: 12.6f}\n".format(
                len(images), 
                self.n_detector_segs, 
                self.nx, 
                self.ny, 
                self.pixel_size_x, 
                self.pixel_size_y
                )
            )
            out.writelines(self.header)
            for k in images:
                out.write(k + '\n')
                out.writelines(self[k].lines[1:])

    def values(self):
        for i in self:
            yield self[i]

    def items(self):
        for i in self:
            yield i, self[i]

    def __iter__(self):
        for i in sorted(self.keys()):
            yield i

class parm():
    """
    A class representing the diffraction parameters associated with an image. 
    Based on the IDXREF output file, XPARM.nXDS.

    Parameters
    ----------
    text : str
        Text of an XPARM.nXDS file associated with a single image. The text
        contains information about unit cell and detector geometry. This is 
        documented at (http://nxds.mpimf-heidelberg.mpg.de/html_doc/nXDS_files.html#XPARM.nXDS).
        Here is an example:

            x_00009.tiff
                   1.306780       0.001248       0.001672       0.765237
              0.929059E+01  0.494593E+02  0.831171E+01
              0.602566E+02 -0.126880E+02  0.819555E+01
              0.111012E+02  0.933118E+01 -0.675948E+02
                 959.999695     959.999695      80.025002
                   1.000000       0.000000       0.000000
                   0.000000       1.000000       0.000000
                   0.000000       0.000000       1.000000
    Attributes
    ----------
    lines : list
        List of lines in the input text
    A : numpy.ndarray
        3-vector describing the direction of the unit cell A-vector
        in the lab frame (X, Y, Z)
    B : numpy.ndarray
        3-vector describing the direction of the unit cell B-vector
        in the lab frame (X, Y, Z)
    C : numpy.ndarray
        3-vector describing the direction of the unit cell C-vector
        in the lab frame (X, Y, Z)
    vert_axis_name : str
        The name of the axis which is closest to vertical
    sign : str
        "+" or "-", the direction of vert_axis_name
    vert_axis_name : str
        The name of the axis which is closest to vertical
    sign : str
        "+" or "-", the direction of vert_axis_name
    degrees : float
        the number of degrees from vertical of vert_axis_name
    phi : float
        Estimated angle between A and the detector x-axis
    imagenumber : int
        Number of image in series. Determined from the filename in the input text. 
    """
    def __init__(self, text):
        #TODO: throw an error if len(lines) is wrong
        self.lines = text.splitlines(True)
        self.imagenumber = int(re.search(r"[0-9]+\..*?$", self.lines[0].strip()).group().split('.')[0])
        self.A = np.array(map(float, self.lines[2].split()))
        self.B = np.array(map(float, self.lines[3].split()))
        self.C = np.array(map(float, self.lines[4].split()))
        self.score = None
        self.vert_axis_name = None
        self.sign           = None
        #TODO: change self.degrees to self.deflection -- fix dependent functions
        self.degrees        = None
        self.phi = None
        self.x,self.y = None,None
        self._findvertical()

    def _findvertical(self):
        a  = self.A/np.linalg.norm(self.A)
        b  = self.B/np.linalg.norm(self.B)
        c  = self.C/np.linalg.norm(self.C)
  
        pa = 180.*np.arccos(np.dot(a, [0., 1., 0.]))/np.pi
        pb = 180.*np.arccos(np.dot(b, [0., 1., 0.]))/np.pi
        pc = 180.*np.arccos(np.dot(c, [0., 1., 0.]))/np.pi
  
        na = np.abs(pa-180.)
        nb = np.abs(pb-180.)
        nc = np.abs(pc-180.)
  
        sign = np.argmin([min(pa,pb,pc), min(na,nb,nc)])
        ax   = np.argmin([min(pa,na), min(pb,nb), min(pc,nc)])
        self.vert_axis_name = ['A', 'B', 'C'][ax]
        self.sign = ['+', '-'][sign]
        self.degrees = [pa, pb, pc][ax]
        self.x = np.dot(a, [1., 0., 0.])
        self.y = np.dot(a, [0., 0., 1.])
        self.phi = np.angle(np.complex(self.x, self.y), deg=True)
        return self

    def align(self, ref, symops=None, idxambops=None):
        ref = np.array([
            ref.A, #/np.linalg.norm(ref.A), 
            ref.B, #/np.linalg.norm(ref.B),
            ref.C, #/np.linalg.norm(ref.A),
            ])

        orientations = [
            np.array([ self.A,  self.B,  self.C]),
        ]


        tmpvar=[]
        if symops is not None:
            for k,op in symops.items():
                for orientation in orientations:
                    tmpvar.append(op(orientation))
        orientations = orientations + tmpvar

        #print "SYMOPS: {}".format(symops)

        tmpvar = []
        if idxambops is not None:
            for orientation in orientations:
                for k,op in idxambops.items():
                    tmpvar.append(np.matmul(op.rot_mat, orientation))
            orientations=tmpvar

        dot = lambda x: np.sum(x * ref)
        self.score=max(map(dot, orientations))
        self.A, self.B, self.C = orientations[np.argmax(map(dot, orientations))]
        self = self._findvertical()
        return self

    def flip_axes(self, ref=None):
        """
        Flip the sign of vert_axis_name by remapping [X, Y, Z] --> [-X, -Y, Z] or [X, -Y, -Z].
        """
        if ref is None:
            ref = self.flip_vert('X')
        pX = parm(''.join(self.lines))
        pX.flip_vert('X')
        pZ = parm(''.join(self.lines))
        pZ.flip_vert('Z')
        
        if ref._dot(pX) > ref._dot(pZ):
            self = pX
        else:
            self = pZ
        return self

    def roty(self, deltaphi):
        deltaphi = np.pi*deltaphi/180.
        rot_mat = np.array([
            [ np.cos(deltaphi), 0, np.sin(deltaphi)],
            [                0, 1,                 0],
            [-np.sin(deltaphi), 0, np.cos(deltaphi)],
            ])
        self.rotate_axes(rot_mat)

    def rotate_axes(self, rot_mat):
        self.A = np.matmul(rot_mat, self.A)
        self.B = np.matmul(rot_mat, self.B)
        self.C = np.matmul(rot_mat, self.C)
        self._findvertical()

    def _dot(self, parm):
        return np.dot(self.A, parm.A) + \
               np.dot(self.B, parm.B) + \
               np.dot(self.C, parm.C) 

    def flip_vert(self, invert_axis=None):
        #X-->-X, Y-->-Y, Z-->Z
        if invert_axis is None:
            invert_axis = 'X'
        if invert_axis == 'X':
            self.A = self.A*[-1., -1., 1.]
            self.B = self.B*[-1., -1., 1.]
            self.C = self.C*[-1., -1., 1.]
        elif invert_axis == 'Z':
            self.A = self.A*[ 1., -1.,-1.]
            self.B = self.B*[ 1., -1.,-1.]
            self.C = self.C*[ 1., -1.,-1.]
        else:
            print "well shit"
        self._findvertical()
        self._update_lines()
	return self

    def _update_lines(self):
        fun = lambda x: "{: 0.6f}E{:+03d}".format(x/10**(np.floor(np.log10(np.abs(x)))+1), int(np.floor(np.log10(np.abs(x))) +1))
        self.lines[2] = " "+" ".join(map(fun, self.A)) + "\n"
        self.lines[3] = " "+" ".join(map(fun, self.B)) + "\n"
        self.lines[4] = " "+" ".join(map(fun, self.C)) + "\n"

    def __str__(self):
        return """parm instance:
            A: {}
            B: {}
            C: {}
            Axis {}: {} degrees to vertical""".format(self.A, self.B, self.C, self.vert_axis_name, self.degrees)
  
class spotnxds(dict):
    def __init__(self, spotfile=None):
        if spotfile is None:
            spotfile = 'SPOT.nXDS'
        f = get_file(spotfile)
        self.directory = f.readline().strip()
        self._populate(f.read())
        f.close()

    def _populate(self, text):
        currkey = None
        for line in text.split('\n'):
            if len(line.split()) == 1:
                currkey = line
                self[currkey] = []
            elif len(line.split()) >= 3:
                self[currkey].append(line)
        self.update({k:spotlist(v) for k,v in self.items()})

    def __iter__(self):
        for i in sorted(self.keys()):
            yield i

    def write(self, outFN):
        with open(outFN, 'w') as out:
            out.write(self.directory + '\n')
            for i in self:
                out.write(self[i].text())

    def copy(self):
        return copy(self)

class spotlist():
    def __init__(self, lines):
        self.nbg     = int(lines[0].split()[0])
        self.nstrong = int(lines[0].split()[1])
        self.nspots  = int(lines[0].split()[2])
        self.data = pd.read_csv(
            StringIO('\n'.join(lines[1:])),
            sep    ="\s*",
            engine ='python',
            names  = [
                'n',
                'X',
                'Y',
                'I',
                'H',
                'K',
                'L',
                ],
        )

    def text(self):
        text = '{: 9d}{: 8d}{: 8d}\n'.format(self.nbg, self.nstrong, self.nspots)
        text = text + '\n'.join(['{: 3d}{: 9.1f}{: 9.1f}{: 8f}.{: 7d}{: 7d}{: 7d}'.format(i.n, i.X, i.Y, i.I, i.H, i.K, i.L) for i in self.data.itertuples()])
        return text.strip()

class xdsinp(dict):
    """
    A simple dictionary based class that parses and writes XDS input files. 
    
    Parameters
    ----------
        input_file : str file, optional
            The XDS input file to read in. Either supply the filename as a string or a file object. If empty, 
            returns and xdsinp object with no keys or values.
    """
    def __init__(self, input_file=None):
        if input_file is None:
            text = ''
        else:
            text = get_file(input_file).read()
        text = re.sub(r"!.+?\n", "", text, flags=re.DOTALL)
        self._paramlist = xds_params
        self.regex = re.compile(r"({})".format('|'.join(map(re.escape, self._paramlist))))
        self._parse(text)

    def _parse(self, text):
        m = re.search(self.regex, text)
        if m is not None:
            n = re.search(self.regex, text[m.end():])
            if n is not None:
                self[m.group()] = text[m.end(): m.end()+n.start()]
                self._parse(text[n.start():])
            else:
                self[m.group()] = text[m.end():]
    def text(self):
        return "\n".join(("{}{}".format(k, self[k]) for k in self._paramlist if k in self))

    def write(self, outFN=None):
        """
        Args:
            outFN (str): output filename to save. defaults to 'XDS.INP'
        """
        if outFN is None:
            outFN = 'XDS.INP' 
        with open(outFN, 'w') as out:
            out.write(self.text())

class nxdsinp(xdsinp):
    """
    a simple dictionary based class that parses and writes XDS input files. 
    
    Args:
        inFN (str or file, optional): the XDS input file to read in. Either supply the filename as a string or a file object. Defaults to 'XDS.INP'
    """
    def __init__(self, input_file=None):
        if input_file is None:
            text = ''
        else:
            text = get_file(input_file).read()
        text = re.sub(r"!.+?\n", "", text, flags=re.DOTALL)
        self._paramlist = nxds_params
        self.regex = re.compile(r"({})".format('|'.join(map(re.escape, self._paramlist))))
        self._parse(text)
    def write(self, outFN=None):
        """
        Args:
            outFN (str): output filename to save. defaults to 'nXDS.INP'
        """
        if outFN is None:
            outFN = "nXDS.INP"
        with open(outFN, 'w') as out:
            out.write(self.text())

class dataset():
    def __init__(self, imageFN=None):
        self.imageFN = imageFN #The first image path
        self.imlist = []
        self.pattern= None
        self.dirname= None
        self._populate()

    def _populate(self):
        suffix = re.search(r"[0-9]+\..*?$", self.imageFN)
        ext    = re.search(r"\.(([0-9]|[A-z])*?)$", self.imageFN)
        self.pattern = self.imageFN[:suffix.start()] + "[0-9]"*(ext.start() - suffix.start()) + ext.group()
        self.imlist = [image(i) for i in sorted(glob(self.pattern))]
        self.dirname= self.imlist[0].dirname

    def __len__(self):
        return len(self.imlist)

    def __iter__(self):
        for i in self.imlist:
            yield i

    def generate_xdsin(self, **kw):
        verbose = kw.get('verbose', False)
        STDOUT  = open(devnull, 'w')
        if verbose:
            STDOUT = None

        if exists("XDS.INP"):
            with open("XDS.INP", "r") as f:
                backup = f.read()
        else:
            backup = None
        call(["generate_XDS.INP", self.pattern], stdout=STDOUT, stderr=STDERR)
        xdsin = xdsinp("XDS.INP")
        if backup is not None:
            with open("XDS.INP", "w") as f:
                f.write(backup)
        return xdsin

    def generate_nxdsin(self):
        xdsin = self.generate_xdsin()
        nxdsin = nxdsinp()
        nxdsin.update(xdsin)
        return nxdsin

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
    index : int
        Index of the image in the rotation series. Determined from the numbering in image_path
    """
#TODO: add support for hkl arrays. modify uncorrectedhkl to iterate by returning 'image' objects. This is major structure change and should go on a new branch. However, it will vastly simplify the way the data are processed. The dream is to be able to go like:
#for ref_image in ref_dataset:
#   ref_image.integrate()
#   for dataset in datasets:
#       for image in dataset:
#           if image.index == ref_image.index:
#               image.integrate(xds_params=image.integration_params)
#or something like that. So that we can hide away all the XDS parameters and such
    def __init__(self, image_path):
        self.path = image_path
        self.filename = re.search(r"[^\/]*?$", image_path).group() 
        if "/" in self.path:
            self.dirname = re.match(r".*\/", self.path).group()
        else:
            self.dirname = "./"
        self.index = int(re.search(r"[0-9]+\..*?$", image_path).group().split('.')[0])

    def __str__(self):
        return """xds.image object
            path    : {}
            filename: {}
            dirname : {}
            index   : {}""".format(self.path, self.filename, self.dirname, self.index)


xds_format_strings = {
    "H"          : "{: 6d}",
    "K"          : "{: 6d}",
    "L"          : "{: 6d}",
    "IOBS"       : "{:.3E}",
    "SIGMA(IOBS)": "{:.3E}",
    "XD"         : "{: 8.1f}",
    "YD"         : "{: 8.1f}",
    "ZD"         : "{: 8.1f}",
    "RLP"        : "{: 10.5f}",
    "PEAK"       : "{: 4d}",
    "CORR"       : "{: 4d}",
    "PSI"        : "{: 8.2f}",
    "ISEG"       : "{:d)",
}

class xds_ascii():
    def __init__(self, hklin):
        self.header = []
        self.fields = []
        with get_file(hklin) as f:
            line = f.readline()
            while line[0] == '!':
                if "!ITEM_" in line:
                    name = line.split("_")[-1].split("=")[0]
                    idx  = int(line.strip().split("=")[-1])
                    while len(self.fields) < val:
                        self.fields.append(None)
                    self.fields[idx-1] = name
                self.header.append(line)
                line = f.readline()

        self.fields = tuple(self.fields)
        if None in self.fields:
            raise ValueError("Missing field descriptor in XDS_ASCII file")
        self._format = lambda x: "".join([xds_format_strings[i] for i in self.fields]).format(x) 
        self.data = pd.read_csv(f, sep=r"\s*", comment='!', names=self.fields, header=None, engine='python')

    def write(self, outFN=None):
        if outFN is None:
            outFN = "nXDS_ASCII.HKL"
        out.write(''.join(self.header))
        with open(outFN, 'w') as out:
            out.write('\n'.join(map(self._format, self.data)))
            out.write("!END_OF_DATA")

class uncorrectedhkl():
    def __init__(self, hklin):
        f = get_file(hklin)
        header  = re.search(r".*?\n!END_OF_HEADER", f.read(), re.DOTALL).group()
        image_names            = re.search(r"(?<=!IMAGE_NAMES ::).*?(?=\n! \n)", header, flags=re.DOTALL).group()
        diffraction_parameters = re.search(r"(?<=!DIFFRACTION_PARAMETERS ::).*?(?=\n! \n)", header, flags=re.DOTALL).group()
        image_control          = re.search(r"(?<=!IMAGE_CONTROL :: ).*?(?=\n!NUMBER)", header, flags=re.DOTALL).group()

        image_names            = re.sub(r"[!,]", "", image_names)
        diffraction_parameters = re.sub(r"[!,]", "", diffraction_parameters)
        image_control          = re.sub(r"[!,]", "", image_control)

        self.spec    = re.search(r"!NUMBER_OF_ITEMS.*?!END_OF_HEADER", header, re.DOTALL).group()
        self.header = re.search(r".*?(?=\n!IMAGE_NAMES ::)", header, re.DOTALL).group()

        self.imagedata = pd.read_csv(StringIO(image_names), sep="\s*", index_col=0, engine='python')

        diffraction_parameter_names = [
            "Image#",
            "wavelength",
            "incident_beam_x",
            "incident_beam_y",
            "incident_beam_z",
            "A_x",
            "A_y",
            "A_z",
            "B_x",
            "B_y",
            "B_z",
            "C_x",
            "C_y",
            "C_z",
            "ORGX",
            "ORGY",
            "F",
        ]
        self.imagedata = self.imagedata.join(
            pd.read_csv(StringIO(diffraction_parameters), 
                sep="\s*", 
                skiprows=1, 
                names=diffraction_parameter_names, 
                index_col=0,
                engine='python'
            )
        )

        image_control_parameter_names = [
            "Image#",
            "#reflections",
            "#background_pixels",
            "#spot_pixels",
            "#strong_spots",
            "#spot_profiles",
            "#overloaded",
            "beam_divergence",
            "sigma1",
            "sigma2",
            "rho",
            "reflecting_range_esd_1",
            "reflecting_range_esd_2",
        ]
        self.imagedata = self.imagedata.join(
            pd.read_csv(StringIO(image_control), 
                sep="\s*", 
                skiprows=1, 
                names=image_control_parameter_names, 
                index_col=0,
                engine='python'
            )
        )

        self.fields = self.spec.split('\n')[-3][1:].split(',')
        self.fields = [i if i!='SIGMA' else 'SIGMA(IOBS)' for i in self.fields]
        f.seek(0)
        self.data = pd.read_csv(f, sep=r"\s*", comment='!', index_col=['H', 'K', 'L'], names=self.fields, header=None, engine='python')
        f.close()

    def write_xds_ascii(self, outFN):
        with open(outFN, 'w') as out:
            out.write("!FORMAT=XDS_ASCII    MERGE=FALSE    FRIEDEL'S_LAW=FALSE\n") #TODO: find out of friedel's law can ever be true here
            out.write("!Generated by PIXI")
            fields = ['H', 'K', 'L', 'IOBS', 'SIGMA(IOBS)']
            out.write("!NUMBER_OF_ITEMS_IN_EACH_DATA_RECORD={:d}\n".format(len(fields)))
            for i,field in enumerate(fields, 1):
                out.write("!ITEM_{}={:d}\n".format(field, i))
            format = lambda x: "".join([xds_format_strings[i] for i in fields]).format(x) 
            out.write('\n'.join(map(format, self.data)))
            out.write("!END_OF_DATA")

SYMOPS = symops()
IDXAMBOPS = symops(dirname(realpath(__file__)) + "/idxambops.lib")
