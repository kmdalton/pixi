import pandas as pd
from _xds_inp import xds_params, nxds_params
from subprocess import call
from os.path import exists
from glob import glob
from StringIO import StringIO
import re
import numpy as np

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
        self.space_group = int(lines[2].strip())
        length,self.n_detector_segs,self.nx,self.ny,self.pixel_size_x,self.pixel_size_y = lines[1].split()

        length,self.n_detector_segs,self.nx,self.ny = map(int, (length,self.n_detector_segs,self.nx,self.ny))
        self.pixel_size_x,self.pixel_size_y = map(float, (self.pixel_size_x,self.pixel_size_y))

        while len(lines) > 5:
            p = [lines.pop(5) for i in range(9)]
            self[p[0].strip()] = parm(''.join(p))
        self.header = lines[2:]

    def align_parms(self):
        """
        Flip unit cell axes to align the closest axis with the +Y direction of the detector. 
        
        See Also
        --------
        parm.flip_axes
        """
        for k,v in self.items():
            if v.sign == '-':
                self[k] = v.flip_axes()


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
                out.writelines(self[k].lines)

    def values(self):
        for i in self:
            yield self[i]

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
    """
    def __init__(self, text):
        #TODO: throw an error if len(lines) is wrong
        self.lines = text.splitlines(True)
        self.A = np.array(map(float, self.lines[2].split()))
        self.B = np.array(map(float, self.lines[3].split()))
        self.C = np.array(map(float, self.lines[4].split()))
        self.vert_axis_name = None
        self.sign           = None
	#TODO: change self.degrees to self.deflection -- fix dependent functions
        self.degrees        = None
	self.phi = None
        self._findvertical__()
  
    def _findvertical__(self):
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
        self.phi = 180.*np.arccos(np.dot(a, [1., 0., 0.]))/np.pi

    def flip_axes(self):
        """
        Flip the sign of vert_axis_name by remapping [X, Y, Z] --> [-X, -Y, Z].
        """
	pX = parm(''.join(self.lines))
	pX.flip_vert('X')
	pZ = parm(''.join(self.lines))
	pZ.flip_vert('Z')
	if pX.phi < pZ.phi:
            self = pX
            print 'x'
        else:
            self = pZ
            print 'z'
        return self

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
        self._findvertical__()
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
        self.imlist = sorted(glob(self.pattern))
	print self.imlist
        #self.imlist = [re.search(r"(?<=(\/)[^\/]*?$", i).group() for i in self.imlist]
        self.imlist = [re.search(r"[^\/]*?$", i).group() for i in self.imlist]
        if "/" in self.imageFN:
            self.dirname = re.match(r".*\/", self.imageFN).group()
        else:
            self.dirname = "./"

    def __len__(self):
        return len(self.imlist)

    def __iter__(self):
        for i in self.imlist:
            yield i

    def generate_xdsin(self):
        if exists("XDS.INP"):
            with open("XDS.INP", "r") as f:
                backup = f.read()
        else:
            backup = None
        call(["generate_XDS.INP", self.pattern])
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
        self.data = pd.read_csv(f, sep=r"\s*", comment='!', names=self.fields, header=None, engine='python')
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

