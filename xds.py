import pandas as pd
from subprocess import call
from os.path import exists
from glob import glob
from StringIO import StringIO
import argparse,re
import numpy as np

def get_file(textfile):
    """Helper for parser classes. 
    Parameters
    ----------
    textfile : str or file

    Returns
    -------
    f : file object
    """
    if isinstance(textfile, str):
        f = open(textfile, 'r')
    elif textfile is None:
        f = open("XDS.INP", 'r')
    else:
        f = input_file
    return f


class xparm(dict):
  def __init__(self, inFN):
    f = get_file("nXDS.INP") if input_file is None else get_file(input_file)
    lines = open(inFN, 'r').readlines()
    f.close()
    while len(lines) > 5:
        p = [lines.pop(5) for i in range(9)]
        self[p[0].strip()] = parm(p)
    self.header = lines

  def align_parms(self):
    for v in self.values():
        if v.sign == '-':
            v.flip_axes()

  def write(self, outFN):
    with open(outFN, 'w') as out:
      out.writelines(self.header)
      for k in sorted(self):
        out.writelines(self[k].lines)

class parm():
  def __init__(self, lines):
    #TODO: throw an error if len(lines) is wrong
    self.lines = lines
    self.A = np.array(map(float, lines[2].split()))
    self.B = np.array(map(float, lines[3].split()))
    self.C = np.array(map(float, lines[4].split()))
    self.vert_axis_name = None
    self.sign           = None
    self.degrees        = None
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

  def flip_axes(self):
    #X-->-X, Y-->-Y, Z-->Z
    self.A = self.A*[-1., -1., 1.]
    self.B = self.B*[-1., -1., 1.]
    self.C = self.C*[-1., -1., 1.]
    self._findvertical__()
    self._update_lines()

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



#All valid nXDS params as per the Oct 14, 2017 version
nxds_params = [
    #Job control
    "JOB=",
    "MAXIMUM_NUMBER_OF_PROCESSORS=",
    "MAXIMUM_NUMBER_OF_JOBS=",
    "CLUSTER_NODES=",
    "IMAGE_DIRECTORY=",
    "IMAGE_LIST=",
    "SECONDS=",
    "VERBOSE=",
    
    #XYCORR",
    "X-GEO_CORR=",
    "Y-GEO_CORR=",
    "BRASS_PLATE_IMAGE=",
    "HOLE_DISTANCE=",
    "MXHOLE=",
    "MNHOLE=",
    "ROFF=",
    "TOFF=",
    "STOE_CALIBRATION_PARAMETERS=",
    
    #INIT
    "BACKGROUND_RANGE=",
    "DARK_CURRENT_IMAGE=",
    "OFFSET=",
    "TRUSTED_REGION=",
    "UNTRUSTED_RECTANGLE=",
    "UNTRUSTED_ELLIPSE=",
    "UNTRUSTED_QUADRILATERAL=",
    "MINIMUM_FRACTION_OF_BACKGROUND_REGION=",
    "NBX=",
    "NBY=",
    
    #COLSPOT
    "STRONG_PIXEL=",
    "MINIMUM_NUMBER_OF_SPOTS=",
    "MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT=",
    "SPOT_MAXIMUM-CENTROID=",
    
    #POWDER
    "POWDER_CENTER_CORRECTION=",
    
    #IDXREF
    "REFINE(IDXREF)=",
    "RGRID=",
    "SEPMIN=",
    "SEPMAX=",
    "CLUSTER_RADIUS=",
    "MAXIMUM_NUMBER_OF_DIFFERENCE_VECTOR_CLUSTERS=",
    "INTEGER_ERROR=",
    "NUMBER_OF_TESTED_BASIS_ORIENTATIONS=",
    "INDEX_ERROR=",
    "INDEX_MAGNITUDE=",
    "INDEX_QUALITY=",
    "MINIMUM_FRACTION_OF_INDEXED_SPOTS=",
    
    #INTEGRATE
    "INCLUDE_RESOLUTION_RANGE=",
    "EXCLUDE_RESOLUTION_RANGE=",
    "BEAM_DIVERGENCE=",
    "BEAM_DIVERGENCE_E.S.D.=",
    "REFLECTING_RANGE_E.S.D.=",
    "NUMBER_OF_PROFILE_GRID_POINTS_ALONG_ALPHA/BETA=",
    "MINPK=",
    "CUT=",
    "PROFILE_FITTING=",
    "SIGNAL_PIXEL=",
    "BACKGROUND_PIXEL=",
    "MINIMUM_EWALD_OFFSET_CORRECTION=",
    "MINIMUM_ZETA=",
    "MAXIMUM_ERROR_OF_SPOT_POSITION=",
    
    #CORRECT
    "REFERENCE_DATA_SET=",
    "POSTREFINE=",
    "USE_REFERENCE_IN_POSTREFINEMENT=",
    "FRIEDEL'S_LAW=",
    "MERGE=",
    "MAX_CELL_AXIS_ERROR=",
    "MAX_CELL_ANGLE_ERROR=",
    "REJECT_ALIEN=",
    
    #Crystal",
    "SPACE_GROUP_NUMBER=",
    "UNIT_CELL_CONSTANTS=",
    
    #Rotation axis
    "ROTATION_AXIS=",
    "OSCILLATION_RANGE=",
    
    #Incident beam
    "X-RAY_WAVELENGTH=",
    "INCIDENT_BEAM_DIRECTION=",
    "FRACTION_OF_POLARIZATION=",
    "POLARIZATION_PLANE_NORMAL=",
    "AIR=",
    
    #Detector hardware
    "DETECTOR=",
    "NX=",
    "NY=",
    "QX=",
    "QY=",
    "MINIMUM_VALID_PIXEL_VALUE=",
    "OVERLOAD=",
    "SILICON=",
    "SENSOR_THICKNESS=",
    
    #Detector geometry
    "DIRECTION_OF_DETECTOR_X-AXIS=",
    "DIRECTION_OF_DETECTOR_Y-AXIS=",
    "ORGX=",
    "ORGY=",
    "DETECTOR_DISTANCE=",
    "DEFAULT_REFINE_SEGMENT=",
    "MINIMUM_NUMBER_OF_REFLECTIONS/SEGMENT=",
    "SEGMENT=",
    "REFINE_SEGMENT=",
    "DIRECTION_OF_SEGMENT_X-AXIS=",
    "DIRECTION_OF_SEGMENT_Y-AXIS=",
    "SEGMENT_ORGX=",
    "SEGMENT_ORGY=",
    "SEGMENT_DISTANCE=",
]

#All valid nXDS params as per the November 11, 2017 version
xds_params = [
    #Job control
    "JOB=",
    "MAXIMUM_NUMBER_OF_JOBS=",
    "MAXIMUM_NUMBER_OF_PROCESSORS=",
    "CLUSTER_NODES=",
    "SECONDS=",
    "NUMBER_OF_IMAGES_IN_CACHE=",
    "TEST=",
    
    #Detector hardware
    "DETECTOR=",
    "NX=",
    "NY=",
    "QX=",
    "QY=",
    "OVERLOAD=",
    "MINIMUM_VALID_PIXEL_VALUE=",
    "SILICON=",
    "SENSOR_THICKNESS=",
    
    #Detector distortions
    "ROFF=",
    "TOFF=",
    "STOE_CALIBRATION_PARAMETERS=",
    "BRASS_PLATE_IMAGE=",
    "HOLE_DISTANCE=",
    "MXHOLE=",
    "MNHOLE=",
    "X-GEO_CORR=",
    "Y-GEO_CORR=",
    
    #Detector noise
    "DARK_CURRENT_IMAGE=",
    "OFFSET=",
    "GAIN=",
    
    #Trusted detector region",
    "TRUSTED_REGION=",
    "UNTRUSTED_RECTANGLE=",
    "UNTRUSTED_ELLIPSE=",
    "UNTRUSTED_QUADRILATERAL=",
    "VALUE_RANGE_FOR_TRUSTED_DETECTOR_PIXELS=",
    "INCLUDE_RESOLUTION_RANGE=",
    "EXCLUDE_RESOLUTION_RANGE=",
    "MINIMUM_ZETA=",
    
    #Detector geometry",
    "DIRECTION_OF_DETECTOR_X-AXIS=",
    "DIRECTION_OF_DETECTOR_Y-AXIS=",
    "ORGX=",
    "ORGY=",
    "DETECTOR_DISTANCE=",
    "SEGMENT=",
    "REFINE_SEGMENT=",
    "DIRECTION_OF_SEGMENT_X-AXIS=",
    "DIRECTION_OF_SEGMENT_Y-AXIS=",
    "SEGMENT_ORGX=",
    "SEGMENT_ORGY=",
    "SEGMENT_DISTANCE=",
    
    #Data images",
    "NAME_TEMPLATE_OF_DATA_FRAMES=",
    "LIB=",
    "DATA_RANGE=",
    "EXCLUDE_DATA_RANGE=",
    "SPOT_RANGE=",
    "BACKGROUND_RANGE=",
    "MINIMUM_FRACTION_OF_BACKGROUND_REGION=",
    
    #Rotation axis
    "ROTATION_AXIS=",
    "OSCILLATION_RANGE=",
    "STARTING_ANGLE=",
    "STARTING_FRAME=",
    "STARTING_ANGLES_OF_SPINDLE_ROTATION=",
    "TOTAL_SPINDLE_ROTATION_RANGES=",
    "RESOLUTION_SHELLS=",
    
    #Incident beam
    "X-RAY_WAVELENGTH=",
    "INCIDENT_BEAM_DIRECTION=",
    "FRACTION_OF_POLARIZATION=",
    "POLARIZATION_PLANE_NORMAL=",
    "AIR=",
    
    #Crystal",
    "SPACE_GROUP_NUMBER=",
    "UNIT_CELL_CONSTANTS=",
    "UNIT_CELL_A-AXIS=",
    "UNIT_CELL_B-AXIS=",
    "UNIT_CELL_C-AXIS=",
    "REIDX=",
    "FRIEDEL'S_LAW=",
    "MAX_CELL_AXIS_ERROR=",
    "MAX_CELL_ANGLE_ERROR=",
    "TEST_RESOLUTION_RANGE=",
    "MIN_RFL_Rmeas=",
    "MAX_FAC_Rmeas=",
    
    #Spot finding",
    "NBX=",
    "NBY=",
    "BACKGROUND_PIXEL=",
    "STRONG_PIXEL=",
    "MAXIMUM_NUMBER_OF_STRONG_PIXELS=",
    "MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT=",
    "SPOT_MAXIMUM-CENTROID=",
    
    #Indexing",
    "RGRID=",
    "SEPMIN=",
    "CLUSTER_RADIUS=",
    "INDEX_ERROR=",
    "INDEX_MAGNITUDE=",
    "INDEX_QUALITY=",
    "INDEX_ORIGIN=",
    "MAXIMUM_ERROR_OF_SPOT_POSITION=",
    "MAXIMUM_ERROR_OF_SPINDLE_POSITION=",
    "MINIMUM_FRACTION_OF_INDEXED_SPOTS=",
    
    #Refinement",
    "REFINE(IDXREF)=",
    "REFINE(INTEGRATE)=",
    "REFINE(CORRECT)=",
    "DEFAULT_REFINE_SEGMENT=",
    "MINIMUM_NUMBER_OF_REFLECTIONS/SEGMENT=",
    
    #Peak profiles",
    "REFLECTING_RANGE=",
    "REFLECTING_RANGE_E.S.D.=",
    "BEAM_DIVERGENCE=",
    "BEAM_DIVERGENCE_E.S.D.=",
    "RELRAD=",
    "RELBET=",
    "NUMBER_OF_PROFILE_GRID_POINTS_ALONG_ALPHA/BETA=",
    "NUMBER_OF_PROFILE_GRID_POINTS_ALONG_GAMMA=",
    "CUT=",
    "DELPHI=",
    "MINPK=",
    "PROFILE_FITTING=",
    "SIGNAL_PIXEL=",
    
    #Correction factors",
    "STRICT_ABSORPTION_CORRECTION=",
    "PATCH_SHUTTER_PROBLEM=",
    "CORRECTIONS=",
    "MINIMUM_I/SIGMA=",
    "NBATCH=",
    "REFLECTIONS/CORRECTION_FACTOR=",
    "REFERENCE_DATA_SET=",
    "FIT_B-FACTOR_TO_REFERENCE_DATA_SET=",
    "WFAC1=",
    "REJECT_ALIEN=",
    "DATA_RANGE_FIXED_SCALE_FACTOR=",
]

class xdsinp(dict):
    """
    a simple dictionary based class that parses and writes XDS input files. 
    
    Args:
        input_file (str or file, optional): the XDS input file to read in. Either supply the filename as a string or a file object. Defaults to 'XDS.INP'
    """
    def __init__(self, input_file=None):
        #TODO: This pattern for all 
        f = get_file("XDS.INP") if input_file is None else get_file(input_file)
        if isinstance(input_file, str):
            f = open(input_file, 'r')
        elif input_file is None:
            f = open("XDS.INP", 'r')
        else:
            f = input_file
        self._paramlist = xds_params
        text = ''.join([re.sub(r'!.+', '', i).strip() for i in f])
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
        f = get_file("nXDS.INP") if input_file is None else get_file(input_file)
        self._paramlist = nxds_params
        text = ''.join([re.sub(r'!.+', '', i).strip() for i in f])
        self.regex = re.compile(r"({})".format('|'.join(map(re.escape, self._paramlist))))
        self._parse(text)

class dataset():
    def __init__(self, imageFN=None):
        self.imageFN = imageFN #The first image path
        self.imlist = []
        self.pattern= None
        self.dirname= None
        self._populate()

    def _populate(self):
        suffix = re.search(r"[0-9]+\..+$", self.imageFN)
        ext    = re.search(r"\..+$", self.imageFN)
        self.pattern = self.imageFN[:suffix.start()] + "[0-9]"*(ext.start() - suffix.start()) + ext.group()
        self.imlist = sorted(glob.glob(self.pattern))
        if "/" in self.imageFN:
            self.dirname = re.sub(r"\/[^\/]+$", "", self.imageFN)
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
        return nxdsinp(StringIO(xdsin.text()))

xds_datatypes = {
    "H"          : int,
    "K"          : int,
    "L"          : int,
    "IOBS"       : float,
    "SIGMA(IOBS)": float,
    "XD"         : float,
    "YD"         : float,
    "ZD"         : float,
    "RLP"        : float,
    "PEAK"       : int,
    "CORR"       : int,
    "PSI"        : float,
}



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
        self.data = pd.read_csv(f, sep=r"\s*", comment='!', names=self.fields, header=None)

    def write(self, outFN=None):
        out.write(''.join(self.header))
        with open(outFN, 'w') as out:
            out.write('\n'.join(map(self._format, self.data)))
