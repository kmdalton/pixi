from subprocess import call
from os.path import exists
from glob import glob
from StringIO import StringIO
import argparse,re


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
        if isinstance(input_file, str):
            f = open(input_file, 'r')
        elif input_file is None:
            f = open("XDS.INP", 'r')
        else:
            f = input_file
        self._paramlist = xds_params
        text = ''.join([re.sub(r'!.+', '', i).strip() for i in f])
        self.regex = re.compile(r"({})".format('|'.join(map(re.escape, self._paramlist)))
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
        if isinstance(input_file, str):
            f = open(input_file, 'r')
        elif input_file is None:
            f = open("XDS.INP", 'r')
        else:
            f = input_file
        self._paramlist = nxds_params
        text = ''.join([re.sub(r'!.+', '', i).strip() for i in f])
        self.regex = re.compile(r"({})".format('|'.join(map(re.escape, self._paramlist)))
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
#TODO: populate self.dirname with the directory the images are in. 

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

