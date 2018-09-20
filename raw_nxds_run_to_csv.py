import re
import pandas as pd
from xds import uncorrectedhkl


hklFN = "INTEGRATE.HKL"

compression = 'bz2'

hkl = uncorrectedhkl(hklFN)
hkl.data.to_csv("reflections.csv." + compression, compression = compression)

#We need to aggregate metadata now
directories = hkl.imagedata.file_name.apply(lambda x: re.match(r"^.*/", x).group()).drop_duplicates()
metadata = None
for directory in directories:
    rundata  = pd.read_csv(directory + 'metadata.tsv', delim_whitespace=True)
    rundata['file_name'] = directory + rundata['filename']
    del(rundata['filename'])
    metadata = pd.concat((metadata, rundata))


metadata = hkl.imagedata.reset_index().set_index('file_name').join(metadata.set_index('file_name'))
metadata.reset_index().set_index('Image#').to_csv("imagedata.csv." + compression, compression = compression)
