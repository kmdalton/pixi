import re
from glob import glob


class dataset():
    def __init__(self, imageFN=None):
        self.imageFN = imageFN #The first image path
        self.imlist = []
        self.pattern= None
        self._populate()

    def _populate(self):
        suffix = re.search(r"[0-9]+\..+$", self.imageFN)
        ext    = re.search(r"\..+$", self.imageFN)
        self.pattern = self.imageFN[:suffix.start()] + "[0-9]"*(ext.start() - suffix.start()) + ext.group()
        self.imlist = sorted(glob.glob(self.pattern))

    def __len__(self):
        return len(self.imlist)

    def __iter__(self):
        for i in self.imlist:
            yield i

