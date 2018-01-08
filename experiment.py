class experiment()
    """
    What's the point of this class? It's a container for other classees that compose the experiment. It will have the root methods for batch integration. Experiments contain crystals. Crystals contain images. 
    """
    def __init__(self):
        pass

class crystal():
    """
    A crystal is a list of phi series.
    """
    def __init__(self):
        pass

class image():
    """
    A single diffraction images

    Attributes
    ----------
    directory : str
        Name of the image directory
    filename : str
        Name of the image file
    index : integer
        Index of the image in the rotation series. This is just determined from the file name. Don't overthink it
    """
    def __init__(self, filename):
        self.
