__author__ = "Zoe Wefers"
__credits__ = ["Ryan Huang"]
__version__ = "1.0.0"
__maintainer__ = "Zoe Wefers"
__email__ = "zoe.wefers@mail.mcgill.ca"

"""Purpose: Object definitions for various cellular compartments."""

import ui

class Dendrite:

    def __init__(self, img_name):
        self.img_name = img_name
        self.color = None
        self.num = None
        self.skeleton = None
        self.print = None
        self.synap_print = None
        self.outline = None
        self.synap_outline = None
        self.spot_coords = {} #{chan: [list of coordinates]} for all mrna_channels (+synapse channel)
    
    def name(self):
        if self.color is None and self.num is None:
            name = "dendrite from " + self.img_name
        elif self.num is not None and self.color is None:
            name = "dendrite #" + str(self.num) + " from " + self.img_name
        elif self.num is None:
            name = self.color + " dendrite from " + self.img_name
        else:
            name = self.color + " dendrite (#" + str(self.num) + ") from " + self.img_name
        return name


class Soma:
    def __init__(self, img_name):
        self.img_name = img_name
        self.color = None
        self.num = None
        self.print = None
        self.outline = None
        self.nucleus = None
        self.spot_coords = {} #{chan: [] for chan in mrna_channels}
 
    def name(self):
        if self.color is None and self.num is None:
            name = "soma from " + self.img_name
        elif self.num is not None and self.color is None:
            name = "soma #" + str(self.num) + " from " + self.img_name
        elif self.num is None:
            name = self.color + " soma from " + self.img_name
        else:
            name = self.color + " soma (#" + str(self.num) + ") from " + self.img_name
        return name

class Nucleus:
    def __init__(self, img_name):
        self.img_name = img_name
        self.num = None
        self.outline = None
        self.print = None

    def name(self):
        if self.num is not None:
            name = "nucleus #" + str(self.num) + " from " + self.img_name
        else:
            name = "nucleus from " + self.img_name
        return name