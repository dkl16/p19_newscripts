
#
# ALl the other things that aren't standard python packages.
#

from starter1 import *
import yt
import looper
reload(looper)
import trackage
reload(trackage)
#from core_dump import *


import looper
reload(looper)
import trackage
reload(trackage)
import tracks_read_write as trw
reload(trw)
from davetools import *

#so we can import things from sub directories
path_list= ["./data_tools", "./pdf_tools", "./testing",\
            "./tools", "./track_tools", "./trash"]
for directory in path_list:
    if directory not in os.path:
        os.path += [directory]


