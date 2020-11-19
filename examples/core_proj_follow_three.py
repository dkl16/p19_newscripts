from starter2 import *

reload(looper)
import three_loopers as TL
reload(loop_apps)
import loop_tools
reload(loop_tools)
kill = []

for looper666 in [TL.looper3]:
    for i,nc in enumerate(looper666.core_list):
        if nc not in looper666.target_indices.keys():
            kill.append(i)
    for i in kill[::-1]:
        looper666.core_list.pop(i)
    loop_apps.core_proj_follow(looper666,axis_list=[0])
