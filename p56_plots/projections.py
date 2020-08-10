"""
This displays what all the functions in loop_apps.py do.
"""
from starter2 import *
import xtra_energy
reload(loop_apps)
#core_list=[0,1, 10, 27,79, 16, 70, 44]
core_list = this_looper.core_list
frame_list=[0] #range(0,130,10)
fields=['density']  

#Look for app_test.h5 and read from that.
#If not, make a new looper and write app_test.h5.
if 'this_looper' not in dir():

    directory = '/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128/GravPotential'
    this_looper = looper.core_looper(directory= directory, derived=[xtra_energy.add_force_terms], sim_name = 'u05', out_prefix = 'plots_to_sort/test', target_frame = 125, frame_list = frame_list, core_list = core_list, fields_from_grid=['x','y','z']+fields)
    this_looper.get_target_indices(h5_name='datasets_small/u05_0125_peaklist.h5',
                                     bad_particle_list='bad_particles.h5')
    this_looper.get_tracks()

if 0:
    #project the entire domain, with particles form cores in the list plotted
    #All cores together.
    if 0:
        this_looper.get_target_indices(h5_name='datasets_small/u05_0125_peaklist.h5',
                                         bad_particle_list='bad_particles.h5')
        this_looper.get_tracks()
    loop_apps.proj_cores2(this_looper,axis_list=[0],core_list=core_list,field='density')

if 1:
    import  datasets_small.u05_core_speeds as sped
    color_dict = {}
    this_looper.frame_list=list(range(0,125,10))+[125]
    for n, regime in enumerate([ [sped.fast_cores,'r'], [sped.ok_cores,'g'], 
                                [sped.slow_cores,'b'], [sped.small_cores, 'c']]):
        color = regime[1]
        this_dict = dict( zip(regime[0], [color]*len(regime[0])))
        color_dict.update(this_dict)

    this_looper.out_prefix='plots_to_sort/proj_with_regimes'
    loop_apps.proj_with_species(this_looper,axis_list=[0],core_list=core_list,field='density', color_dict=color_dict)

