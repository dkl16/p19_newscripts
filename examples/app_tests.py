"""
This displays what all the functions in loop_apps.py do.
"""
from starter2 import *
import xtra_energy
reload(loop_apps)
core_list=[79, 80]
frame_list=[40]
fields=['density']  
if 'this_looper' not in dir():

    directory = '/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128/GravPotential'
    this_looper = looper.core_looper(directory= directory,
                                     derived=[xtra_energy.add_test_energies],
                                     sim_name = 'u05',
                                     out_prefix = 'plots_to_sort/test',
                                     target_frame = 125,
                                     frame_list = frame_list,
                                     core_list = core_list,
                                     fields_from_grid=['x','y','z']+fields
                                  )
    this_looper.get_target_indices(h5_name='datasets_small/u05_0125_peaklist.h5',
                                     bad_particle_list='bad_particles.h5')
    this_looper.get_tracks()

if 0:
    #print the centroid of each core.
    loop_apps.print_centroid(this_looper)
if 0:
    #project the entire domain, with particles form cores in the list plotted
    #All cores together.
    loop_apps.proj_cores(this_looper,axis_list=[0],core_list=core_list,field='density')

if 0:
    #this takes a list of particle indices.  
    #mostly for debugging purposes.
    loop_apps.select_particles(this_looper,axis_list=[0],
                               these_particles=this_looper.target_indices[ core_list[0] ])

if 1:
    #Plot each core individually, zoomed in to the core itself.
    this_looper.core_list=[79]
    reload(looper)
    loop_apps.core_proj_follow(this_looper,field='density',axis_list=[0])
