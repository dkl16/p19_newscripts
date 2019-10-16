from starter2 import *

import xtra_energy
if 'this_looper' not in dir():

    directory = '/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128/GravPotential'
    this_looper = looper.core_looper(directory= directory,
                                     derived=[xtra_energy.add_test_energies],
                                     sim_name = 'u05',
                                     out_prefix = 'plots_to_sort/test',
                                     target_frame = 125,
                                     frame_list = [40],
                                     core_list = [79],
                                     fields_from_grid=['x','y','z']
                                  )
    this_looper.get_target_indices(h5_name='datasets_small/u05_0125_peaklist.h5',
                                     bad_particle_list='bad_particles.h5')
    #this_looper.get_tracks()

loop_apps.proj_onecore(this_looper,axis_list=[0],core_list=[79],field='density')
