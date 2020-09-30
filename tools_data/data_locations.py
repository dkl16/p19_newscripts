"""
A container for code portability.
Put
export machine = my_machine_name
in your .bashrc
Shipsterns and Mullaghmore are two of my computers, copy that.
"""
from starter2 import *
#enzo_directory = '/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128/GravPotential'
sim_u05 = '/archive2/dcollins4096/Paper19/u05-r4-l4-128-Beta0.2'
sim_u10 = '/archive2/dcollins4096/Paper19/u10_r4_l4_128-Beta2/GravPotential'
sim_u11 = '/archive2/dcollins4096/Paper19/u11_r4_l4_128-Beta20/GravPotential'
sims={'u05':sim_u05,'u10':sim_u10,'u11':sim_u11}

peaks_u05 = 'datasets_small/u05_0125_peaklist.h5'
peaks_u10 = 'datasets_small/u10_0082_peaklist.h5'
peaks_u11 = 'datasets_small/u11_0088_peaklist.h5'
peak_list = {'u05':peaks_u05,'u10':peaks_u10,'u11':peaks_u11}


target_frames={'u05':125,'u10':82,'u11':88}

u05_sixteen_frame ='/scratch1/dcollins/Paper19/Datasets/track_indfix_sixteenframe/*h5'
u10_every_ten = "/scratch1/dcollins/Paper19/Datasets/u10_every_ten/u10_all_primitives_primitives_c*_nXXX0.h5"
u11_every_ten = "/scratch1/dcollins/Paper19/Datasets/u11_every_ten/u11_all_primitives_primitives_c*_nXXX0.h5"
every_ten = {'u05':u05_sixteen_frame,'u10':u10_every_ten,'u11':u11_every_ten}



bad_particles_u05='datasets_small/u05_bad_particles.h5'
bad_particles={'u05':bad_particles_u05,'u10':None,'u11':None}

n_particles={'u05':'datasets_small/u05_n_particles.txt',
             'u10':'datasets_small/u10_n_particles.txt',
             'u11':'datasets_small/u11_n_particles.txt'}



#sixteen_frame ='/scratch1/dcollins/Paper19/Datasets/track_indfix_sixteenframe'
snapshot_base ='/scratch1/dcollins/Paper19/Datasets/'
output_directory = "./plots_to_sort/"
if 'machine' in os.environ:
    if os.environ['machine'] in ['mullaghmore']:
        #enzo_directory = '/Users/dcollins/scratch/P19/u05-r4-l4-128/GraviPotential'
        u05_snapshot_location = '/home/dcollins/scratch/Paper19/track_index_fix'

