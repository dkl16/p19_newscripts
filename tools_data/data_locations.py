"""
A container for code portability.
Put
export machine = my_machine_name
in your .bashrc
Shipsterns and Mullaghmore are two of my computers, copy that.
"""
from starter2 import *

#
# guess the computer
#

machine = os.environ.get('machine',None)
if machine is None:
    if 'machine' in os.environ:
        machine = os.environ['machine']
    elif os.path.exists("/scratch1"):
        machine = 'Nazare'
    elif  os.path.exists("/data/cb1"):
        machine = 'cloudbreak'
    else:
        print("data_locations.py Bad error: cannot detrmine machine")

output_directory = "./plots_to_sort"
sim_u14=None        
if machine == 'Nazare':
    sim_u05 = '/archive2/dcollins4096/Paper19/u05-r4-l4-128-Beta0.2'
    sim_u10 = '/archive2/dcollins4096/Paper19/u10_r4_l4_128-Beta2/GravPotential'
    sim_u11 = '/archive2/dcollins4096/Paper19/u11_r4_l4_128-Beta20/GravPotential'

    u05_every_ten ='/scratch1/dcollins/Paper19/Datasets/track_indfix_sixteenframe/*h5'
    u10_every_ten = "/scratch1/dcollins/Paper19/Datasets/u10_every_ten/u10_all_primitives_primitives_c*_nXXX0.h5"
    u11_every_ten = "/scratch1/dcollins/Paper19/Datasets/u11_every_ten/u11_all_primitives_primitives_c*_nXXX0.h5"


elif machine == 'cloudbreak':
    sim_u05 = '/data/cb1/Projects/P19_CoreSimulations/u05-r4-l4-128-Beta0.2/GravPotential'
    sim_u10 = '/data/cb1/Projects/P19_CoreSimulations/u10_r4_l4_128-Beta2/GravPotential'
    sim_u11 = '/data/cb1/Projects/P19_CoreSimulations/u11_r4_l4_128-Beta20/GravPotential'
    sim_u14 = '/data/cb1/Projects/P19_CoreSimulations/u14_sphere'

    u05_every_ten = '/data/cb1/Projects/P19_CoreSimulations/CoreSets/u05_every_ten/*h5'
    u10_every_ten = "/data/cb1/Projects/P19_CoreSimulations/CoreSets/u10_every_ten/u10_all_primitives_primitives_c*_nXXX0.h5"
    u11_every_ten = "/data/cb1/Projects/P19_CoreSimulations/CoreSets/u11_every_ten/u11_all_primitives_primitives_c*_nXXX0.h5"

    sim_u101 = '/data/cb1/Projects/P19_CoreSimulations/u05-r4-l4-128-Beta0.2/GravPotential'
    sim_u102 = '/data/cb1/Projects/P19_CoreSimulations/u10_r4_l4_128-Beta2/GravPotential'  
    sim_u103 = '/data/cb1/Projects/P19_CoreSimulations/u11_r4_l4_128-Beta20/GravPotential' 
    u101_every_ten = "/data/cb1/Projects/P19_CoreSimulations/CoreSets/u101_every_ten/u101_all_primitives_primitives_c*_nXXX0.h5"
    u102_every_ten = "/data/cb1/Projects/P19_CoreSimulations/CoreSets/u102_every_ten/u102_all_primitives_primitives_c*_nXXX0.h5"
    u103_every_ten = "/data/cb1/Projects/P19_CoreSimulations/CoreSets/u103_every_ten/u103_all_primitives_primitives_c*_nXXX0.h5"

    #The 200 series continues u05, u10, and u11 to 1tff, 
    #and ties the analysis to cores found at 1 tff.
    sim_u201 = '/data/cb1/Projects/P19_CoreSimulations/u05-r4-l4-128-Beta0.2/GravPotential'
    sim_u202 = '/data/cb1/Projects/P19_CoreSimulations/u202-Beta2/GravPotential'  
    sim_u203 = '/data/cb1/Projects/P19_CoreSimulations/u203-Beta20/GravPotential' 

    u201_every_ten = "/data/cb1/Projects/P19_CoreSimulations/CoreSets/u201/*h5"
    u202_every_ten = "/data/cb1/Projects/P19_CoreSimulations/CoreSets/u202/*h5"
    u203_every_ten = "/data/cb1/Projects/P19_CoreSimulations/CoreSets/u203/*h5"
else:
    sim_u05 = None
    sim_u10 = None
    sim_u11 = None

    u05_every_ten = None
    u10_every_ten = None
    u11_every_ten = None

    sim_u101 = None
    sim_u102 = None
    sim_u103 = None
    u101_every_ten = None
    u102_every_ten = None
    u103_every_ten = None

    energy_vorticity_files = {'u05':"/data/cb1/Projects/P19_CoreSimulations/CoreSets/energy_vorticity/u05/u05_energy_vorticity_primitives_c????_nXXX0.h5",
                              'u10':"/data/cb1/Projects/P19_CoreSimulations/CoreSets/energy_vorticity/u10/u10_energy_vorticity_primitives_c????_nXXX0.h5",
                              'u11':"/data/cb1/Projects/P19_CoreSimulations/CoreSets/energy_vorticity/u11/u11_energy_vorticity_primitives_c????_nXXX0.h5"}
    energy_vorticity_files = {'u05':"./u05_energy_vorticity_primitives_c????_nXXX0.h5",
                              'u10':"./u10_energy_vorticity_primitives_c????_nXXX0.h5",
                              'u11':"./u11_energy_vorticity_primitives_c????_nXXX0.h5"}


sims={'u05':sim_u05,'u10':sim_u10,'u11':sim_u11,'u101':sim_u101,'u102':sim_u102,'u103':sim_u103}
sims.update({'u201':sim_u201,'u202':sim_u202,'u203':sim_u203})
sims['u14']=sim_u14


peaks_u05 = 'datasets_small/u05_0125_peaklist.h5'
peaks_u10 = 'datasets_small/u10_0082_peaklist.h5'
peaks_u11 = 'datasets_small/u11_0088_peaklist.h5'
peaks_u101 = 'datasets_small/u101_0080_peaklist.h5'
peaks_u102 = 'datasets_small/u102_0080_peaklist.h5'
peaks_u103 = 'datasets_small/u103_0080_peaklist.h5'
peaks_u14  = 'datasets_small/u14_0025_peaklist.h5'
peak_list = {'u05':peaks_u05,'u10':peaks_u10,'u11':peaks_u11, 'u101':peaks_u101,'u102':peaks_u102,'u103':peaks_u103}
peak_list['u14']=peaks_u14


target_frames={'u05':125,'u10':82,'u11':88,'u101':80,'u102':80,'u103':80,'u14':25}

every_ten = {'u05':u05_every_ten,'u10':u10_every_ten,'u11':u11_every_ten, 'u101':u101_every_ten,'u102':u102_every_ten,'u103':u103_every_ten}
every_ten.update({'u201':u201_every_ten,'u202':u202_every_ten,'u203':u203_every_ten})



bad_particles_u05='datasets_small/u05_bad_particles.h5'
bad_particles={'u05':bad_particles_u05,'u10':None,'u11':None}

n_particles={'u05':'datasets_small/u05_n_particles.txt',
             'u10':'datasets_small/u10_n_particles.txt',
             'u11':'datasets_small/u11_n_particles.txt',
             'u101':'datasets_small/u101_n_particles.txt',
             'u102':'datasets_small/u102_n_particles.txt',
             'u103':'datasets_small/u103_n_particles.txt',
             'u14':'datasets_small/u14_n_particles.txt'}

