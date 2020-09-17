"""
This just pulls partilce data from enzo and stores it.
Changing
core_list 
frame_list
fields
changes what gets extracted.
"""
from starter2 import *
import xtra_energy
import data_locations as dl
reload(dl)
#
# set sim
#
this_simname = 'u10'
all_nonzero = looper.get_all_nonzero(dl.n_particles[this_simname])

output_base = "%s_cores"%this_simname
if 1:
    """this set of parameters extracts all primitive quantities"""
    core_list = all_nonzero
    frame_list = [0]#+list(range(10,130,10))+[125]
    fields = ['x','y','z','density','velocity_magnitude','magnetic_field_strength', 'velocity_divergence']
    fields += ['velocity_x','velocity_y','velocity_z']
    fields += ['magnetic_field_%s'%s for s in 'xyz']
    output_base = "primitive_test"
    derived=[]

if 0:
    """this set of parameters extracts all primitive quantities"""
    core_list = all_nonzero.astype('int')
    frame_list = [0,1]+list(range(10,130,10))+[125]
    fields = ['x','y','z','velocity_magnitude','magnetic_field_strength', 'velocity_divergence']
    fields += ['velocity_x','velocity_y','velocity_z']
    fields += ['magnetic_field_%s'%s for s in 'xyz']
    output_base = "all_primitives"
    derived=[]

if 0:
    """This set extracts magnetic work"""
    core_list = all_nonzero.astype('int')
    frame_list =[0,1]+list(range(10,130,50))+[125]
    fields = ['mag_work']
    derived=[xtra_energy.add_force_terms]
    output_base = 'mag_work_only'

if 0:
    #Pull a whole list for each core, saving each core to its own file
    for core in core_list:  
        output_name = '%s_c%04d.h5'%(output_base,core)
        if os.path.exists(output_name):
            print("File exists, skipping "+output_name)
            continue
        this_looper = looper.core_looper(directory= dl.sims[this_simname],
                                         sim_name = this_simname,
                                         out_prefix = this_simname,
                                         target_frame = dl.target_frames[this_simname],
                                         frame_list = frame_list,
                                         core_list =  [core],# core_list,
                                         fields_from_grid=fields,
                                         derived = derived
                                      )
        this_looper.get_target_indices(h5_name=dl.peak_list[this_simname],
                                         bad_particle_list=dl.bad_particles[this_simname])
        this_looper.get_tracks()
        this_looper.save(output_name)

if 1:
    #Pull all cores at once. 
    output_name = '%s_many.h5'%(output_base)
    if not os.path.exists(output_name):
        this_looper = looper.core_looper(directory= dl.sims[this_simname],
                                         sim_name = this_simname,
                                         out_prefix = this_simname,
                                         target_frame = dl.target_frames[this_simname],
                                         frame_list = frame_list,
                                         core_list =  core_list,
                                         fields_from_grid=fields,
                                         derived = derived
                                      )
        this_looper.get_target_indices(h5_name=dl.peak_list[this_simname],
                                         bad_particle_list=dl.bad_particles[this_simname])
        this_looper.get_tracks()
        #this_looper.save(output_name)
