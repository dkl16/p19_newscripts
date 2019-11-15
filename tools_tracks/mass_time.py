from starter2 import *
import data_locations as dl
import davetools
reload(davetools)

plt.close('all')


file_list=glob.glob('%s/*h5'%dl.sixteen_frame)
output_prefix=''
file_list=glob.glob('frame_2_20_30*h5')[:5]
output_prefix='f_2_20_30'
#for debug purposes you may want a reduced list 
#file_list=file_list[:3]    

if 'this_looper' not in dir():
    this_looper=looper.core_looper(directory=dl.enzo_directory)
    for nfile,fname in enumerate(file_list):
        this_looper.load_loop(fname)
        print( "Reading file %d of %d"%(nfile,len(file_list)))
    thtr = this_looper.tr
    thtr.sort_time()
    all_cores = np.unique(thtr.core_ids)

core_list=all_cores[:1]
rm = rainbow_map(len(all_cores))

if 'rho_extents' not in dir():
    rho_extents=davetools.extents()
    r_extents=davetools.extents()
    for nc,core_id in enumerate(all_cores):
        ms = trackage.mini_scrubber(thtr,core_id)
        if ms.nparticles == 1:
            continue
        density = thtr.c([core_id],'density')
        rho_extents(density)
        r_extents(ms.r)

dx = 1./2048
nx = 1./dx
for nc,core_id in enumerate(core_list):

    for nf,frame in enumerate(thtr.frames):
        x =np.floor(thtr.c([core_id],'x')/dx)[:,nf]#or whatever the number of zones is
        y =np.floor(thtr.c([core_id],'y')/dx)[:,nf]
        z =np.floor(thtr.c([core_id],'z')/dx)[:,nf]
        density = thtr.c([core_id],'density')[:,nf]
        cell_volume = thtr.c([core_id],'cell_volume')[:,nf]
        index = x + nx*(y * nx*z)
        ar = np.argsort(index)
        rs = np.argsort(ar)
        isorted=index[ar]
        mask = np.ones_like(density,dtype='bool')
        mask[1:] = isorted[1:]-isorted[:-1] != 0
        mask2 = mask[ rs]
        mass = (density[mask2]*cell_volume[mask2]).sum()
