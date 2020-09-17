from starter2 import *
import matplotlib.image as mpimg

import data_locations as dl
reload(dl)
plt.close('all')

file_list=glob.glob('%s/*h5'%dl.sixteen_frame)
out_prefix = 'CHANGE_THIS_PREFIX'

out_prefix='tracks'
file_list = file_list[:2]


if 'this_looper' not in dir():
    this_looper=looper.core_looper(directory=dl.enzo_directory)
    for nfile,fname in enumerate(file_list):
        this_looper.load_loop(fname)
        print( "File %d of %d"%(nfile,len(file_list)))
    thtr = this_looper.tr
    thtr.sort_time()
all_cores = np.unique(thtr.core_ids)
rm = rainbow_map(len(all_cores))

tsorted = thtr.times
core_list = all_cores
if 'collapse_times' not in dir() or True:
    collapse_times=[]
    used_cores=[]
    for core_id in core_list:
        plt.close('all')
        used_cores.append(core_id)
        n0=0

        fig,ax=plt.subplots(1,1)
        ms = trackage.mini_scrubber(thtr,core_id)
        density = thtr.c([core_id],'density')
        cell_volume = thtr.c([core_id],'cell_volume')

        tmap=rainbow_map(ms.ntimes)
        norm = mpl.colors.Normalize()
        norm.autoscale( np.log10(density[:,n0]))
        cmap = mpl.cm.jet
        color_map = mpl.cm.ScalarMappable(norm=norm,cmap=cmap)


        for npart in list(range(ms.nparticles))[:10]:
            c = color_map.to_rgba(np.log10(density[npart,n0]))
            ax.plot( thtr.times, density[npart,:],c=c,linewidth=.1)#linestyle=':')
        ax.plot(tsorted, density.mean(axis=0)[:],c='k')

        #collapse time
        if (density>1e3).any():
            first_collapse_index = np.where( (density>1e3).sum(axis=0) >0 )[0][0] 
            t_collapse = thtr.times[first_collapse_index]
            rho_col = density[:,first_collapse_index].max()
            #ax.scatter( t_collapse, rho_col ,marker='*',s=5)
            ax.plot( [t_collapse-0.005, t_collapse+0.005],[rho_col,rho_col],c=[0.5]*3)
            ax.plot( [t_collapse     , t_collapse     ],[rho_col/5,5*rho_col],c=[0.5]*3)
            collapse_times.append(t_collapse)

        t0 = thtr.times[0]
        t1 = thtr.times[-1]
        rho0 = (density[:,0]*cell_volume[:,0]).sum()/cell_volume[:,0].sum()
        rho1 = density[:,0].max() 
        alpha = 1.8
        G=1620./(4*np.pi)
        tff_global = np.sqrt(3*np.pi/(32*G*1))

        tff_local = np.sqrt(3*np.pi/(32*G*rho0))
        rhot = rho0*(1-(tsorted/tff_local)**2)**-alpha

        tff_max = np.sqrt(3*np.pi/(32*G*rho1))
        rho_max = rho1*(1-(tsorted/tff_max)**2)**-alpha
        #tc =t1*(1-(rho1/rho0)**(-1./alpha))**-0.5
        #rho_c = 3*np.pi/(32*G*tc**2)

        ok = np.isnan(rhot)==False
        okmax = np.isnan(rho_max)==False
        ax.plot( tsorted[ok], rhot[ok], c='r',label=r'$tc/tff = %0.2e$'%(tc/tff_local))
        ax.plot( tsorted[okmax], rho_max[okmax], c='g',label=r'$tc/tff = %0.2e$'%(tc/tff_local))
        ax.plot( [tff_max,tff_max],[0.1,1e5],c='r')
        ax.plot( [tff_local,tff_local],[0.1,1e5],c='r')
        ax.plot( [tff_global,tff_global],[0.1,1e5],c='r')

        ax.legend(loc=0)
        ylim = [thtr.track_dict['density'].min(), thtr.track_dict['density'].max()]
        xlim = [0, 0.05]

        axbonk(ax,xscale='linear',yscale='log',xlabel='t',ylabel=r'$\rho$', ylim=ylim,xlim=xlim)
        oname = "%s/%s_density_6_c%04d"%(dl.output_directory,out_prefix,core_id)
        fig.savefig(oname)
        print("Saved "+oname)

if 1:
    fig, ax = plt.subplots(1,1) 
    ax.hist( collapse_times, histtype='step')
    axbonk(ax, xlabel=r'$t_c$', ylabel = 'N')
    oname = "%s/%s_tc_hist"%(dl.output_directory,out_prefix)
    fig.savefig( oname)
    print(oname)

