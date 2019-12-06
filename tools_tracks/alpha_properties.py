from starter2 import *
import data_locations as dl
import davetools
reload(davetools)
import scatter_fit
reload(scatter_fit)

plt.close('all')


file_list=glob.glob('%s/*h5'%dl.sixteen_frame)
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

core_list=all_cores
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

asort =  np.argsort(thtr.times)
if (asort != sorted(asort)).any():
    print("Warning: times not sorted.")
do_all_plots=False
if 'alpha' not in dir():
    alpha = []
    mass = []
    density_0=[]
    mean_div = []
    alpha_dict={}

    for nc,core_id in enumerate(core_list):

        #miniscrubber computes distance, r^2, several other quantities
        ms = trackage.mini_scrubber(thtr,core_id)
        if ms.nparticles == 1:
            continue

        density = thtr.c([core_id],'density')


        #Compute alpha
        #compute mean initial density.
        this_r = ms.r.flatten()
        min_r = 1./2048
        this_r[ this_r < min_r ] = min_r    
        fit = scatter_fit.scatter_fit(plt,this_r,density.flatten())
        alpha.append(fit['fit'][0])
        alpha_dict[core_id]=alpha[-1]
        density_0.append( density[:,0].mean()) #this is not good

        #compute initial mass.
        for nf in [0]:
            dx=1./2048
            nx=2048
            x =np.floor(thtr.c([core_id],'x')/dx)[:,nf]
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
            mass.append((density[mask2]*cell_volume[mask2]).sum())

        if do_all_plots:
            tmap=rainbow_map(ms.ntimes)
            fig, axd1=plt.subplots(1,1)
            for n_count,n_time in enumerate(asort):
                time=thtr.times[n_time]
                if time == 0:
                    continue
                c=tmap(n_count,ms.nparticles)
                this_r=ms.r[:,n_time]+0

                axd1.scatter(this_r,density[:,n_time],c=c,label=thtr.times[n_time],s=0.1)
            axd1.set_title(r'$\alpha = %0.3f$'%alpha[-1])
            davetools.axbonk(axd1,xscale='log',yscale='log',xlabel='r',ylabel=r'$\rho$',
                             xlim=r_extents.minmax, ylim=rho_extents.minmax)
            outname = '%s/density_radius_c%04d'%(dl.output_directory,core_id)
            fig.savefig(outname)
            print("saved "+outname)
            plt.close(fig)

if 'alpha_proxy' not in dir() or True:
    do_all_plots=True
    div_v=[]
    v_over_r=[]
    divv_extents=davetools.extents()
    divv_extents( thtr.track_dict['velocity_divergence'])
    proxy_extents=davetools.extents()
    for nc,core_id in enumerate(core_list):
        #miniscrubber computes distance, r^2, several other quantities
        ms = trackage.mini_scrubber(thtr,core_id)
        if ms.nparticles == 1:
            continue
        density = thtr.c([core_id],'density')
        div_v = thtr.c([core_id],'velocity_divergence')
        mag_v = thtr.c([core_id],'velocity_magnitude')
        divv_extents(div_v)
        if do_all_plots:
            tmap=rainbow_map(ms.ntimes)
            fig, axes=plt.subplots(2,2)
            ax00=axes[0][0]; ax01=axes[0][1]
            ax10=axes[1][0]; ax11=axes[1][1]

            for n_count,n_time in enumerate(asort):
                time=thtr.times[n_time]
                if time == 0:
                    continue
                c=tmap(n_count,ms.nparticles)
                this_r=ms.r[:,n_time]+0

                ax00.scatter(this_r,div_v[:,n_time],c=c,label=thtr.times[n_time],s=0.1)
                ax01.scatter(this_r,mag_v[:,n_time],c=c,label=thtr.times[n_time],s=0.1)


                alpha_proxy = div_v[:,n_time]/mag_v[:,n_time]*ms.r[:,n_time]
                proxy_extents(alpha_proxy)
                ax10.scatter(this_r,alpha_proxy.flatten(),c=c,label=thtr.times[n_time],s=0.1)
            davetools.axbonk(ax00,xscale='log',yscale='linear',xlabel='r',ylabel=r'$\nabla\cdot v$',
                             xlim=r_extents.minmax, ylim=rho_extents.minmax)
            ax00.set_yscale('symlog',linthreshy=1)
            ax00.set_ylim(-5e5,5e5)

            davetools.axbonk(ax01,xscale='log',yscale='log',xlabel='r',ylabel=r'$|v|$',
                             xlim=r_extents.minmax, ylim=rho_extents.minmax)
            ax01.set_yscale('symlog',linthreshy=1)
            ax01.set_ylim(0,100)

            ax10.plot([1e-3,1e-1], [alpha_dict[core_id],alpha_dict[core_id]],c='k')
            ax10.plot([1e-3,1e-1], [-0.5,-0.5],c=[0.5]*4)
            ax10.plot([1e-3,1e-1], [-3,-3],c=[0.5]*4)
            davetools.axbonk(ax10,xscale='log',yscale='linear',xlabel='r',ylabel=r'$\nabla\cdot v/(v/r)$',
                             xlim=r_extents.minmax, ylim=rho_extents.minmax)
            ax10.set_yscale('symlog',linthreshy=0.1)
            ax10.set_ylim(proxy_extents.minmax)

            outname = '%s/divergence_proxy_c%04d'%(dl.output_directory,core_id)
            fig.savefig(outname)
            print("saved "+outname)
            plt.close(fig)


if 0:
    mass=nar(mass)
    fig, axes=plt.subplots(1,2)
    ax0 = axes[0]; ax1=axes[1]
    ax0.scatter(density_0,alpha)
    ax1.scatter(mass,alpha)
    davetools.axbonk(ax0,ylabel=r'$\alpha$', xlabel=r'$\langle \rho \rangle$', xscale='log',yscale='linear')
    davetools.axbonk(ax1,ylabel=r'$\alpha$', xlabel=r'$M(t=0)$', xscale='log',yscale='linear')
    ax1.set_xlim(mass.min(),mass.max())
    outname = '%s/alpha_mass.png'%dl.output_directory
    fig.savefig(outname)
    print(outname)


