from starter2 import *
import matplotlib.image as mpimg

from scipy.optimize import curve_fit
import data_locations as dl
reload(dl)
plt.close('all')

if 'this_simname' not in dir():
    this_simname = 'u11'
file_list=glob.glob(dl.every_ten[this_simname])

out_prefix="%s_"%this_simname
#file_list = file_list[:2]

if 'this_looper' not in dir():
    this_looper=looper.core_looper(directory=dl.sims[this_simname])
    for nfile,fname in enumerate(file_list):
        this_looper.load_loop(fname)
        print( "File %d of %d"%(nfile,len(file_list)))
    thtr = this_looper.tr
    thtr.sort_time()
all_cores = np.unique(thtr.core_ids)
rm = rainbow_map(len(all_cores))

#if 'collapse_times' not in dir():
class trial():
    def __init__(self):
        self.collapse_times=[]
        self.used_cores=[]
        self.tff_local_list=[]
        self.tff_harm_list=[]
        self.rho_harmonic = []
        self.rho_0_list = []
        self.rho_c_list = []
        self.tc=[]
    def count_particles(self,the_looper):
        cores=[]
        nparticles=[]
        for core in the_looper.target_indices:
            cores.append(core)
            nparticles.append(the_looper.target_indices[core].size)
        return {'core':nar(core),'nparticles':nar(nparticles)}
    def run(self,this_looper,do_all_plots=True,core_list=None):
        collapse_times=[]

        thtr = this_looper.tr
        thtr.sort_time()
        all_cores = np.unique(thtr.core_ids)
        rm = rainbow_map(len(all_cores))

        tsorted = thtr.times
        if core_list is None:
            core_list = all_cores
        for core_id in core_list:
            print('go ', core_id)
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

            #collapse time
            if (density>1e3).any():
                first_collapse_index = np.where( (density>1e3).sum(axis=0) >0 )[0][0] 
                t_collapse = thtr.times[first_collapse_index]
                rho_col = density[:,first_collapse_index].max()
                #ax.scatter( t_collapse, rho_col ,marker='*',s=5)
                ax.plot( [t_collapse-0.005, t_collapse+0.005],[rho_col,rho_col],c=[0.5]*3)
                ax.plot( [t_collapse     , t_collapse     ],[rho_col/5,5*rho_col],c=[0.5]*3)
                collapse_times.append(t_collapse)
            else:
                rho_col = -1
                t_collapse = -1
                collapse_times.append(-1)

            self.collapse_times=collapse_times
            t0 = thtr.times[0]
            t1 = thtr.times[-1]
            rho0 = (density[:,0]*cell_volume[:,0]).sum()/cell_volume[:,0].sum()
            rho1 = density[:,0].max() 
            alpha = 1.8
            G=1620./(4*np.pi)
            tff_global = np.sqrt(3*np.pi/(32*G*1))

            tff_local = np.sqrt(3*np.pi/(32*G*rho0))
            self.tff_local_list.append(tff_local)

            rhot = rho0*(1-(tsorted/tff_local)**2)**-alpha

            tff_max = np.sqrt(3*np.pi/(32*G*rho1))
            rho_max = rho1*(1-(tsorted/tff_max)**2)**-alpha

            #this was not particularly useful
            rho_harm =  np.exp((np.log(density[:,0])*cell_volume[:,0]).sum()/cell_volume[:,0].sum()) 
            self.rho_harmonic.append(rho_harm)
            self.rho_0_list.append(rho0)
            self.tff_harm_list.append( np.sqrt(3*np.pi/(32*G*rho_harm)))



            if do_all_plots:
                if ms.nparticles<100:
                    this_slice=slice(None)
                else:
                    this_slice=slice(None,None,10)

                for npart in list(range(ms.nparticles))[this_slice]:
                    c = color_map.to_rgba(np.log10(density[npart,n0]))
                    ax.plot( thtr.times, density[npart,:],c=c,linewidth=.1)#linestyle=':')
                rho_mean = density.mean(axis=0)[:]
                ax.plot(tsorted, rho_mean,c='k')
                #ok = np.isnan(rhot)==False
                #okmax = np.isnan(rho_max)==False
                #ax.plot( tsorted[ok], rhot[ok], c='r',label=r'$t_{\rm{ff}}$')
                #ax.plot( tsorted[okmax], rho_max[okmax], c='g',label=r'$t_{\rm{ff,max}}$')
                #ax.plot( [tff_max,tff_max],[0.1,1e5],c='r')
                #ax.plot( [tff_local,tff_local],[0.1,1e5],c='r')
                #ax.plot( [tff_global,tff_global],[0.1,1e5],c='r')

                #tc =t1*(1-(rho1/rho0)**(-1./alpha))**-0.5
                #rho_c = 3*np.pi/(32*G*tc**2)
                #rhot = rho0*(1-(tsorted/tff_local)**2)**-alpha

                if rho_col > 0:
                    tc =t_collapse*(1-(rho_col/rho0)**(-1./alpha))**-0.5
                    rho_c = 3*np.pi/(32*G*tc**2)
                    rho_tc = rho_c*(1-(tsorted/tc)**2)**-alpha
                    ax.plot( tsorted, rho_tc, c='g')
                    print("stuff rhoc/rho0 %0.2f"%(rho_c/rho0))
                    rho_tff = rho0*(1-(tsorted/tff_local)**2)**-alpha
                    ax.plot( tsorted, rho_tff, c='b')
                    #rho_c = 3*np.pi/(32*G*tc**2)
                    self.tc.append(tc)
                    self.rho_c_list.append(rho_c)
                else:
                    self.tc.append(-1)
                    self.rho_c_list.append(-1)
                def rho_ff_2(t,rho0):
                    alpha = 1.8
                    rho_max = rho0*(1-(t/np.sqrt(3*np.pi/(32*G*rho0)))**2)**-alpha
                    return rho_max

                def rho_ff(t,rho0):
                    alpha = 1.8
                    rho_max = rho0*(1-(t/np.sqrt(3*np.pi/(32*G*rho0)))**2)**-alpha
                    return rho_max
                fits, cov = curve_fit(rho_ff, thtr.times, np.log(rho_mean), [rho0])
                rho_ff_fit = rho_ff(tsorted, fits[0])
                ok = np.where(np.isnan(rho_ff_fit)==False)
                #ax.plot( tsorted[ok], rho_ff_fit[ok] ,c='r')
                self.fits=fits


                ax.legend(loc=0)
                ylim = [thtr.track_dict['density'].min(), thtr.track_dict['density'].max()]
                xlim = [0, 0.05]

                axbonk(ax,xscale='linear',yscale='log',xlabel='t',ylabel=r'$\rho$', ylim=ylim,xlim=xlim)
                oname = "%s/%s_density_6_c%04d"%(dl.output_directory,out_prefix,core_id)
                fig.savefig(oname)
                print("Saved "+oname)

if 'other_looper' not in dir():
    other_looper=looper.core_looper(directory=dl.sims['u05'])
    file_list=glob.glob(dl.every_ten['u05'])
    for nfile,fname in enumerate(file_list):
        other_looper.load_loop(fname)
        print( "File %d of %d"%(nfile,len(file_list)))

if 'third_looper' not in dir():
    third_looper=looper.core_looper(directory=dl.sims['u11'])
    file_list=glob.glob(dl.every_ten['u11'])
    for nfile,fname in enumerate(file_list):
        third_looper.load_loop(fname)
        print( "File %d of %d"%(nfile,len(file_list)))


tool1=trial()
tool1.run(this_looper,do_all_plots=True)#,core_list=[328,181,305])

if 'tool1' not in dir():
    tool1=trial()
    tool1.run(this_looper,do_all_plots=False)
    tool2=trial()
    tool2.run(other_looper,do_all_plots=False)
    tool3=trial()
    tool3.run(third_looper,do_all_plots=False)

stuff1 = tool1.count_particles( this_looper)
stuff2 = tool2.count_particles( other_looper)
stuff3 = tool3.count_particles( third_looper)

if 0:
    G=1620./(4*np.pi)
    tff_global = np.sqrt(3*np.pi/(32*G*1))
    fig, ax = plt.subplots(1,1) 
    ct1=nar(tool1.collapse_times)/tff_global
    ct2=nar(tool2.collapse_times)/tff_global
    ct3=nar(tool3.collapse_times)/tff_global
    ok1=ct1>0
    ok2=ct2>0
    ok3=ct3>0
    ax.hist( ct1[ok1], histtype='step',label='u05')
    ax.hist( ct2[ok2], histtype='step',label='u10')
    ax.hist( ct3[ok3], histtype='step',label='u11')
    axbonk(ax, xlabel=r'$t_c/t_{\rm{ff}}$', ylabel = 'N')
    oname = "%s/%s_tc_hist"%(dl.output_directory,out_prefix)
    fig.savefig( oname)
    plt.close(fig)
    print(oname)

if 0:
    fig, ax = plt.subplots(1,1) 
    ax.scatter(stuff1['nparticles'][ok1],ct1[ok1],label='u05')
    ax.scatter(stuff2['nparticles'][ok2],ct2[ok2],label='u10')
    ax.scatter(stuff3['nparticles'][ok3],ct3[ok3],label='u11')
    axbonk(ax,xlabel='N',ylabel=r'$t_c/t_{\rm{ff}}$',xscale='log')
    ax.legend(loc=0)
    oname = "%s/%s_np_tc_hist"%(dl.output_directory,out_prefix)
    fig.savefig( oname)
    plt.close(fig)
    print(oname)

if 0:

    fig, ax = plt.subplots(1,1) 

    ax.scatter(stuff1['nparticles'][ok1],ct1[ok1]*tff_global/nar(tool1.tff_local_list)[ok1],label='u05')
    ax.scatter(stuff2['nparticles'][ok2],ct2[ok2]*tff_global/nar(tool2.tff_local_list)[ok2],label='u10')
    ax.scatter(stuff3['nparticles'][ok3],ct3[ok3]*tff_global/nar(tool3.tff_local_list)[ok3],label='u11')
    axbonk(ax,xlabel='N',ylabel=r'$t_c/t_{\rm{ff,local}}$',xscale='log')
    ax.legend(loc=0)
    oname = "%s/%s_np_tc_local"%(dl.output_directory,out_prefix)
    fig.savefig( oname)
    plt.close(fig)
    print(oname)

if 0:

    fig, ax = plt.subplots(1,1) 

    ax.scatter(stuff1['nparticles'][ok1],ct1[ok1]*tff_global/nar(tool1.tff_harm_list)[ok1],label='u05')
    ax.scatter(stuff2['nparticles'][ok2],ct2[ok2]*tff_global/nar(tool2.tff_harm_list)[ok2],label='u10')
    ax.scatter(stuff3['nparticles'][ok3],ct3[ok3]*tff_global/nar(tool3.tff_harm_list)[ok3],label='u11')
    axbonk(ax,xlabel='N',ylabel=r'$t_c/t_{\rm{ff,harmonic}}$',xscale='log')
    ax.legend(loc=0)
    oname = "%s/%s_np_tc_harmonic"%(dl.output_directory,out_prefix)
    fig.savefig( oname)
    plt.close(fig)
    print(oname)

if 0:

    fig, ax = plt.subplots(1,1) 

    rho_a =min(tool1.rho_0_list)
    rho_b =max(tool1.rho_0_list)
    ax.plot( [rho_a,rho_b],[rho_a,rho_b], c='k')
    ax.plot( [2*rho_a,2*rho_b],[rho_a,rho_b], c='k')
    ax.scatter(tool1.rho_0_list, tool1.rho_harmonic,label='u05',marker='.',s=.9,c='r')
    ax.scatter(tool2.rho_0_list, tool2.rho_harmonic,label='u10',marker='.',s=.9,c='g')
    ax.scatter(tool3.rho_0_list, tool3.rho_harmonic,label='u11',marker='.',s=.9,c='b')
    axbonk(ax,xlabel=r'$\rho_0$',ylabel=r'$\rho_{\rm{harm}}$',xscale='log',yscale='log')
    ax.legend(loc=0)
    oname = "%s/%s_rho_harmonic"%(dl.output_directory,out_prefix)
    fig.savefig( oname)
    plt.close(fig)
    print(oname)
