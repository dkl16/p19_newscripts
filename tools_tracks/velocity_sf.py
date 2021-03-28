"""

THIS ALMOST WORKS.
I need to tweak the structure function in tools_tracks/sf2 in some
odd ways to make it work, though.  There's a wierd extra 5.

"""




from starter2 import *
import matplotlib.image as mpimg

from scipy.optimize import curve_fit
import data_locations as dl
reload(dl)
plt.close('all')


class trial():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.cores_used=[]
        self.r=[]
        self.v=[]



    def by_frame(self,frame,core_list=None):
        self.r=[]
        self.v=[]

        thtr = self.this_looper.tr
        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list = all_cores

        thtr.sort_time()

        tsorted = thtr.times
        self.core_list=core_list
        fig,ax=plt.subplots(1,1)
        rmin, rmax = 1./2048, 0.4
        vmin, vmax = 0, 8
        nx=ny=32
        self.rbins = np.logspace(np.log10(rmin), np.log10(vmax),nx+1)
        self.vbins = np.linspace(0, vmax, ny+1)
        self.hist = np.zeros([nx,ny])
        def cen(arr):
            return 0.5*(arr[1:]+arr[:-1])
        self.TheX = np.r_[(ny)*[cen(self.rbins)]].transpose()
        self.TheY = np.r_[(nx)*[cen(self.vbins)]]
        self.x_del = (self.rbins[1:]-self.rbins[:-1])
        self.y_del = (self.vbins[1:]-self.vbins[:-1])
        self.x_del.shape = (self.x_del.size,1)
        self.dv = 1./(self.x_del*self.y_del)

        for core_id in core_list:
            ms = trackage.mini_scrubber(thtr,core_id)
            tmap=rainbow_map(ms.ntimes)
            if ms.nparticles < 15:
                continue
            print('go ', core_id)
            self.cores_used.append(core_id)
            n0=0
            rmin = 1./2048
            nt = np.where(thtr.frames == frame)[0][0]

            density = thtr.c([core_id],'density')[:,nt]
            cell_volume = thtr.c([core_id],'cell_volume')[:,nt]
            this_r = ms.r[:,nt]
            this_r[ this_r < rmin] = rmin
            asort = np.argsort(this_r)
            unsort = np.argsort(asort)
            rsort = this_r[asort]
            vr = ms.vr_rel[:,nt]
            vrs = vr[asort]
            dv = cell_volume[asort]
            vmean = np.sqrt(np.cumsum(vrs**2*dv)/np.cumsum(dv))
            self.r.append(rsort)
            self.v.append(vmean)
            ax.plot(rsort,vmean, c=tmap(nt))

            this_hist, xedge, yedge= np.histogram2d(rsort, vmean, bins=[self.rbins,self.vbins])
            self.hist+= this_hist.astype('float')
            #rquan = np.digitize(rsort, rbins)
            #vquan = np.digitize(vmean, vbins)
            #self.hist[rquan,vquan] += 1
        

        axbonk(ax,xlabel='r',ylabel=r'$\sigma_v(r)$',xscale='log',xlim=[rmin,0.4], ylim=[0,10])
        #ax.set_yscale('symlog',linthresh=10)

        outname = "%s/%s_vmeans_cXXXX_n%04d.png"%(dl.output_directory, self.this_looper.out_prefix, frame)
        print(outname)
        fig.savefig(outname)
        fig,ax=plt.subplots(1,1)

        ax.pcolormesh(self.rbins, self.vbins, self.hist)
        fig.savefig("%s/%s_hist_cXXXX_n%04d.png"%(dl.output_directory, self.this_looper.out_prefix, frame))




        plt.close('all')
    def plot(self,frame):
        pdf = self.hist/(self.dv)
        pdf /= (pdf*self.dv).sum()
        fig,ax=plt.subplots(1,1)
        cmap = copy.copy(mpl.cm.get_cmap("viridis"))
        cmap.set_under('w')
        minmin = pdf[pdf>0].min()
        norm = mpl.colors.LogNorm(vmin=minmin,vmax=pdf.max())
        for r,v in zip(self.r,self.v):
            ax.plot(r,v,c=[0.5,0.5,0.5,0.3],lw=0.1)
        for r,v in zip(self.r,self.v):
            ax.scatter(r[0],v[0],c='k')
        ax.plot( self.rbins, [1.0]*self.rbins.size,c=[0.5]*3)
        ploot=ax.pcolormesh(self.TheX, self.TheY, pdf,cmap=cmap,norm=norm,alpha=0.5)
        axbonk(ax,yscale='linear',xscale='log', xlim=[1./2048,0.4], ylim=[0,9])
        fig.colorbar(ploot,ax=ax)
        DO_SF = False
        if DO_SF:
            import sf2
            reload( sf2)
            if 'msf' not in dir() or True:
                msf = sf2.make_sf(TL.looper1,0)
                msf.bin_kludged()
            ax.plot( msf.binned2[1], msf.binned2[2],c='k')
        ax.set_title('TAKE 2')
        fig.savefig("%s/%s_hist_cXXXX_n%04d.png"%(dl.output_directory, self.this_looper.out_prefix, frame))
        plt.close(fig)


    def run(self,do_all_plots=True,core_list=None):

        thtr = self.this_looper.tr
        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list = all_cores

        thtr.sort_time()

        tsorted = thtr.times
        self.core_list=core_list
        fig,ax=plt.subplots(1,1)
        for core_id in core_list:
            ms = trackage.mini_scrubber(thtr,core_id)
            tmap=rainbow_map(ms.ntimes)
            if ms.nparticles < 3:
                continue
            print('go ', core_id)
            self.cores_used.append(core_id)
            n0=0
            ax.clear()
            rmin = 1./2048
            for nt, time in enumerate(tsorted):
                density = thtr.c([core_id],'density')[:,nt]
                cell_volume = thtr.c([core_id],'cell_volume')[:,nt]

                this_r = ms.r[:,nt]
                this_r[ this_r < rmin] = rmin
                asort = np.argsort(this_r)
                unsort = np.argsort(asort)
                rsort = this_r[asort]
                vr = ms.vr_rel[:,nt]
                vrs = vr[asort]
                dv = cell_volume[asort]
                vmean = np.sqrt(np.cumsum(vrs**2*dv)/np.cumsum(dv))
                ax.plot(rsort,vmean, c=tmap(nt))
            axbonk(ax,xlabel='r',ylabel=r'$\sigma_v(r)$',xscale='log',xlim=[rmin,0.4], ylim=[0,10])
            #ax.set_yscale('symlog',linthresh=10)

            outname = "%s/%s_vmeans_c%04d.png"%(dl.output_directory, self.this_looper.out_prefix, core_id)
            print(outname)
            fig.savefig(outname)


        plt.close(fig)





if 'do_all_plots' not in dir():
    do_all_plots = False


import three_loopers as TL

if 1:
    run1 = trial(TL.looper1)
#run1.run()#core_list=[10,11,31])

for frame in range(100):
    run1.by_frame(frame)
    run1.plot(frame)


