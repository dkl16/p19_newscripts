
from starter2 import *
import matplotlib.image as mpimg
from scipy.spatial import ConvexHull
import matplotlib.patches as patches
import data_locations as dl
reload(dl)
plt.close('all')


class hull_tool():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.cell_volumes=[]
        self.hull_volumes=[]
        self.good_cores=[]

    def run(self,do_2d_plots=False,core_list=None):

        thtr = self.this_looper.tr
        all_cores = np.unique(thtr.core_ids)
        rm = rainbow_map(len(all_cores))
        if core_list is None:
            core_list = all_cores
        if do_2d_plots:
            fig_many, ax_many = plt.subplots(2,2,figsize=(8,8))
            ax_all = [ax_many[0][0], ax_many[0][1], ax_many[1][0]]
        for core_id in core_list:
            ms = trackage.mini_scrubber(thtr,core_id)
            if ms.r.shape[0] <= 4:
                continue
            self.good_cores.append(core_id)
            asort =  np.argsort(thtr.times)
            delta=0.1

            mask = slice(None)
            for it,nt in enumerate([0]):#asort):

                this_x,this_y,this_z=ms.this_x[mask,nt],ms.this_y[mask,nt], ms.this_z[mask,nt]
                this_p = [this_x,this_y,this_z]

                points_3d = np.array(list(zip(this_x,this_y,this_z)))
                hull_3d = ConvexHull(points_3d)
                self.hull_volumes.append(hull_3d.volume)
                self.cell_volumes.append( thtr.c([core_id],'cell_volume')[mask,nt].sum())


                if do_2d_plots:
                    #for ax in ax_all:
                    #    ax.clear()
                    #    ax.set_aspect('equal')
                    #    ax.plot([0,1,1,0,0],[0,0,1,1,0])

                    for LOS in [0,1,2]:
                        x = [0,2,0][LOS]
                        y = [1,1,2][LOS]
                        xlab='xyz'[x]
                        ylab='yzx'[y]



                        ax_all[LOS].scatter(this_p[x], this_p[y],s=0.1)

                        points_2d = np.array(list(zip(this_p[x],this_p[y])))
                        hull_2d = ConvexHull(points_2d)
                        vert_x = points_2d[hull_2d.vertices,0]
                        vert_y = points_2d[hull_2d.vertices,1]
                        vert_x = np.concatenate([vert_x,vert_x[0:1]])
                        vert_y = np.concatenate([vert_y,vert_y[0:1]])
                        ax_all[LOS].plot(vert_x, vert_y, 'k')

                        x_min = min([this_p[x].min(), -delta])
                        x_max = max([this_p[x].max(), 1+delta])
                        y_min = min([this_p[y].min(), -delta])
                        y_max = max([this_p[y].max(), 1+delta])


                        axbonk(ax_all[LOS],xlabel=xlab,ylabel=ylab,xlim=[x_min,x_max],ylim=[y_min,y_max])
                        #ax_many.set_title(title)
                    outname = '%s/%s_hull_3d_t_c%04d_n%04d.png'%(dl.output_directory,self.this_looper.out_prefix,core_id,it)
                    fig_many.savefig(outname)
                    print("Wrote "+outname)

import three_loopers as tl
if 'ht1' not in dir() or clobber:
    ht1 = hull_tool(tl.looper1)
    ht1.run(do_2d_plots=True)#,core_list=[10,11])
if 'ht2' not in dir() or clobber:
    ht2 = hull_tool(tl.looper2)
    ht2.run(do_2d_plots=True)
if 'ht3' not in dir() or clobber:
    ht3 = hull_tool(tl.looper3)
    ht3.run(do_2d_plots=True)

if 1:
    import p56_plots.density_AC as AC
    reload(AC)
    if 'a1' not in dir():
        a1 = AC.ac_thing('u05'); a1.plot()
        a2 = AC.ac_thing('u10'); a2.plot()
        a3 = AC.ac_thing('u11'); a3.plot()
        acs={'u05':a1,'u10':a2,'u11':a3}

if 0:
    fig,ax=plt.subplots(1,1)
    for nrun,ht in enumerate([ht1,ht2,ht3]):
        c='rgb'[nrun]
        hull_lengths = nar(ht.hull_volumes)**(1./3)
        vals, bins = np.histogram(hull_lengths)
        bc = 0.5*(bins[1:]+bins[:-1])
        ax.plot(bc,vals/vals.sum(),color=c,label=ht.this_looper.out_prefix)
        axbonk(ax,xlabel=r'$\rm{Hull\ Length}$',ylabel=r'$\rm{N}$')

        ac = acs[ht.this_looper.out_prefix]
        ax.plot(ac.binned[1],ac.binned[2], c=c)
    ax.legend(loc=0)
    fig.savefig('plots_to_sort/hull_lengths.pdf')
    plt.close('fig')
    

if 0:
    fig,ax=plt.subplots(1,1)
    if 'ext_hull' not in dir():
        ext_hull=extents()
    for nrun,ht in enumerate([ht1,ht2,ht3]):
        odd = nar(ht.hull_volumes) < nar(ht.cell_volumes)
        not_odd = nar(ht.hull_volumes) >= nar(ht.cell_volumes)
        ax.scatter(ht.hull_volumes,ht.cell_volumes,c='k')
        ax.scatter(nar(ht.hull_volumes)[odd],nar(ht.cell_volumes)[odd],c='r')
        ax.set_aspect('equal')
        #axbonk(ax,xlabel=r'$\rm{Hull\ Volume}$',ylabel=r'$\rm{Cell\ Volume}$',xlim=[0,0.07],ylim=[0,0.07])
        ext_hull(nar(ht.hull_volumes))
        ext_hull(nar(ht.cell_volumes))
        ax.plot(ext_hull.minmax, ext_hull.minmax,c='g')
        axbonk(ax,xlabel=r'$\rm{Hull\ Volume}$',ylabel=r'$\rm{Cell\ Volume}$',
               xlim=ext_hull.minmax,ylim=ext_hull.minmax, xscale='log',yscale='log')
        fig.savefig('plots_to_sort/%s_volume_fractions.pdf'%ht.this_looper.out_prefix)

if 0:
    fig,ax=plt.subplots(1,1)
    #ax.hist(hull_lengths,histtype='step',color='k',normed=True)
    vals, bins = np.histogram(hull_lengths)
    bc = 0.5*(bins[1:]+bins[:-1])
    db = (bins[1:]-bins[:-1])
    ax.plot(bc,vals/vals.sum())
    ax.plot(AC.binned[1],AC.binned[2]/AC.binned[2][0],c='r')
    rect=patches.Rectangle((0,0),AC.L,AC.ACb[0],facecolor=[0.8]*3)
    ax.add_patch(rect)
    fig.savefig('plots_to_sort/%s_sizes.pdf'%this_simname)


