"""
Applications of looper tools.

Decorators:  The @looper.frame_loop and @looper.core_loop decorate functions.
@looper.frame_loop
def this_code(this_looper,other_args):
    ...
will cause this_code to be executed on each frame in this_looper.frame_list.

@looper.core_loop   
def this_code(this_looper):
    stuff

will cause the this_code to be executed for each snapshot on each frame.
"""
from starter2 import *


@looper.core_loop   #the decorator takes care of the loop
def print_centroid(looper,snapshot): #here's a function, needs to take at least these arguments)
    print("Core %d frame %d centroid (%s)"%(
          snapshot.core_id, snapshot.frame, str( snapshot.R_centroid)))

@looper.frame_loop
def proj_cores(self, axis_list=[0,1,2],core_list=[], field='density'):
    for axis in axis_list:
        center = self.ds.arr(nar([0.5]*3),'code_length')
        self.proj = self.ds.proj(field,axis,center=center)
        self.proj = yt.SlicePlot(self.ds, axis=axis, fields=[field], center=center)
        radius_from_core = []
        core_label = ""
        for nc,core_number in enumerate(core_list):
            core_label += "c%04d_"%core_number
            self.proj.annotate_select_particles(1.0, col='r', indices=self.target_indices[core_number])
            outname = '%s_full_particles_%sn%04d'%(self.out_prefix,core_label,self.current_frame)
        print( self.proj.save(outname))

@looper.frame_loop
def select_particles(looper,these_particles=None,axis_list=[0,1,2]):
    for axis in axis_list:
        proj = looper.ds.proj('density',axis,center='c')
        pw = proj.to_pw(center = 'c',width=(1.0,'code_length'), origin='domain')
        pw.set_cmap('density','gray')
        radius_from_core = []
        core_label = ""
        pw.annotate_select_particles(1.0, col='r', indices=these_particles, p_size=5)
        outname = '%s_select_n%04d'%(looper.out_prefix,looper.current_frame)
        print(outname)
        print( pw.save(outname))

@looper.core_loop
def core_proj_follow(looper,snapshot, field='density', axis_list=[0,1,2], color='r'):
    for ax in axis_list:
        scale_min = snapshot.ds.arr(0.05,'code_length')
        scale = max([2*snapshot.R_mag.max(), scale_min])
        sph = snapshot.ds.sphere(snapshot.R_centroid,scale)
        proj = snapshot.ds.proj(field,ax,center=snapshot.R_centroid, data_source = sph)
        pw = proj.to_pw(center = snapshot.R_centroid,width=(1.0,'code_length'), origin='domain')
        pw.zoom(1./scale.v)
        pw.set_cmap(field,'gray')
        pw.annotate_sphere(snapshot.R_centroid,snapshot.R_mag.max(), circle_args={'color':color} ) #R_mag.max())
        pw.annotate_text(snapshot.R_centroid,
                         "%d"%snapshot.core_id,text_args={'color':color}, 
                         inset_box_args={'visible':False},
                         coord_system='data')
        pw.annotate_select_particles(1.0, col=color, indices=snapshot.target_indices)
        outname = "%s_c%04d_n%04d_centered"%(looper.out_prefix,snapshot.core_id,snapshot.frame)
        print(pw.save(outname))
