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
reload(looper)


def cbump(value):
    """for bumping the centroid value to between 0,1.
    """
    output = (value+0).v
    output[output>1] = output[output>1]-1
    output[output<0] = output[output<0]+1
    return output
@looper.core_loop   #the decorator takes care of the loop
def print_centroid(looper,snapshot): #here's a function, needs to take at least these arguments)
    print("Core %d frame %d centroid (%s)"%(
          snapshot.core_id, snapshot.frame, str( snapshot.R_centroid)))

@looper.frame_loop
def proj_cores(self, axis_list=[0,1,2],core_list=[], field='density'):
    """Full projections of the data, with core particles marked."""
    for axis in axis_list:
        center = self.ds.arr(nar([0.5]*3),'code_length')
        self.proj = self.ds.proj(field,axis,center=center)
        self.proj = yt.ProjectionPlot(self.ds, axis=axis, fields=[field], center=center)
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
    if snapshot.R_centroid is None:
        snapshot.get_all_properties()
    ds = snapshot.get_ds()
    for ax in axis_list:
        center = ds.arr(snapshot.R_centroid,'code_length')
        Rmax = snapshot.R_mag.max()
        scale_min = ds.arr(0.05,'code_length')
        scale = max([2*Rmax, scale_min])
        sph = ds.sphere(center,scale)
        proj = ds.proj(field,ax,center=center, data_source = sph) 
        pw = proj.to_pw(center = center,width=(1.0,'code_length'), origin='domain')
        pw.zoom(1./scale.v)
        pw.set_cmap(field,'gray')
        pw.annotate_sphere(center,Rmax, circle_args={'color':color} ) #R_mag.max())
        pw.annotate_text(center,
                         "%d"%snapshot.core_id,text_args={'color':color}, 
                         inset_box_args={'visible':False},
                         coord_system='data')
        pw.annotate_select_particles(1.0, col=color, indices=snapshot.target_indices)
        outname = "%s_c%04d_n%04d_centered"%(looper.out_prefix,snapshot.core_id,snapshot.frame)
        print(pw.save(outname))

@looper.core_loop
def core_circle(looper,snapshot, axis_list=[0,1,2]):
    for axis in axis_list:
        proj = looper.ds.proj('density',axis,center='c')
        pw = proj.to_pw(center = 'c',width=(1.0,'code_length'), origin='domain')
        pw.set_cmap('density','gray')
        radius_from_core = []
        core_label = ""
        pw.annotate_select_particles(1.0, col='r', indices=looper.target_indices[looper.core_id])
        pw.annotate_sphere(snapshot.R_centroid,snapshot.R_mag.max() ) #R_mag.max())
        outname = '%s_fullcircle_cores_%sn%04d'%(looper.out_prefix,looper.core_id,looper.current_frame)
        print( pw.save(outname))

@looper.frame_loop
def proj_cores_with_annotations(self, axis_list=[0,1,2], color_dict={}):
    """Full projections, with particles plotted.  Also plots the radius and number for each core.
    Colors can be passed in through a dictionary, indexed by core_id"""
    for axis in axis_list:
        proj = self.ds.proj('density',axis,center='c')
        pw = proj.to_pw(center = 'c',width=(1.0,'code_length'), origin='domain')
        pw.set_cmap('density','gray')
        radius_from_core = []
        core_label = ""
        for nc,core_number in enumerate(self.target_indices):
            this_snapshot = self.make_snapshot(self.current_frame,core_number)
            if this_snapshot.R_centroid is None:
                this_snapshot.get_all_properties()
            center = this_snapshot.R_centroid
            color = color_dict.get(core_number,'r')
            core_label += "c%04d_"%core_number
            pw.annotate_text(cbump(center),
                             "%d"%core_number,text_args={'color':color}, 
                             inset_box_args={'visible':False},
                             coord_system='data')
            pw.annotate_select_particles(1.0, col=color, indices=self.target_indices[core_number])
        outname = '%s_full_particles_%sn%04d'%(self.out_prefix,core_label,self.current_frame)
        print( pw.save(outname))
