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


@looper.frame_loop
def proj_onecore(self, axis_list=[0,1,2],core_list=[], field='density'):
    for axis in axis_list:
        #center=self.ds.arr(nar([0.4,0.7,0.8]),'code_length')
        center = self.ds.arr(nar([0.5]*3),'code_length')
        #self.ds.periodicity = (True,True,True)
        self.proj = self.ds.proj(field,axis,center=center)
        self.proj = yt.SlicePlot(self.ds, axis=axis, fields=[field], center=center)


        #pw = self.proj.to_pw(center = center,width=(1.0,'code_length'), origin='domain')
        #pw.set_cmap(field,'gray')
        radius_from_core = []
        core_label = ""
        for nc,core_number in enumerate(core_list):
            core_label += "c%04d_"%core_number
            self.proj.annotate_select_particles(1.0, col='r', indices=self.target_indices[core_number])
            outname = '%s_full_particles_%sn%04d'%(self.out_prefix,core_label,self.current_frame)
        print( self.proj.save(outname))
