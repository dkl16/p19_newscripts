

from starter2 import *
import product
reload(product)
import make_page
reload(make_page)

basedir = "/Users/dcollins/RESEARCH3/Paper19c_u10u11/0000_more_good_plots/"
for this_simname in ['u05','u10','u11']:
    glob1 = "%s/density_time/%s/%s_density_6_c????.png"%(basedir,this_simname,this_simname)
    reg1 = re.compile(r'%s/density_time/%s/%s_density_6_c(\d\d\d\d).png'%(basedir,this_simname,this_simname))
    p1 = product.product('rho_t', regexp=reg1, myglob=glob1, parameters=['core_id'],style='single',width=400)
    p1.get_frames()

    p13 = product.product('Tc_raw',fname='browser/%s_ct.h5'%this_simname,field='collapse_times',style='value',width=400)


    product_list=[p13,p1]
    cl=make_page.make_page(product_list, core_list=None,htmlname='browser/output_%s.html'%(this_simname))


