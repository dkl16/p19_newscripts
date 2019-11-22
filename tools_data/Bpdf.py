import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colorbar as cb
import matplotlib.colors as colors
import numpy as np
import yt
from importlib import reload
plt.clf()

directory = '/Users/dcollins/scratch/P49d/ca03_turb_weak'
frame=339

if 'print_reload_reminder' not in dir():
    #If you do
    #>>> import Bpdf as bp
    #it will do stuff, and the stuff it did persists, so you can do
    #>>> reload(bp)
    #and re-run, you can change how you plot without re-building the objects.
    #It's faster
    print("from importlib import reload")
    print_reload_reminder = False

#
# Create yt derived quantities.
#
def costhetaz(field,data):
    return data['magnetic_field_z']/data['magnetic_field_strength']
def thetaz(field,data):
    #arccos is behaving strangely with YTarrays
    #The zeros_like ensures we get the right type of object.
    output = np.zeros_like( data['costhetaz'])
    output[:] = np.arccos(data['costhetaz'][:].v) 
    return output

#This allows us to add these fields to any dataset
def add_cos_theta(obj):
    obj.add_field('costhetaz',costhetaz,units='dimensionless')
    obj.add_field('thetaz',thetaz,units='dimensionless')


#This conditional tests to see if region has been defined.
#if not, define it.  This process is slow so don't repeat if you don't have to.
if 'region' not in dir():
    ds = yt.load('%s/DD%04d/data%04d'%(directory,frame,frame))
    add_cos_theta(ds)
    region = ds.all_data()

#set up bins.
Nbins=64
theta_bins=np.linspace(0,np.pi,2*Nbins)
field_min = region['magnetic_field_strength'].min()
field_max = region['magnetic_field_strength'].max()
field_bins=np.logspace(np.log10(field_min),np.log10(field_max),Nbins)
bins={'thetaz':theta_bins,'magnetic_field_strength':field_bins,'magnetic_field_z':field_bins}

#
# Produce joint and marginalized distributions.
# Formally these are not PDFs since we don't explicityly normalize,
# BUT since the total volume is 1, they are in fact PDFs.
#
if 'joint' not in dir():
    bin_fields = ['magnetic_field_strength','thetaz']
    joint = yt.create_profile(region, bin_fields=bin_fields,fields=['cell_volume'], weight_field=None, override_bins=bins)
    pp = yt.PhasePlot.from_profile(joint)


if 'prof_mag' not in dir():
    prof_mag  = yt.create_profile(region, bin_fields=['magnetic_field_strength'],fields=['cell_volume'], weight_field=None, override_bins=bins)
if 'prof_theta' not in dir():
    prof_theta  = yt.create_profile(region, bin_fields=['theta_z'],fields=['cell_volume'], weight_field=None, override_bins=bins)
if 'prof_bz' not in dir():
    prof_bz  = yt.create_profile(region, bin_fields=['magnetic_field_z'],fields=['cell_volume'], weight_field=None, override_bins=bins)


#
# Now make plots with the profiles.
# 
#

if 1:
    #Marginalized distributions
    #set up plot tools
    plt.close('all')
    fig,ax=plt.subplots(1,3)
    ax0=ax[0];ax1=ax[1]; ax2=ax[2]

    #the joint probability 
    ph = joint['cell_volume']
    ph_theta = np.sum(ph,axis=0) #marginalize over magnetic field

    thetabins = prof_theta.x_bins
    thetacen = 0.5*(thetabins[1:]+thetabins[:-1])

    #
    ax0.plot(thetacen,prof_theta['cell_volume'],c='k',label=r'$P(\theta)$')
    ax0.plot(thetacen,ph_theta,c='r',label= r'$\int P(B,\theta) dB$')
    ax0.set_title('theta')
    #ax0.legend(loc=2)

    b_bins = prof_bz.x_bins
    b_cen = 0.5*(b_bins[1:]+b_bins[:-1])

    ax1.plot(b_cen,prof_mag['cell_volume'],'k',label=r'$P(B)$')
    ax1.plot(b_cen,np.sum(ph,axis=1),'r',label=r'$\int P(B,\theta) d\theta')
    ax1.set_title('B')
    ax1.set_xscale('log')




    dbin = b_cen[1:]-b_cen[:-1]
    pbz=prof_bz['cell_volume']
    dpbz=(pbz[1:]-pbz[:-1])/dbin
    ax2.plot(b_cen, b_cen*pbz)
    ax2.plot(b_cen[1:], -b_cen[1:]*dpbz,c='r')
    ax2.plot(b_cen, prof_mag['cell_volume'],c='k')
    ax2.set_xscale('log')
    plt.savefig('%s/marginalized_distributions.png'%plot_directory)
        
if 1:
    fig2, axes2 = plt.subplots(1,2)
    ax0=axes2[0]; ax1=axes2[1] #this is in case I change the number of plots in the figure.
    p_theta=prof_theta['cell_volume']
    p_mag = prof_mag['cell_volume']

    #The product joint distribution P(B)P(theta)
    p_theta.shape = (p_theta.size,1)
    ph2 = p_theta*p_mag

    ph=np.transpose(joint['cell_volume'])
    min_val = ph[ph>0].min()
    max_val = ph.max()
    norm = colors.Normalize(vmin=min_val,vmax=max_val)
    extents = [prof_mag.x_bins.min(), prof_mag.x_bins.max(),0,1]
    pl=ax0.imshow(ph, norm=norm, interpolation='nearest',origin='lower',extent=extents,aspect='auto')
    cb=fig2.colorbar(pl,ax=ax0)
    cb.cmap.set_under('w')

    norm = colors.Normalize(vmin=min_val,vmax=max_val)
    pl1=ax1.imshow(ph2, norm=norm, interpolation='nearest',origin='lower')
    cb1=fig2.colorbar(pl1,ax=ax1)
    cb1.cmap.set_under('w')

    fig2.savefig('%s/joint_test.png'%plot_directory)


