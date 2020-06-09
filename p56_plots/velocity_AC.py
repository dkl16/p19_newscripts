from starter2 import *
import scipy.signal
import matplotlib.patches as patches
plt.close('all')
figsize = None #(12,12)
import tools.radial_binner as rb
reload(rb)
from scipy.optimize import curve_fit
def powerlaw(r,rho0, r0, alpha):
    return alpha*np.log10(r/r0) + np.log10(rho0)
def powerlaw2(r, r0,alpha):
    rhosquared=1
    return alpha*np.log10(r/r0) + np.log10(rhosquared)
axis=0
def make_k_freqs(nk,real=False, d=1):
    ny = nk
    nz = nk
    if real:
        nz = nk//2+1

    k_freq = np.zeros([3,nk,ny,nz])
    k1=np.fft.fftfreq(nk,d=d)
    #kx, ky, kz = np.meshgrid(k1,k1,k1)
    #k_freq[0,...]=kx
    #k_freq[1,...]=ky
    #k_freq[2,...]=kz
    x = np.repeat(k1,nk*nk)
    x.shape = (nk,nk,nk)
    y = np.repeat(k1,nk*nk)
    y.shape = (nk,nk,nk)
    y=y.swapaxes(0,1)
    z = np.repeat(k1,nk*nk)
    z.shape = (nk,nk,nk)
    z=z.swapaxes(0,2)
    if real:
        x = x[:,:,:nz]
        y = y[:,:,:nz]
        z = z[:,:,:nz]
    k_freq[0,...]=x
    k_freq[1,...]=y
    k_freq[2,...]=z
    return k_freq



if 'vx' not in dir():
    directory = '/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128/GravPotential'
    ds = yt.load("%s/DD%04d/data%04d"%(directory,0,0))
    left=[0.0]*3
    resolution = ds['TopGridDimensions'] 
    cg=ds.covering_grid(0,left,resolution)#,num_ghost_zones=num_ghost_zones)

    #rho1 = cg['density'].v #-1 #[:40,:40,:40]
    #rho2=rho1
    vx = cg['velocity_x']
    vy = cg['velocity_y']
    vz = cg['velocity_z']

if 1:
    fig2,axes2=plt.subplots(2,2,figsize=figsize)
    a20,a21=axes2[0]
    a22,a23=axes2[1]

    fig3,ax2=plt.subplots(1,1,figsize=figsize)

if 1:
    if 'ACx' not in dir():
        print('Correlate')
        ACx=scipy.signal.correlate(vx,vx,mode='same',method='fft')
        ACy=scipy.signal.correlate(vy,vy,mode='same',method='fft')
        ACz=scipy.signal.correlate(vz,vz,mode='same',method='fft')


    a20.imshow( (ACx.real).sum(axis=axis) , cmap='Greys')
    a21.imshow( (ACy.real).sum(axis=axis) , cmap='Greys')
    a22.imshow( (ACz.real).sum(axis=axis) , cmap='Greys')

if 1:    
    
    bins = np.arange(0.5*dx,1,dx)
    db = bins[1:]-bins[:-1]
    dx = 1./128
    x0 = -0.5
    rall=np.mgrid[0.5*dx+x0:1+x0:dx,0.5*dx+x0:1+x0:dx,0.5*dx+x0:1+x0:dx]

if 1:
    print('binning')
    v2a = ACx+ACy+ACz
    v2 = 2*9-2*v2a
    binned2=rb.rb( rall, v2/v2.size,bins=bins)
    a23.plot( binned2[0], binned2[1])
    fig2.savefig('plots_to_sort/vel.png')

if 0:
    #a24.plot( binned[0],binned[1],c='r',label='binned ACfft')
    a24.plot( binned2[0],binned2[1],c='g',label=r'$AC_\rho(r)$')

    ok = binned2[1]>0
    #fits = np.polyfit(np.log10(binned2[0][ok]), np.log10(binned2[1][ok]),2)
    #a24.plot( binned2[0][ok], binned2[0][ok]**fits[-1])

if 0:
    """fit to power laws.  Clearly this is not a great thing to do."""
    popt, pcov = curve_fit(powerlaw, binned2[0][ok][:20],np.log10(binned2[1][ok][:20]), p0=[1,1,-2])
    fit_rho0, fit_r0, fit_alpha = popt

    popt2, pcov2 = curve_fit(powerlaw2, binned2[0][ok][:20],np.log10(binned2[1][ok][:20]), p0=[1,-2])
    fit_r02, fit_alpha2 = popt2

    rrr = binned2[0]
    acac = binned2[1]
    a24.plot( rrr, 10**powerlaw(rrr, fit_rho0, fit_r0, fit_alpha),label='powerlaw')
    a24.plot( rrr, 10**powerlaw2(rrr, fit_r02, fit_alpha2),label='tweak2')
    a24.plot( rrr, 10**powerlaw2(rrr, L, fit_alpha2),label='tweak2')

if 0:
    """actual correlation length"""

    AC = binned2[1]
    L = np.sum(AC*db)/AC[0]

if 0:
    """do we want the rectangle?"""
    rect=patches.Rectangle((0,0),L,AC[0],facecolor=[0.2]*3)
    a24.add_patch(rect)

if 0:
    axbonk(a24,xlabel='$r$', ylabel=r'$\rm{AC}_\rho(r)$',
           #yscale='log',xscale='log')
           yscale='linear',xscale='linear')
    #a24.set_yscale('log')
    #x0a24.set_xscale('log')
    

if 0:
    a24.legend(loc=0)

fig2.savefig('plots_to_sort/p56_convolve_4.png')
fig3.savefig('plots_to_sort/p56_convolve_5.png')
print('saved')

