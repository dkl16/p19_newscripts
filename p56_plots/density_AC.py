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



if 0:
    if 0:   
        """this is what you want to do."""
        directory = '/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128/GravPotential'
        ds = yt.load("%s/DD%04d/data%04d"%(directory,0,0))
        prefix='u05'
    if 1:
        """a test"""
        directory = "/scratch1/dcollins/Paper19/u21_sphere_large"
        ds = yt.load("%s/DD%04d/data%04d"%(directory,0,0))
        prefix='u21'
    if 1:
        """a smaller test test"""
        directory = "/scratch1/dcollins/Paper19/u22_sphere_128"
        ds = yt.load("%s/DD%04d/data%04d"%(directory,0,0))
        prefix='u22'

    left=[0.0]*3
    resolution = ds['TopGridDimensions'] 
    cg=ds.covering_grid(0,left,resolution)#,num_ghost_zones=num_ghost_zones)

    rho1 = cg['density'].v -1 #[:40,:40,:40]
    rho2=rho1
else:
    dx = 1./128
    R_sphere = 0.25
    x,y,z = np.mgrid[0:1:dx,0:1:dx, 0:1:dx]
    r = np.sqrt( (x-0.5)**2+(y-0.5)**2+(z-0.5)**2)
    rho1 = np.zeros_like(x)
    rho1[ r<R_sphere] = 0.5
    rho2=rho1
    prefix="test"

if 1:
    fig2,axes2=plt.subplots(2,2,figsize=figsize)
    a20,a21=axes2[0]
    a23,a26=axes2[1]

    fig3,a24=plt.subplots(1,1,figsize=figsize)

    a20.imshow( rho2.sum(axis=axis), cmap='Greys' )

if 1:
    if 'AC3d' not in dir() or True:
        print('Correlate')
        AC3d=scipy.signal.correlate(rho1,rho2,mode='same',method='fft')
        AC3d=np.roll(AC3d, AC3d.shape[0]//2,axis=0)
        AC3d=np.roll(AC3d, AC3d.shape[1]//2,axis=1)
        AC3d=np.roll(AC3d, AC3d.shape[2]//2,axis=2)
        print('rolled')
              

    print( 'fft')
    rhohat = np.fft.ifftn(rho1)
    rhohatdag = rhohat.conj()
    rho_prod = rhohat*rhohatdag
    print( 'fft2')
    AC3dft = np.fft.fftn(rho_prod)

    k_array = make_k_freqs( AC3dft.shape[0],real=False)
    kmag = np.sqrt(k_array[0,...]**2 + k_array[1,...]**2 + k_array[2,...]**2)
    ktt=kmag.flatten()


    a21.imshow( (AC3d.real).sum(axis=axis) , cmap='Greys')
    proj = AC3dft.real.sum(axis=axis)
    a23.imshow(proj,cmap='Greys')

if 1:    
    bins = np.fft.fftfreq(AC3dft.shape[0])
    ok = bins > 0
    bins = np.r_[0,bins[ok]]
    db = bins[1:]-bins[:-1]
    print('binning')
    binned=rb.rb( k_array, AC3dft.real,bins=bins)
    binned2=rb.rb( k_array, AC3d.real/AC3d.size,bins=bins)

if 1:
    a24.plot( binned[0],binned[1],c='r',label='binned ACfft')
    a24.plot( binned2[0],binned2[1],c='g',label=r'$AC_\rho(r)$')

    r_coord = binned[0]
    h = R_sphere-r_coord
    h[ h<0] = 0
    V = np.pi*h**2/3*(3*R_sphere-h)
    AC_m = 2*V
    a24.plot( r_coord, AC_m)

    ok = binned2[1]>0

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

if 1:
    """actual correlation length"""

    AC = binned2[1]
    L = np.sum(AC*db)/AC[0]

if 1:
    """do we want the rectangle?"""
    rect=patches.Rectangle((0,0),L,AC[0],facecolor=[0.2]*3)
    a24.add_patch(rect)

if 1:
    axbonk(a24,xlabel='$r$', ylabel=r'$\rm{AC}_\rho(r)$',
           #yscale='log',xscale='log')
           yscale='linear',xscale='linear')
    #a24.set_yscale('log')
    #x0a24.set_xscale('log')
    

if 0:
    a24.legend(loc=0)

fig2.savefig('plots_to_sort/p56_convolve_4_%s.png'%prefix)
fig3.savefig('plots_to_sort/p56_convolve_5_%s.png'%prefix)
print('saved')

