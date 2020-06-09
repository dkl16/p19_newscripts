from starter2 import *
import scipy.signal
import matplotlib.patches as patches
from PIL import Image
plt.close('all')
if 0:
    size=128
    dx = 1./size
    x = np.mgrid[0:1:dx]
    rho = np.zeros(size)
    rho[ x > 0.25] = 1
    rho[ x > 0.75] = 0

if 1:
    rho = np.array(Image.open("circ.tif"))


    #AC=scipy.signal.correlate(rho,rho,mode='same',method='fft')
    AC = np.fft.fftn(rho)
    AC = AC*np.conj(AC)
    AC = np.fft.ifftn(AC)
    AC = np.fft.fftshift(AC)
    AC = np.real(AC)
    AC = AC/AC.max()

    fig,ax = plt.subplots(1,1)
    size=64
    dx = 1./size
    x,y = np.mgrid[0:1:dx, 0:1:dx]
    rmag = np.sqrt((x-0.5)**2 + (y-0.5)**2)

    #ax.plot(rmag[33,33:],  AC[33,33:],c='k')
    ax.plot(  AC[33,33:],c='k')
    r1 = np.arange(31)/32 #rmag[33,33:]
    y = 2.0*(np.arccos(r1)- r1*np.sqrt(1.0-r1**2))/np.pi;
    ax.plot( y)

    fig.savefig('plots_to_sort/AC1.png')

if 0:
    size=64
    dx = 1./size
    x,y = np.mgrid[0:1:dx, 0:1:dx]
    rmag = np.sqrt((x-0.5)**2 + (x-0.5)**2)
    rho = np.zeros_like(x)
    R_sph =  0.25
    rho[ rmag < R_sph] = 1


    import radial_binner as rb
    cen,bins = rb.rb2(rmag, AC/AC.size)


    fig,ax = plt.subplots(1,1)
    #ax.plot(rmag, AC/AC.size,c='k')
    ax.plot(cen,bins,c='r')
    ax.plot(np.linspace(0,1,31),AC[33,33:]/AC.max(),c='g')
    d = cen
    d[ d>R_sph]=R_sph
    theta = 2*np.arccos(d/R_sph)

    A = 2*R_sph**2/2*(theta-np.sin(theta))
    A = 2*R_sph**2*(np.arccos(cen) - cen*np.sqrt(1-cen**2))
    ax.plot(cen,A,c='g')
    fig.savefig('plots_to_sort/ac_test.png')
    print('wrote plots_to_sort/ac_test.png')

