from starter2 import *
import scipy.signal
import matplotlib.patches as patches

if 0:
    size=128
    dx = 1./size
    x = np.mgrid[0:1:dx]
    rho = np.zeros(size)
    rho[ x > 0.25] = 1
    rho[ x > 0.75] = 0

    AC=scipy.signal.correlate(rho,rho,mode='same',method='fft')

    fig,ax = plt.subplots(1,1)
    ax.plot(x, AC)
    fig.savefig('plots_to_sort/ac_test.png')

if 1:
    size=128
    dx = 1./size
    x,y = np.mgrid[0:1:dx, 0:1:dx]
    rmag = np.sqrt((x-0.5)**2 + (x-0.5)**2)
    rho = np.zeros_like(x)
    R_sph =  0.25
    rho[ rmag < R_sph] = 1

    AC=scipy.signal.correlate(rho,rho,mode='same',method='fft')

    import radial_binner as rb
    cen,bins = rb.rb2(rmag, AC/AC.size)


    fig,ax = plt.subplots(1,1)
    ax.plot(rmag, AC/AC.size,c='k')
    ax.plot(cen,bins,c='r')
    d = cen
    d[ d>R_sph]=R_sph
    theta = 2*np.arccos(d/R_sph)

    A = 2*R_sph**2/2*(theta-np.sin(theta))
    A = 2*R_sph**2*(np.arccos(cen) - cen*np.sqrt(1-cen**2))
    ax.plot(cen,A,c='g')
    fig.savefig('plots_to_sort/ac_test.png')
    print('wrote plots_to_sort/ac_test.png')

