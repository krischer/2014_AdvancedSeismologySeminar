__author__ = "Fabian Lindner, Florian Woelfl, 2014"


# Script for Pseudospectral method with Fourier integrals
# based on Pseudospectral lecture and Matlab code ps_fourier.m of H. Igel

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

# Parameters
nx = 1024  # space dimension
isx = 350  # source loacation
nt = 4000  # number of time steps
c = 2000.  # acoustic velocity
dx = 20.  # grid spacing in m
dt = .0015  # time increment
it0 = 30  # source time
ir = 500  # receiver location
idisp = 30  # display frequency
ag = 0. * dx  # gaussian in space for source
FD = 1  # fd solution is plotted if set to 1


# initalization of grid
x = np.arange(nx)
x = x * dx  # removed "/1000", should be in meters now
# initialization of space dependent fields
p = np.zeros(nx)  # pressure
p_fd = np.zeros(nx)
lp2 = np.zeros(nx)  # pressure 2nd derivative
lp2_fd = np.zeros(nx)
pold = np.zeros(nx)  # pressure at previous time step
pold_fd = np.zeros(nx)
pnew = np.zeros(nx)  # pressure at next time step
pnew_fd = np.zeros(nx)
# initialization of empty sesimogram
seis = np.zeros(nt) 
seis_fd = np.zeros(nt)
time = np.linspace(0, nt * dt, num=nt)


# Stability condition
print 'Stability criterion: %f (pref. < 0.2)' % (c*dt/dx)


# Source time function
src = np.zeros(nt)
# Ricker wavelet
freq = input("Please, give the source frequency in Hz: ")
t0 = 4 / (freq * np.pi)
t = np.arange(1, nt + 1) * dt
src = (np.sqrt(np.pi) / 2) * (freq * np.pi * (t - t0)) * np.exp(-(freq * np.pi * (t - t0))**2)
src = np.diff(src)
src = src / np.max(src)
# append one entry as one is lost in differentiation
src = np.append(src, 0.0)
# spatial distribution of source
gauss = np.zeros(nx)
if ag == 0.:
    gauss[isx] = 1
else:
    gauss = np.exp(-1/ag**2 * (x - x[isx]) ** 2)


# function for calculating the 1st spatial derivative
def sder1d(f,dx):
    # initialize k
    kmax = np.pi / dx
    dk = kmax / (nx / 2)
    k = np.arange(float(nx))
    k[:nx/2] = k[:nx/2] * dk 
    k[nx/2:] = k [:nx/2] - kmax
    k = 1j * k

    # FFT and IFFT
    ff = np.fft.fft(f)
    ff = k * ff
    df = np.real(np.fft.ifft(ff))
    return df
# function for calculating the 2nd spatial derivative

def sder2d(f,dx):
    # initialize k
    kmax = np.pi / dx
    dk = kmax / (nx / 2)
    k = np.arange(float(nx))
    k[:nx/2] = k[:nx/2] * dk 
    k[nx/2:] = k [:nx/2] - kmax
    k = (1j * k)**2

    # FFT and IFFT
    ff = np.fft.fft(f)
    ff = k * ff
    df = np.real(np.fft.ifft(ff))
    return df


# initialize plots
plt.ion()
fig = plt.figure(figsize=(15,10))
gs = gridspec.GridSpec(2,2,width_ratios=[1,1],height_ratios=[1,2],hspace=0.4)
ax1 = plt.subplot(gs[0])
length = 2.5 * 1./float(freq)
frac = int((length/time[-1]) * nt)
ax1.plot(time[:frac], src[:frac])
ax1.set_title('Source time function')
ax1.set_xlabel('[s]')

ax2 = plt.subplot(gs[1])
ax2.plot(x[isx-10:isx+10],gauss[isx-10:isx+10])
ax2.set_title('Spation distribution of source')
ax2.set_xlabel('[m]')

ax3 = plt.subplot(gs[2])
upd31, = ax3.plot(x,p)
leg1, = ax3.plot(x[isx],0,'r*',markersize=11)
leg2, = ax3.plot(x[ir],0,'k^',markersize=8)
leg3, = ax3.plot([x[0],x[-1]],[0,0], 'k')
if FD == 1:
   upd32, = ax3.plot(x,p_fd,'g')
   ax3.legend((leg1,leg2,upd31,upd32), ('Source','Receiver','PS','FD'), loc='upper right',numpoints=1)
else:
   ax3.legend((leg1,leg2), ('Source','Receiver'), loc='upper right',numpoints=1)
lim = 1./float(freq) * 10000
ax3.set_ylim(-lim,lim)
ax3.set_xlim(x[0],x[-1])
ax3.set_title('Propagating wave')
ax3.set_xlabel('Distance [m]')

ax4 = plt.subplot(gs[3])
upd41, = ax4.plot(time,seis)
upd42, = ax4.plot([0],[0],'b|',markersize=15)
if FD == 1:
    upd43, = ax4.plot(time,seis_fd,'g')
ax4.set_ylim(-lim,lim)
ax4.set_xlim(time[0],time[-1])
ax4.set_title('Seismogram at receiver')
ax4.set_xlabel('Time [s]')
plt.show()


# Begin time extrapolation
print 'Time exptrapolation ...'
for it in xrange(nt):
    # pseudospectral - fourier
    # laplacian with spectral derivative
    lp2 = sder2d(p,dx)
    # calculate two times first derivative
    #lp1 = sder1d(p,dx)
    #lp2 = sder1d(lp1,dx)
    # add source and compute new pressure field
    pnew = 2 * p - pold + c**2 * dt**2 * (lp2 + gauss * src[it])

    pold = p
    p = pnew

    # output
    seis[it] = p[ir]

    # Calculate also FD solution
    # second spatial derivative
    for i in xrange(2,nx-2):
        lp2_fd[i] = -1./12.*p_fd[i+2] + 4./3.*p_fd[i+1] - 5./2.*p_fd[i] + 4./3.*p_fd[i-1] - 1./12.*p_fd[i-2]
    lp2_fd /= dx**2
    # add source and compute new pressure field
    pnew_fd = 2 * p_fd - pold_fd + c**2 * dt**2 * (lp2_fd + gauss * src[it])

    pold_fd = p_fd
    p_fd = pnew_fd

    # output
    seis_fd[it] = p_fd[ir]
    

    # update data for plot
    if (it % idisp) is 0:
        upd31.set_ydata(p)
        upd41.set_ydata(seis)
        upd42.set_data(time[it],seis[it])
        # plot fd solution only if FD is set to 1
        if FD == 1:
            upd32.set_ydata(p_fd)
            upd43.set_ydata(seis_fd)
        fig.canvas.draw()
