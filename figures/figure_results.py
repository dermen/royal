from pylab import *
import h5py
from matplotlib.gridspec import GridSpec

# read in the simulation
ff = np.load( '/mnt/data/cxs_simulations/corz.npz' )
c1 = ff['c1']
c2 = ff['c2']
c12= ff['c12']

# read in the data
f = h5py.File('/mnt/ssrl/silver/12dec2012/ssrl-silver-all.hdf5' ) 
cor1 = array( f['cor1'] )
cor2 = array( f['cor2'] )
cor12 = array( f['cor12'] )
cor1_mean = cor1.mean(0)
cor2_mean = cor2.mean(0)
cor12_mean = cor12.mean(0)
cor1_std = cor1.std(0)
cor2_std = cor2.std(0)
cor12_std = cor12.std(0)
cor1_err = cor1_std / sqrt( cor1.shape[0] ) 
cor2_err = cor2_std / sqrt( cor2.shape[0] ) 
cor12_err = cor12_std / sqrt( cor12.shape[0] )

# analytical peaks
peaks = [0.97, 1.248, 1.603, 1.945, 2.24]


#ticks = map( lambda x: str( int( x * 180 / pi)), peaks )
#gca().get_xaxis().set_ticks(peaks)
#gca().get_xaxis().set_ticklabels(ticks)

figure( 1, figsize=(5,7) )
gs = GridSpec( 2,3 )
ax1 = subplot( gs[0,:-1] ) # simulation
ax2 = subplot( gs[1,:-1] ) # data
ax3 = subplot( gs[:,-1]  ) # peak widths 


# plot the simulation
lw=2
color = 'b'
offset_1 = 0.005
offset_2 = 0
offset_12 = 0.018
phis = linspace( 0,2*pi,c1.shape[0])
ax1.plot(phis , c1 + offset_1 ,lw=lw,color=color)
ax1.plot(phis , c2 + offset_2 ,lw=lw,color=color)
ax1.plot(phis , c12 +offset_12,lw=lw,color=color)


lw=2   # plot analytical
alpha=0.2
ls = '-'
color='b'
for peak in peaks:
    ax1.plot(ones(2)*(peak),linspace(-1,1,2),color=color,ls=ls,lw=lw,alpha=alpha )

ax1.get_yaxis().set_ticks([])
ax1.get_xaxis().set_ticks([])
ax1.set_xlim(0.86,2.37)
ax1.set_ylim(-0.0039,0.032)
ax1.text(pi*60/180,0.032-0.005, 'A', color='r',fontsize=18)


# plot the data
offset_1 = 0.0005
offset_2 = -0.00005
offset_12 = 0.00085
color = 'b'
lw = 2
ls='d'
ms=4
phis = linspace( 0,2*pi,cor1_mean.shape[0] ) 
ax2.plot( phis, cor1_mean + offset_1,'d',color=color,ms=ms)
ax2.plot( phis, cor2_mean+offset_2,'d',color=color,ms=ms)
ax2.plot( phis, cor12_mean+offset_12,'d',color=color,ms=ms)

lw=2   # plot analytical
alpha=0.2
ls = '-'
color='b'
for peak in peaks:
    ax2.plot(ones(2)*(peak),linspace(-1,1,2),color=color,ls=ls,lw=lw,alpha=alpha )

ticks = map( lambda x: str( int( x * 180 / pi)), peaks )
ax2.get_xaxis().set_ticks(peaks)
ax2.get_xaxis().set_ticklabels(ticks)
ax2.set_xlim(0.86,2.37)
ax2.set_ylim(-0.00037,0.0013)
ax2.get_yaxis().set_ticks([])
ax2.set_xlabel(r"$\Delta$ (deg.)", fontsize=18)
ax2.text(pi*60/180,0.0013-0.0002, 'B', color='r',fontsize=18)

# plot the peak widths
ax3.get_yaxis().set_ticks([])
#gca().get_xaxis().set_ticks([])

indL_sim = int( 2.19*4096/ (2*pi) )
indR_sim = int( 2.31*4096/ (2*pi) )
indL = int( 2.19*4320/(2*pi) )
indR = int( 2.31*4320/(2*pi) )

cor = cor12[:,indL:indR]
cor_mean = cor.mean(0)
cor_std = cor.std(0)
cor_err = cor_std / sqrt( cor.shape[0] )

sim = c12[indL_sim:indR_sim] 


cor_phis = linspace( 0,2*pi, cor12_mean.shape[0]  )[indL:indR]
sim_phis = linspace( 0,2*pi, c12.shape[0]  )[indL_sim:indR_sim ]

sim = sim * cor_mean.max() / sim.max()

alpha = 0.4
print cor_phis.shape, cor_mean.shape, cor_err.shape
ax3.fill_between( cor_phis, cor_mean,  cor_mean + 1.96*cor_err, color=color,alpha=alpha )
ax3.fill_between( cor_phis, cor_mean,  cor_mean - 1.96*cor_err, color=color,alpha=alpha )
ax3.plot( cor_phis, cor_mean,'d', color=color, lw=lw)
# shift sim by one dat point to align peaks to better compare widths ( 1 dat point = 360 / 4096 degrees, so no big deal )
ax3.plot( sim_phis[1:], sim[:-1] ,'b--',color=color,lw=3 )
ticks = [2.23, 2.245, 2.258]
tick_labels = map( lambda x: '%.1f' % (x * 180 / pi), ticks )
ax3.get_xaxis().set_ticks(ticks)
ax3.get_xaxis().set_ticklabels(tick_labels)

ax3.set_ylim(-7.083e-5,0.0004385 ) 
ax3.set_xlim(2.225,2.262)
ax3.set_xlabel(r"$\Delta$ (deg.)", fontsize=18)
ax3.legend(['measured','simulated'], prop={'size':11},loc=9)

ax3.text(2.23,0.00036, 'C', color='r',fontsize=18)

show()
