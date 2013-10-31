from pylab import *
import h5py

f = h5py.File('/mnt/ssrl/silver/12dec2012/ssrl-silver-all.hdf5' ) 

cor1 = array( f['cor1'] )
cor2 = array( f['cor2'] )
cor12 = array( f['cor12'] )

cor1_mean = cor1.mean(0)
cor2_mean = cor2.mean(0)
cor12_mean = cor12.mean(0)

step = 6

cor1_mean = convolve(cor1_mean, ones( (step,))/step )[step-1:][::step]
cor2_mean = convolve(cor2_mean, ones( (step,))/step )[step-1:][::step]
cor12_mean = convolve(cor12_mean, ones( (step,))/step )[step-1:][::step]

cor1_std = cor1.std(0)
cor2_std = cor2.std(0)
cor12_std = cor12.std(0)

# im pretty sure I need a factor of root something here... probs neglig

cor1_std = convolve(cor1_std, ones( (step,))/step )[step-1:][::step]
cor2_std = convolve(cor2_std, ones( (step,))/step )[step-1:][::step]
cor12_std = convolve(cor12_std, ones( (step,))/step )[step-1:][::step]

cor1_err = cor1_std / sqrt( cor1.shape[0] ) 
cor2_err = cor2_std / sqrt( cor2.shape[0] ) 
cor12_err = cor12_std / sqrt( cor12.shape[0] )


figure( 1, figsize=(4.5,7) )

subplot(221)

phis = linspace( 0,2*pi,cor1_mean.shape[0] ) 
offset_1 = 0.0004
offset_2 = -0.00005
offset_12 = 0.00068

#data means
color = 'b'
lw = 2
plot( phis, cor1_mean + offset_1,lw=lw,color=color)
plot( phis, cor2_mean+offset_2,lw=lw,color=color)
plot( phis, cor12_mean+offset_12,lw=lw,color=color)

# data errors
alpha = 0.4
color = 'b'
fill_between( phis, cor1_mean+offset_1,  cor1_mean + 1.96*cor1_err+offset_1, color=color,alpha=alpha )
fill_between( phis, cor1_mean+offset_1,  cor1_mean - 1.96*cor1_err+offset_1, color=color,alpha=alpha )
fill_between( phis, cor2_mean+offset_2,  cor2_mean + 1.96*cor2_err + offset_2,color=color,alpha=alpha )
fill_between( phis, cor2_mean+offset_2,  cor2_mean - 1.96*cor2_err + offset_2,color=color,alpha=alpha )
fill_between( phis, cor12_mean + offset_12, cor12_mean + 1.96*cor12_err + offset_12,color=color,alpha=alpha )
fill_between( phis, cor12_mean+offset_12, cor12_mean - 1.96*cor12_err + offset_12,color=color,alpha=alpha )

# analytical
lw=2
alpha=0.2
ls = '-'
color='b'
plot(ones(2)*(1.248),linspace(-1,1,2),color=color,ls=ls,lw=lw,alpha=alpha )
plot(ones(2)*(1.945),linspace(-1,1,2),color=color,ls=ls,lw=lw,alpha=alpha )
plot(ones(2)*(1.603),linspace(-1,1,2),color=color,ls=ls,lw=lw,alpha=alpha )
plot(ones(2)*(0.97),linspace(-1,1,2),color=color,ls=ls,lw=lw,alpha=alpha )
plot(ones(2)*(2.24),linspace(-1,1,2),color=color,ls=ls,lw=lw,alpha=alpha )

xlim(0.86,2.37)
ylim(-0.00037,0.00105)
gca().get_yaxis().set_ticks([])
gca().get_xaxis().set_ticks([])

subplot(223)

ff = np.load( '/mnt/data/cxs_simulations/corz.npz' )
c1 = ff['c1']
c2 = ff['c2']
c12= ff['c12']

lw=2
color = 'b'
offset_1 = 0.005
offset_2 = 0
offset_12 = 0.018
phis = linspace( 0,2*pi,c1.shape[0])
plot(phis , c1 + offset_1 ,lw=lw,c=color)
plot(phis , c2 + offset_2 ,lw=lw,c=color)
plot(phis , c12 +offset_12,lw=lw,c=color)

#analytical
lw=2
alpha=0.2
ls = '-'
color='b'
plot(ones(2)*(1.248),linspace(-1,1,2),color=color,ls=ls,lw=lw,alpha=alpha )
plot(ones(2)*(1.945),linspace(-1,1,2),color=color,ls=ls,lw=lw,alpha=alpha )
plot(ones(2)*(1.603),linspace(-1,1,2),color=color,ls=ls,lw=lw,alpha=alpha )
plot(ones(2)*(0.97),linspace(-1,1,2),color=color,ls=ls,lw=lw,alpha=alpha )
plot(ones(2)*(2.24),linspace(-1,1,2),color=color,ls=ls,lw=lw,alpha=alpha )

gca().get_yaxis().set_ticks([])
xlim(0.86,2.37)
ylim(-0.0039,0.032)

peaks = [0.97, 1.248, 1.603, 1.945, 2.224]
ticks = map( lambda x: str( int( x * 180 / pi)), peaks )
gca().get_xaxis().set_ticks(peaks)
gca().get_xaxis().set_ticklabels(ticks)

xlim(0.86,2.37)
xlabel(r"$\Delta$ (deg.)", fontsize=18)

subplot( 222 ) 
gca().get_yaxis().set_ticks([])
gca().get_xaxis().set_ticks([])

indL_sim = int( 2.19*4096/ (2*pi) )
indR_sim = int( 2.31*4096/ (2*pi) )
indL = int( 2.19*4320/(2*pi) )
indR = int( 2.31*4320/(2*pi) )

zoom12 = cor12.mean(0)[indL:indR]
zoom12_sim = c12[indL_sim:indR_sim] 


# take off 3 pixels for alignment ( we will mention discrepancy in paper, but this fig is to illustrate difference in width and works best if both peaks are aligned
lw = 2
color = 'b'
plot( (zoom12 / zoom12.max())[3:],color=color, lw=lw,ls='-'  )
plot( zoom12_sim / zoom12_sim.max(),color=color,lw=lw,ls='--' )
ylim(-0.1,1.025)
xlim(20,50)
legend(['measured','simulated'], prop={'size':11},loc=4)


subplot( 224)

phis = linspace( 0,2*pi,cor1_mean.shape[0] ) 
offset_1 = 0.0004
offset_2 = -0.00005
offset_12 = 0.00068

color = 'b'
lw = 2
plot( phis, cor1_mean + offset_1,lw=lw,color=color)
plot( phis, cor2_mean+offset_2,lw=lw,color=color)
plot( phis, cor12_mean+offset_12,lw=lw,color=color)

# data errors
alpha = 0.4
color = 'b'
fill_between( phis, cor1_mean+offset_1,  cor1_mean + 1.96*cor1_err+offset_1, color=color,alpha=alpha )
fill_between( phis, cor1_mean+offset_1,  cor1_mean - 1.96*cor1_err+offset_1, color=color,alpha=alpha )
fill_between( phis, cor2_mean+offset_2,  cor2_mean + 1.96*cor2_err + offset_2,color=color,alpha=alpha )
fill_between( phis, cor2_mean+offset_2,  cor2_mean - 1.96*cor2_err + offset_2,color=color,alpha=alpha )
fill_between( phis, cor12_mean + offset_12, cor12_mean + 1.96*cor12_err + offset_12,color=color,alpha=alpha )
fill_between( phis, cor12_mean+offset_12, cor12_mean - 1.96*cor12_err + offset_12,color=color,alpha=alpha )


#ticks = [str(int( pi/2 * 180 / pi ) ), str( int( ( 0.1 + 6.18 )* 90 / pi ) ) ,  str(int( 6.18 * 180 / pi ))]
ticks = ["90", "180", "270"]
gca().get_xaxis().set_ticks(ticks)
ylim(-0.00037,0.00105)
gca().get_yaxis().set_ticks([])
gca().get_xaxis().set_ticks([ pi/2, pi, 3*pi/2])
gca().get_xaxis().set_ticklabels(ticks)

xlim(0.1,2*pi-0.1)
xlabel(r"$\Delta$ (deg.)", fontsize=18)


show()
