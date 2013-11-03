from odin.xray import parse
from sys import argv,exit
from scipy.interpolate import RectBivariateSpline as RBS
from scipy.interpolate import UnivariateSpline as US
import numpy as np
import pylab as plt
import h5py
import os
import collections as coll
import matplotlib.cm as cm
from matplotlib.gridspec import GridSpec
from pylab import *
def make_ticklabels_invisible(fig):
    for i, ax in enumerate(fig.axes):
        #ax.text(0.5, 0.5, "ax%d" % (i+1), va="center", ha="center")
        for tl in ax.get_xticklabels() + ax.get_yticklabels():
            tl.set_visible(False)

def spline_eval(spline, interp_points ):
    interped_vals = []
    for i in xrange( interp_points.shape[0] ):
        interped_vals.append(  spline( *interp_points[i] )[0][0] )
    interped_vals=  np.array ( interped_vals )
    return interped_vals

def common_mode ( panel  ):
    intens_hist = np.histogram( panel.flatten(), bins = 10000 )
    intens_offset = (intens_hist[1][ np.argmax(intens_hist[0])] +\
                    intens_hist[1][ np.argmax(intens_hist[0]) +1  ] ) / 2.
    panel -= intens_offset
    return panel

def depolarize(panel,r_vals,phi_vals, wavelen,pixsize,detdist, pol_factor):
    theta_vals = np.arctan2( r_vals*pixsize ,detdist  )  #/ 2.
    norm  = pol_factor*( 1 - (np.sin(theta_vals)**2) * (np.cos( phi_vals )**2) )
    norm  += (1-pol_factor) * ( 1 - (np.sin(theta_vals)**2) * (np.sin( phi_vals )**2) )
    panel /= norm
    return panel

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

def window_rms(a1, window_size):
    a2 = np.power(a1,2)
    window = np.ones(window_size)/float(window_size)
    return np.sqrt(np.convolve(a2, window, 'valid'))

cbf = parse.CBF( '/mnt/ssrl/silver/12dec2012/cbf/A5/p3d-p7s-188-10/p3deg_p7s_26_1_00360.cbf' ) 
img = cbf.intensities
img2 = copy( img ) 
wavelen = cbf.wavelength
pixsize = cbf.pixel_size[0]
detdist = cbf.path_length
    
# scan parameters

a,b   = 1231.25, 1264.55 # A5 / A6
q     = 351.95 # q1 = 411.6
panels=[22,23,24,  ## DETECTOR DIST = 188
        27,   29,
        32,   34,
        37,38,39]
masked_1, masked_2 = 4,6 # these are the asic indices that are partially shadowed by the beamstop

# DETECTOR PARAMETERS
X,Y        = np.meshgrid( np.arange(2463),np.arange( 2527) ) 
edge       = 3
R          = np.sqrt( (X-a)**2 + ( Y-b )**2 )
PHI        = np.arctan2( ( Y-b) , ( X-a ) )
num_phi    = 4320
phi_values = np.linspace( -np.pi, np.pi, num_phi  )
ring       = np.array( [ [ q*np.sin(phi)+b, q*np.cos( phi) + a] for phi in phi_values ] )


## PANEL PARAMETERS
pan_ind_y    = np.array([0, 212, 424, 636, 848,
                    1060,1272,1484,1696,
                    1908,2120,2332])
pan_ind_x    = np.array([0, 494, 988, 1482, 1976])
pan_ind_y   += edge
pan_ind_x   += edge                                              
ydim, xdim   = 195-edge,487-edge
###############

num_pan = len( panels ) 
pan_array = np.zeros( ( num_pan, ydim, xdim)  )
x_array   = np.zeros_like( pan_array)
y_array   = np.zeros_like( pan_array)
phi_array = np.zeros_like( pan_array)
r_array   = np.zeros_like( pan_array)

for i_panel in xrange( num_pan):
    panel_y = panels[i_panel] / 5
    panel_x = panels[i_panel] % 5 - 1
    ystart = pan_ind_y[ panel_y ]
    xstart = pan_ind_x[ panel_x ]
    
    y_array[i_panel] = Y[ ystart:ystart+ydim, xstart:xstart+xdim  ]
    x_array[i_panel] = X[ ystart:ystart+ydim, xstart:xstart+xdim  ]
    phi_array[i_panel] = PHI[ ystart:ystart+ydim, xstart:xstart+xdim  ]
    r_array[i_panel] = R[ ystart:ystart+ydim, xstart:xstart+xdim  ]

# find the intersection of the panel with desired the points along the ring
#panels_in_order = [27,22,23,24,29,34,39,38,37,32] # in order around the ring
#order_index = map( lambda x: panels.index( x), panels_in_order )
interp_here_array = []
for i in xrange( num_pan):
    if i == masked_1:
        interp_here = ring[ ring[:,1] > x_array[i].min() ]
        interp_here = interp_here[ interp_here[:,1] < x_array[i].max() ]
        interp_here = interp_here[ interp_here[:,0] > y_array[i,:-40].min() ]
        interp_here = interp_here[ interp_here[:,0] < y_array[i,:-40].max() ]
    elif i == masked_2:
        interp_here = ring[ ring[:,1] > x_array[i].min() ]
        interp_here = interp_here[ interp_here[:,1] < x_array[i].max() ]
        interp_here = interp_here[ interp_here[:,0] > y_array[i,15:].min() ]
        interp_here = interp_here[ interp_here[:,0] < y_array[i,15:].max() ]
    else:
        interp_here = ring[ ring[:,1] > x_array[i].min() ]
        interp_here = interp_here[ interp_here[:,1] < x_array[i].max() ]
        interp_here = interp_here[ interp_here[:,0] > y_array[i].min() ]
        interp_here = interp_here[ interp_here[:,0] < y_array[i].max() ]
    interp_here_array.append( interp_here )




pol_intens = np.ones(  num_phi   ) * -1
pol_intens_bin = np.ones(  num_phi  )* -1

figure(1, figsize=( 5,5) )
gs1 = GridSpec( 5,5 )

ax1 = subplot( gs1[:-2,:-2] ) # intensity ring
ax2 = subplot( gs1[:-2,-2:] )  # # hist
ax3 = subplot( gs1[-2,:-2] ) # binary
ax5 = subplot( gs1[-1,:-2] ) # ave binary
ax4 = subplot( gs1[-2:,-2:] ) # det image

subplots_adjust( left= 0.08, bottom=0.11,right=0.93, top=0.95, hspace=0.04,wspace=0.04 )

#   make the rings
intens_vals = []

NN = 0

for i in xrange( num_pan):
    ystart = pan_ind_y[ panels[i] /5 ]
    xstart = pan_ind_x[ panels[i] % 5 -1 ]
    
    pan_array = img[ ystart:ystart+ydim, xstart:xstart+xdim  ]
    pan_array = common_mode( pan_array)
    pan_array = depolarize( pan_array, r_array[i], phi_array[i], wavelen, pixsize,detdist, 0.99 )
    
    if i == masked_1:
        rbs = RBS ( y_array[i,:-40,0], x_array[i,0], pan_array[:-40])
    elif i == masked_2:
        rbs = RBS ( y_array[i,15:,0], x_array[i,0], pan_array[15:])
    else:
        rbs = RBS ( y_array[i,:,0], x_array[i,0], pan_array)

    phis = np.arctan2( interp_here_array[i][:,0 ]-b, interp_here_array[i][:,1]-a )
    vals = spline_eval( rbs, interp_here_array[i] ) 
#           make the ring binary

    NN += len( vals ) 

# the array values are the same down to floating point precision
    phi_inds = np.array( map( lambda x : find_nearest( phi_values, x ) , phis ) ) 
    pol_intens[   phi_inds] = np.copy(vals)
   
    cut_off = vals.mean() + vals.std()
    
    vals[ vals < cut_off  ] = 0
    vals[ vals > 0] = 1
    pol_intens_bin[  phi_inds] = vals
#       make it binary
    xm = np.where( vals == 1 )[0]
    ax1.plot( phi_inds, pol_intens[ phi_inds] ,'bx',ms=2.2 ) 
    ax1.plot( phi_inds[xm], pol_intens[phi_inds[xm]] , 'bx', ms=2.2)
    ax1.plot( [phi_inds[0],phi_inds[-1]] ,  np.ones( 2) * cut_off, c  ='red',lw=2.5 )
    
    ax3.plot( phi_inds, pol_intens_bin[ phi_inds] , 'b', lw = 0.5 ) 
    intens_vals += list(  pol_intens[ phi_inds ] )

Imin = 20500
Imax = 55000
# main
ax1.set_xlim(0,num_phi)
ax1.set_ylim(Imin,Imax)
ticks = [25000,30000,35000,40000,45000,50000]
ax1.get_yaxis().set_ticks(ticks)
ax1.get_yaxis().set_ticklabels( map( lambda x : str(x/1000)+'k' , ticks))
ax1.set_xlim(0,num_phi)
ax1.grid(True)
ax1.get_xaxis().set_ticks( [num_phi*90./360.,num_phi*180./360.,num_phi*270./360.] )
ax1.get_xaxis().set_ticklabels(['90','180','270'])
ax1.xaxis.tick_top()
ax1.text( num_phi * 10/360.,52000, "A", color='r', fontsize=18  )

#setp(ax1.get_xticklabels(),visible=False)


# binary
ax3.get_xaxis().set_ticks( [num_phi*90./360.,num_phi*180./360.,num_phi*270./360.] )
ax3.get_yaxis().set_ticks([1])
ax3.get_yaxis().set_ticklabels(['1'])
ax3.set_xlim(0,num_phi)
ax3.set_ylim(0,2.)
setp(ax3.get_xticklabels(),visible=False)
ax3.set_xlim(0,num_phi)
ax3.xaxis.grid(True)
ax3.text( num_phi * 10/360, 1.51, "C", color='r', fontsize=18  )

# hist
#ax2.get_yaxis().set_ticks([])
#setp(ax2.get_yticklabels(), visible=False)
ax2.hist( intens_vals, orientation='horizontal',bins=40,align='mid', log=True,rwidth=0.5,color='b', label='total pixels' )
ax2.set_ylim(Imin,Imax)
ax2.get_yaxis().set_ticks(ticks)
ax2.get_yaxis().set_ticklabels( map( lambda x : str(x/1000)+'k' , ticks))
ax2.get_xaxis().set_ticks([ 1,10,100 ])
ax2.get_xaxis().set_ticklabels( [r"10$^0$",r"10$^1$",r"10$^2$"] )
ax2.legend(prop={'size':8}, loc=(.35, .87))
ax2.xaxis.tick_top()
ax2.yaxis.tick_right()
ax2.get_xaxis().set_tick_params(which='both',width=1)
ax2.get_xaxis().set_tick_params(which='minor',length=4)
ax2.get_xaxis().set_tick_params(which='major',length=8)
ax2.yaxis.grid(True)
ax2.text(330, 47000,"B", color='r', fontsize=18  )

#ax2.get_xaxis().set_ticks([])

# image
img_masked = np.ma.masked_where( img2 < 0, img2 ) 
#cbar = colorbar(

ax4.imshow( img_masked , cmap='gist_yarg',vmin=-200,vmax=5000,aspect='auto')# , orientation='horizontal', 
 #   ticks = [600,4500],pad = 0.1,fraction=0.05)

#cbar.ax.set_xticklabels( ['low','high'] )
#cbar.ax.tick_params( which='both', size=0 )
more = 200
ax4.set_xlim(a-q-more,a+q + more)
ax4.set_ylim(b-q-more,b+q + more)
ax4.get_xaxis().set_ticks([])
ax4.get_yaxis().set_ticks([])
ax4.text( 1640,1510, "E", color='r', fontsize=18  )


# ave binary
scan = np.load( '/mnt/ssrl/silver/12dec2012/scan.npz')['scan']
ax5.plot( scan.mean(0),'bx',ms=1.5)
#ax5.get_yaxis().set_ticks([])

ax5.get_yaxis().set_ticks([0.2])
ax5.get_yaxis().set_ticklabels(['0.2'])
ax5.get_xaxis().set_ticks([num_phi*90./360.,num_phi*180./360.,num_phi*270./360.])
ax5.get_xaxis().set_ticklabels(['90','180','270'])
ax5.set_ylim(0,.4)
ax5.set_xlim(0,num_phi)
ax5.set_xlabel(r'$\phi$ ( deg.)',fontsize=18)
ax5.xaxis.grid(True)
ax5.text( num_phi * 10/360.,0.295, "D", color='r', fontsize=18  )

plt.show()
