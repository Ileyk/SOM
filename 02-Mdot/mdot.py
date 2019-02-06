import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

import matplotlib.cm as cm
from matplotlib.colors import rgb2hex
# from bokeh import mpl
# from bokeh.plotting import output_file, show, ColumnDataSource, figure, vplot
# from bokeh.models import HoverTool
import pandas as pd
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import FormatStrFormatter
import matplotlib.ticker as ticker
import matplotlib.image as mpimg
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatterMathtext
import holoviews as hv
# import colorcet as cc
from matplotlib.colors import LinearSegmentedColormap

from matplotlib.patches import Rectangle

import os

import subprocess

from parfile import *

# display
# (M2,P)
# (3,4) | (3,8) | (3,12)
# (2,4) | (2,8) | (2,12)
# (1,4) | (1,8) | (1,12)
# => M2[i//3] and P[i%3]
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

hv.notebook_extension("matplotlib")

# x = np.arange(0, np.pi, 0.1)
# y = np.arange(0, 2*np.pi, 0.1)
# X, Y = np.meshgrid(x, y)
# Z = np.cos(X) * np.sin(Y) * 10
# colors = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1)]  # R -> G -> B
# n_bins = [4, 6, 10, 100]  # Discretizes the interpolation into bins
# cmap_name = 'my_list'
# # fig, ax = plt.subplots(1, 1, figsize=(6, 6))
# vmax=25
# cm = LinearSegmentedColormap.from_list(cmap_name, [(0,    'blue'),
#                                               (0.2/vmax, 'red'),
#                                               (0.6/vmax, 'green'),
# 											  (2./vmax, 'black'),
# 											  (6./vmax, 'yellow'),
# 											  (20./vmax, 'orange'),
# 											  (vmax/vmax, 'cyan')], N=7)

# im = ax.imshow(Z, interpolation='nearest', origin='lower', cmap=cm)
# ax.set_title("N bins: %s" % 4)
# fig.colorbar(im, ax=ax)
# plt.show()

# fig.subplots_adjust(left=0.02, bottom=0.06, right=0.95, top=0.94, wspace=0.05)
# for ax in zip(axs.ravel()):
#     # Create the colormap
#     cm = LinearSegmentedColormap.from_list(
#         cmap_name, colors, N=4)
#     # Fewer bins will result in "coarser" colomap interpolation
#     im = ax.imshow(Z, interpolation='nearest', origin='lower', cmap=cm)
#     ax.set_title("N bins: %s" % 4)
#     fig.colorbar(im, ax=ax)
# plt.show()

# def plot_array():

def fmt(x, pos):
	a, b = '{:.1e}'.format(x).split('e')
	b = int(b)
	return r'{} $\cdot $10$^{{{}}}$'.format(a, b)

# def fmtCbar(x, pos):
# 	#a = '{:5.2f}'.format(x) #.split('e')[0]
# 	a = '{:5.1f}'.format(x)#.split('e')[0]
# 	return a
# 	# return r'{}'.format(a)

def fmtCbar(x, pos):
	a, b = '{:.1e}'.format(x).split('e')
	b = int(b)
	if (b!=0 and b!=1 and b!=-1):
		return r'{}E{}'.format(a, b)
	elif (b==0):
		return r'{}'.format(a)
	elif (b==1):
		return '{:2.0f}'.format(x)
	elif (b==-1):
		return '{:3.2f}'.format(x)

font = {'family' : 'normal',
'size'   : fontsize}

#'weight' : 'bold',

plt.rc('font', **font)

om_orb_max=2.
om_orb_min=0.001
N=100
om_orb=np.linspace(om_orb_min,om_orb_max,N)
Jtot=np.zeros(N)
Jtot=0.25*om_orb*(1.+3./(om_orb**(4./3.)))
fig, fig1 = plt.subplots(1,figsize=(1.1*figsize,figsize))
fig1.plot(om_orb,Jtot,linestyle='solid',marker='o',linewidth=3,color=colors[0])
fig1.set_xlim(om_orb_min,om_orb_max)
fig1.set_xlabel(r'$\Omega_{orb}$',fontweight='bold',fontsize=fontsize)
fig1.set_ylabel(r'Total angular momentum', fontweight='bold', fontsize=fontsize)
fig1.grid(which='major', linestyle='dotted', linewidth=2, alpha=0.9)
fig1.grid(which='minor', linestyle='dotted', linewidth=2, alpha=0.9)
fig.tight_layout()
fig.savefig(outputs+'total_angular_momentum.png',bbox_inches='tight')
plt.close()
# - - -
M2=2.
om_0_to_om_orb_0=1.
M1=np.asarray([35.,25.,15.])
om=np.zeros((len(M1),N))
M10=np.asarray(M1[0])
al=np.zeros(len(M1))
C=np.zeros(len(M1))
al=3.*((M10+M2)/(M1+M2))**(2./3.)
C=al*(om_0_to_om_orb_0)**(1./3.)+1.
fig, fig1 = plt.subplots(1,figsize=(1.1*figsize,figsize))
for i in range(len(M1)):
	print M1[i], al[i], C[i]
	om[i][:]=-al[i]*om_orb[:]**(-1./3.)+C[i]
	fig1.plot(om_orb,om[i][:],linestyle='solid',marker='o',linewidth=3,color=colors[i])
fig1.plot(om_orb,om_orb,linestyle='dashed',linewidth=3,color='k')
fig1.plot(1./om_0_to_om_orb_0,1.,marker='o',markersize=10,linewidth=0,color='k')
fig1.text(1./om_0_to_om_orb_0,1.1,r'Starting point', fontsize=fontsize/2, color='k', fontweight='bold', rotation=0, horizontalalignment='center',verticalalignment='center')
fig1.text(0.8, 0.2, r'only d.o.f. is om_0_to_om_orb_0', fontsize=fontsize/2, color='k', fontweight='bold', rotation=0, horizontalalignment='center',verticalalignment='center',transform = fig1.transAxes)
fig1.set_xlim(om_orb_min,om_orb_max)
fig1.set_ylim(0.,np.max(C)/2.)
fig1.text(0.2, 0.9, r'Red is initial mass (35), then green (25) then blue (15)', fontsize=fontsize/2, color='k', fontweight='bold', rotation=0, horizontalalignment='center',verticalalignment='center',transform = fig1.transAxes)
fig1.set_xlabel(r'$\Omega_{orb}$',fontweight='bold',fontsize=fontsize)
fig1.set_ylabel(r'$\omega$', fontweight='bold', fontsize=fontsize)
fig1.grid(which='major', linestyle='dotted', linewidth=2, alpha=0.9)
fig1.grid(which='minor', linestyle='dotted', linewidth=2, alpha=0.9)
fig.tight_layout()
fig.savefig(outputs+'om_as_a_function_of_om_orb_and_M1.png',bbox_inches='tight')
plt.close()
# - - -
M2=2.
M1_min=10.
M1_max=35. # inital mass
M1=np.linspace(M1_min,M1_max,N)
M10=M1[len(M1)-1]
om_c=((M1+M2)/(M10+M2))**(-0.25) # ratio of critial value to initial critical value
fig, fig1 = plt.subplots(1,figsize=(1.1*figsize,figsize))
fig1.plot(M1,om_c,linestyle='solid',marker='o',linewidth=3,color='k')
fig1.set_xlim(M1_max,M1_min)
fig1.set_ylim(0.,1.2*np.max(om_c))
fig1.set_xlabel(r'$M_1$',fontweight='bold',fontsize=fontsize)
fig1.set_ylabel(r'$\Omega_{c}$', fontweight='bold', fontsize=fontsize)
fig1.text(0.4, 0.2, r'As time goes by, mass decreases (so we move to right)', fontsize=fontsize/2, color='k', fontweight='bold', rotation=0, horizontalalignment='center',verticalalignment='center',transform = fig1.transAxes)
fig1.text(0.4, 0.1, r'Conclusion : not a significant increase', fontsize=fontsize/2, color='k', fontweight='bold', rotation=0, horizontalalignment='center',verticalalignment='center',transform = fig1.transAxes)
fig1.grid(which='major', linestyle='dotted', linewidth=2, alpha=0.9)
fig1.grid(which='minor', linestyle='dotted', linewidth=2, alpha=0.9)
fig.tight_layout()
fig.savefig(outputs+'synch_shift_as_mass_is_lost.png',bbox_inches='tight')
