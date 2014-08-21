#generic python modules
import argparse
import operator
from operator import itemgetter
import sys, os, shutil
import os.path

##########################################################################################
# RETRIEVE USER INPUTS
##########################################################################################

#=========================================================================================
# create parser
#=========================================================================================
version_nb = "0.0.1"
parser = argparse.ArgumentParser(prog = 'xvg_plot_energy_colour', usage='', add_help = False, formatter_class = argparse.RawDescriptionHelpFormatter, description =\
'''
*****************************************************
v''' + version_nb + '''
author: Jean Helie (jean.helie@bioch.ox.ac.uk)
git: https://github.com/jhelie/xvg_plot_energy_colour
*****************************************************

[ DESCRIPTION ]

This script plots the evolution of the potential energy and colours it
based on the conformation status of that peptide

[ USAGE ]

Option	      Default  	Description                    
-----------------------------------------------------
-e			: xvg file for energy evolution
-c			: xvg file(s) for conformation status
-o	epot_vs_status	: name of outptut file
--micro			: use microsecond instead of ns for x axis
--comments	@,#	: lines starting with these characters will be considered as comment

Other options
-----------------------------------------------------
--version		: show version number and exit
-h, --help		: show this menu and exit
 
''')

#options
parser.add_argument('-e', nargs=1, dest='energy_xvgfilename', help=argparse.SUPPRESS, required=True)
parser.add_argument('-c', nargs=1, dest='status_xvgfilename', help=argparse.SUPPRESS, required=True)
parser.add_argument('-o', nargs=1, dest='output_file', default=["epot_vs_status"], help=argparse.SUPPRESS)
parser.add_argument('--micro', dest='micro', action='store_true', help=argparse.SUPPRESS)
parser.add_argument('--comments', nargs=1, dest='comments', default=['@,#'], help=argparse.SUPPRESS)

#other options
parser.add_argument('--version', action='version', version='%(prog)s v' + version_nb, help=argparse.SUPPRESS)
parser.add_argument('-h','--help', action='help', help=argparse.SUPPRESS)

#=========================================================================================
# store inputs
#=========================================================================================

args = parser.parse_args()
args.energy_xvgfilename = args.energy_xvgfilename[0]
args.status_xvgfilename = args.status_xvgfilename[0]
args.output_file = args.output_file[0]
args.comments = args.comments[0].split(',')

#=========================================================================================
# import modules (doing it now otherwise might crash before we can display the help menu!)
#=========================================================================================

#generic science modules
try:
	import numpy as np
except:
	print "Error: you need to install the numpy module."
	sys.exit(1)
try:
	import matplotlib as mpl
	mpl.use('Agg')
	import matplotlib.colors as mcolors
	mcolorconv = mcolors.ColorConverter()
	import matplotlib.cm as cm				#colours library
	import matplotlib.ticker
	from matplotlib.ticker import MaxNLocator
	from matplotlib.font_manager import FontProperties
	fontP=FontProperties()
	from matplotlib.collections import LineCollection

except:
	print "Error: you need to install the matplotlib module."
	sys.exit(1)
try:
	import pylab as plt
except:
	print "Error: you need to install the pylab module."
	sys.exit(1)

#=======================================================================
# sanity check
#=======================================================================

if not os.path.isfile(args.energy_xvgfilename):
	print "Error: file " + str(args.energy_xvgfilename) + " not found."
	sys.exit(1)

if not os.path.isfile(args.status_xvgfilename):
	print "Error: file " + str(args.status_xvgfilename) + " not found."
	sys.exit(1)

##########################################################################################
# FUNCTIONS DEFINITIONS
##########################################################################################

#=========================================================================================
# data loading
#=========================================================================================

def load_xvg():															#DONE
	
	global y_min
	global y_max
	global times
	global data_energy
	global data_status
	
	#energy
	#------
	#get file content
	print " -reading energy file... "
	filename = args.energy_xvgfilename
	with open(filename) as f:
		lines = f.readlines()
	
	#determine legends and nb of lines to skip
	tmp_nb_rows_to_skip = 0
	for l_index in range(0,len(lines)):		
		line = lines[l_index]
		if line[0] in args.comments:
			tmp_nb_rows_to_skip += 1
		
	#get data
	tmp_data = np.loadtxt(filename, skiprows = tmp_nb_rows_to_skip)
	nb_rows = np.shape(tmp_data)[0]
	nb_cols = np.shape(tmp_data)[1]
	if nb_cols > 2:
		print "Error: more than 2 columns detected in file " + str(filename)
		sys.exit(1)
	
	#store time
	times = np.zeros((nb_rows,1))
	times[:,0] = tmp_data[:,0]/float(1000)

	#store energy
	data_energy = np.zeros((nb_rows,1))
	data_energy[:,0] = tmp_data[:,1]
	
	#status
	#------
	#get file content
	print " -reading conformation status file..."
	filename = args.status_xvgfilename
	with open(filename) as f:
		lines = f.readlines()
	
	#determine legends and nb of lines to skip
	tmp_nb_rows_to_skip = 0
	for l_index in range(0,len(lines)):
		line = lines[l_index]
		if line[0] in args.comments:
			tmp_nb_rows_to_skip += 1
	
	#get data
	tmp_data = np.loadtxt(filename, skiprows = tmp_nb_rows_to_skip)
	nb_rows = np.shape(tmp_data)[0]
	nb_cols = np.shape(tmp_data)[1]
	if np.shape(tmp_data)[0] != nb_rows:
		print "Error: file " + str(filename) + " has " + str(np.shape(tmp_data)[0]) + " data rows, whereas file " + str(args.energy_xvgfilename) + " has " + str(nb_rows) + " data rows."
		sys.exit(1)
	if nb_cols > 2:
		print "Error: more than 2 columns detected in file " + str(filename)
		sys.exit(1)
	if not np.array_equal(tmp_data[:,0],times[:,0]):
		print "\nError: the first column of file " + str(filename) + " is different than that of " + str(args.energy_xvgfilename) + "."
		sys.exit(1)
	
	#store data
	data_status = np.zeros((nb_rows,1))
	data_status[:,0] = tmp_data[:,1]
	
	#switch to microseconds
	if args.micro:
		times[:,0] = times[:,0]/float(1000)
	
	return

#=========================================================================================
# outputs
#=========================================================================================


# Create a colormap for red, green and blue and a norm to color
# f' < -0.5 red, f' > 0.5 blue, and the rest green

#  - 0: surfacic
#  - 1: TM
#  - 2: TM*
#  - 3: U-shape

def graph_xvg():
	
	print " -plotting data..."
	
	#create line collections
	status = [0,1,2,3]
	labels = ['interfacial','TM','TM*','U-shape']
	colours = ['b','y','g','m']
	cmap = mcolors.ListedColormap(colours)
	norm = mcolors.BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], cmap.N)
	points = np.array([times[:,0],data_energy[:,0]]).T.reshape(-1,1,2)
	segments = np.concatenate([points[:-1], points[1:]], axis=1)	
	lc = LineCollection(segments, cmap = cmap, norm = norm)
	lc.set_array(data_status[:,0])
	lc.set_linewidth(3)
	
	#open files
	filename_svg = os.getcwd() + '/' + str(args.output_file) + '.svg'
	
	#create figure
	fig = plt.figure(figsize=(8, 6.2))
	fig.suptitle("System potential energy")

	#plot line
	ax = plt.gca()
	#ax.add_collection(lc)
	plt.gca().add_collection(lc)
		
	if args.micro:
		plt.xlabel('time (us)')	
	else:
		plt.xlabel('time (ns)')
	plt.ylabel('potential energy (kJ.mol-1)')
	
	#save figure
	plt.xlim(times[:,0].min(), times[:,0].max())
	#ax.set_xlim(0, 5000)
	ax.set_ylim(data_energy[:,0].min(), -91000)
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.xaxis.set_ticks_position('bottom')
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_major_locator(MaxNLocator(nbins=6))
	ax.yaxis.set_major_locator(MaxNLocator(nbins=10))
	ax.xaxis.labelpad = 20
	ax.yaxis.labelpad = 20
	ax.get_xaxis().set_tick_params(direction='out')
	ax.get_yaxis().set_tick_params(direction='out')
	plt.setp(ax.xaxis.get_majorticklabels(), fontsize = "small")
	plt.setp(ax.yaxis.get_majorticklabels(), fontsize = "small")
	fig.savefig(filename_svg, transparent = True)
	plt.close()

	return

##########################################################################################
# MAIN
##########################################################################################

load_xvg()
graph_xvg()

#=========================================================================================
# exit
#=========================================================================================
print "\nFinished successfully! Check result in file '" + str(args.output_file) + ".svg'."
print ""
sys.exit(0)
