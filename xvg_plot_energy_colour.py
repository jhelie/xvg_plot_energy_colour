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
version_nb = "0.0.2"
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

The energy file can contain 2 or 3 columns (in the latter case the last column is
interpreted as a standard deviation) and the conformation status file can contain any
number of columns but only the first 2 will be considered.

The first column (correspondig to time) must be exactly the same in each file

[ USAGE ]

Option	      Default  	Description                    
-----------------------------------------------------
-e			: xvg file for energy evolution
-c			: xvg file(s) for conformation status
-o	epot_vs_status	: name of outptut file
--preset		: use pre-set colours for colours (default = use jet scale based on min-max values in -c file)
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
parser.add_argument('--preset', dest='preset', action='store_true', help=argparse.SUPPRESS)
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
	global data_energy_avg
	global data_energy_max
	global data_energy_min
	global data_status
	global nb_cols
	
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
	tmp_data = tmp_data[1:,:]		#remove first line
	nb_rows = np.shape(tmp_data)[0]
	nb_cols = np.shape(tmp_data)[1]
	if nb_cols > 3:
		print "Error: more than 3 columns detected in file " + str(filename)
		sys.exit(1)
	
	#store time
	times = np.zeros((nb_rows,1))
	times[:,0] = tmp_data[:,0]/float(1000)

	#store energy
	data_energy_avg = np.zeros((nb_rows,1))
	data_energy_avg[:,0] = tmp_data[:,1]
	if nb_cols > 2:
		data_energy_max = np.zeros((nb_rows,1))
		data_energy_min = np.zeros((nb_rows,1))
		data_energy_max[:,0] = data_energy_avg[:,0] + tmp_data[:,2]
		data_energy_min[:,0] = data_energy_avg[:,0] - tmp_data[:,2]
	
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
	tmp_data = tmp_data[1:,:]		#remove first line
	if np.shape(tmp_data)[0] != nb_rows:
		print "Error: file " + str(filename) + " has " + str(np.shape(tmp_data)[0]) + " data rows, whereas file " + str(args.energy_xvgfilename) + " has " + str(nb_rows) + " data rows."
		sys.exit(1)
	if np.shape(tmp_data)[1] > 2:
		print "Warning: more than 2 columns detected in file " + str(filename) + ", only the first 2 will be takin into account."
	if not np.array_equal(tmp_data[:,0],times[:,0]):
		print "\nError: the first column of file " + str(filename) + " is different than that of " + str(args.energy_xvgfilename) + "."
		sys.exit(1)
	
	#store data
	data_status = np.zeros((nb_rows,1))
	data_status[:,0] = tmp_data[:,1]
	
	#post-processing
	#---------------
	#switch to microseconds
	if args.micro:
		times[:,0] = times[:,0]/float(1000)
	
	#shift data down
	y_min = 0
	if nb_cols > 2:
		tmp_min = np.min(data_energy_min[:,0])
		data_energy_avg[:,0] -= tmp_min
		data_energy_max[:,0] -= tmp_min
		data_energy_min[:,0] -= tmp_min	
		y_max = 1.1*np.max(data_energy_max[:,0])
		y_min = 0.9*np.min(data_energy_min[:,0])
	else:
		tmp_min = np.min(data_energy_avg[:,0])
		data_energy_avg[:,0] -= tmp_min
		y_max = 1.1*np.max(data_energy_avg[:,0])
		y_min = 0.9*np.min(data_energy_avg[:,0])
	
	#y_max = 850
	
	return

#=========================================================================================
# outputs
#=========================================================================================


def graph_xvg():
	
	print " -plotting data..."
	
	#create line collections
	if args.preset:
		# Create a colormap for red, green and blue and a norm to color
		# f' < -0.5 red, f' > 0.5 blue, and the rest green
		
		#  - 0: surfacic
		#  - 1: TM
		#  - 2: TM*
		#  - 3: U-shape

		status = [0,1,2,3]
		labels = ['interfacial','TM','TM*','U-shape']
		colours = ['b','y','g','m']
		cmap = mcolors.ListedColormap(colours)
		norm = mcolors.BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], cmap.N)
	else:

		#min / max
		#status_min = int(np.min(data_status[:,0]))
		#status_max = int(np.max(data_status[:,0]))
		status_min = 1
		status_max = 9
		
		#labels
		labels = []
		for c_size in range(status_min, status_max+1):
			labels.append(str(c_size))

		#boundaries (assumes everyone is > 0)
		tmp_boundaries = np.arange(status_min-1, status_max+1)
		tmp_boundaries = tmp_boundaries + 0.5
		
		#colours
		colours = []
		tmp_cmap = cm.get_cmap('jet')
		colours_sizes_value = tmp_cmap(np.linspace(0, status_min, status_max-status_min+1))
		for c_index in range(0, status_max-status_min+1):
			colours.append(colours_sizes_value[c_index])
		cmap = mcolors.ListedColormap(colours)
		norm = mcolors.BoundaryNorm(tmp_boundaries, cmap.N)
		
	points = np.array([times[:,0],data_energy_avg[:,0]]).T.reshape(-1,1,2)
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
	ax.add_collection(lc)
	if nb_cols > 2:
		plt.fill_between(times[:,0], data_energy_min[:,0], data_energy_max[:,0], color = "#262626", edgecolor = "#262626", linewidth = 0, alpha = 0.2)		
	
	if args.micro:
		plt.xlabel('time (us)')	
	else:
		plt.xlabel('time (ns)')
	plt.ylabel('relative potential energy (kJ.mol-1)')
	
	#save figure
	ax.set_xlim(times[:,0].min(), times[:,0].max())
	ax.set_ylim(y_min, y_max)
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.xaxis.set_ticks_position('bottom')
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_major_locator(MaxNLocator(nbins=6))
	ax.yaxis.set_major_locator(MaxNLocator(nbins=10))
	#ax.xaxis.labelpad = 20
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
