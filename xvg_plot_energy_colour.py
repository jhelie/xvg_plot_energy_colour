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

if len(args.particles_xvgfilenames) != len(args.charges_xvgfilenames):
	print "Error: different number of particles files (" + str(len(args.particles_xvgfilenames)) + ") and charges files (" + str(len(args.charges_xvgfilenames)) + ")."
	sys.exit(1)
	
for f in args.particles_xvgfilenames + args.charges_xvgfilenames:
	if not os.path.isfile(f):
		print "Error: file " + str(f) + " not found."
		sys.exit(1)

##########################################################################################
# FUNCTIONS DEFINITIONS
##########################################################################################

labels = ['AM Cter outwards','AM Nter outwards','SMa']
#labels = ['r_1','r_2','r_3','r_1','r_2','r_3']
colours = ['#579D1C','#579D1C','#579D1C','#7E0021','#7E0021','#7E0021']
#colours = ['#1D91C0','#579D1C','#7E0021']		#cyan, green, dark red

#=========================================================================================
# data loading
#=========================================================================================

def load_xvg():															#DONE
	
	global y_min
	global y_max
	global data_density_pops
	global data_density_charge
	global boundaries_charge_min
	global boundaries_charge_max
	global boundaries_pops_min
	global boundaries_pops_max
	tmp_min = np.float("inf")
	
	#particles
	#---------
	for f_index in range(0,len(args.particles_xvgfilenames)):
		#display progress
		progress = '\r -reading particle file ' + str(f_index+1) + '/' + str(len(args.particles_xvgfilenames)) + '                      '  
		sys.stdout.flush()
		sys.stdout.write(progress)
		
		#get file content
		filename = args.particles_xvgfilenames[f_index]
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
		
		#check that each file has the same number of data rows
		if f_index == 0:
			nb_rows = np.shape(tmp_data)[0]
			distances = np.zeros((nb_rows,1))													#distance from cluster
			data_density_pops = np.zeros(len(args.particles_xvgfilenames))
		else:
			if np.shape(tmp_data)[0] != nb_rows:
				print "Error: file " + str(filename) + " has " + str(np.shape(tmp_data)[0]) + " data rows, whereas file " + str(args.xvgfilenames[0]) + " has " + str(nb_rows) + " data rows."
				sys.exit(1)
		#check that each file has the same number of columns
		if f_index == 1:
			nb_cols = np.shape(tmp_data)[1]
		elif f_index ==2:
			if np.shape(tmp_data)[1] != nb_cols:
				print "Error: file " + str(filename) + " has " + str(np.shape(tmp_data)[1]) + " data columns, whereas file " + str(args.xvgfilenames[0]) + " has " + str(nb_cols) + " data columns."
				sys.exit(1)
		#check that each file has the same first column
		if f_index == 0:
			distances[:,0] = tmp_data[:,0]
		else:
			if not np.array_equal(tmp_data[:,0],distances[:,0]):
				print "\nError: the first column of file " + str(filename) + " is different than that of " + str(args.xvgfilenames[0]) + "."
				sys.exit(1)
		
		#store data
		if args.density:
			data_density_pops[f_index] = tmp_data[80,4]
		else:
			data_density_pops[f_index] = -np.log(tmp_data[80,4])
			tmp_min = min([tmp_min, np.nanmin(-np.log(tmp_data[:,4]))])

	if args.density:
		#y_max = max(data_density_pops)*1.1
		y_max = 0.00130
	else:
		data_density_pops -= tmp_min
		y_max = 12
	
	print ''
	
	#charges
	#-------
	for f_index in range(0,len(args.charges_xvgfilenames)):
		#display progress
		progress = '\r -reading charge file ' + str(f_index+1) + '/' + str(len(args.charges_xvgfilenames)) + '                      '  
		sys.stdout.flush()
		sys.stdout.write(progress)
		
		#get file content
		filename = args.charges_xvgfilenames[f_index]
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
		
		#check that each file has the same number of data rows
		if f_index == 0:
			nb_rows = np.shape(tmp_data)[0]
			distances = np.zeros((nb_rows,1))													#distance from cluster
			data_density_charge = np.zeros(len(args.charges_xvgfilenames))
		else:
			if np.shape(tmp_data)[0] != nb_rows:
				print "Error: file " + str(filename) + " has " + str(np.shape(tmp_data)[0]) + " data rows, whereas file " + str(args.charges_xvgfilenames[0]) + " has " + str(nb_rows) + " data rows."
				sys.exit(1)
		#check that each file has the same number of columns
		if f_index == 1:
			nb_cols = np.shape(tmp_data)[1]
		elif f_index ==2:
			if np.shape(tmp_data)[1] != nb_cols:
				print "Error: file " + str(filename) + " has " + str(np.shape(tmp_data)[1]) + " data columns, whereas file " + str(args.charges_xvgfilenames[0]) + " has " + str(nb_cols) + " data columns."
				sys.exit(1)
		#check that each file has the same first column
		if f_index == 0:
			distances[:,0] = tmp_data[:,0]
		else:
			if not np.array_equal(tmp_data[:,0],distances[:,0]):
				print "\nError: the first column of file " + str(filename) + " is different than that of " + str(args.charges_xvgfilenames[0]) + "."
				sys.exit(1)
		
		#store data
		data_density_charge[f_index] = tmp_data[80,4]
	
	#find x and y boundaries
	#-----------------------
	if args.density:
		boundaries_charge_min = 1.01 * max([max(data_density_charge[0:3]),data_density_charge[5],max(data_density_charge[6:9]),max(data_density_charge[12:15])])  #take into account run #3 of AM_zCter as there's no flipflop there
		boundaries_charge_max = 0.99 * min([min(data_density_charge[3:5]),min(data_density_charge[9:12]),min(data_density_charge[15:17])]) #don't take into account run #3 of AM_zCter as there's no flipflop there
		boundaries_pops_min = 1.01 * max([max(data_density_pops[0:3]),data_density_pops[5],max(data_density_pops[6:9]),max(data_density_pops[12:15])])  #take into account run #3 of AM_zCter as there's no flipflop there
		boundaries_pops_max = 0.99 * min([min(data_density_pops[3:5]),min(data_density_pops[9:12]),min(data_density_pops[15:17])]) #don't take into account run #3 of AM_zCter as there's no flipflop there	
	else:
		boundaries_charge_min = 1.01 * max([max(data_density_charge[0:3]),data_density_charge[5],max(data_density_charge[6:9]),max(data_density_charge[12:15])])  #take into account run #3 of AM_zCter as there's no flipflop there
		boundaries_charge_max = 0.99 * min([min(data_density_charge[3:5]),min(data_density_charge[9:12]),min(data_density_charge[15:17])]) #don't take into account run #3 of AM_zCter as there's no flipflop there
		boundaries_pops_max = 0.99 * min([min(data_density_pops[0:3]),data_density_pops[5],min(data_density_pops[6:9]),min(data_density_pops[12:15])])  #take into account run #3 of AM_zCter as there's no flipflop there
		boundaries_pops_min = 1.01 * max([max(data_density_pops[3:5]),max(data_density_pops[9:12]),max(data_density_pops[15:17])]) #don't take into account run #3 of AM_zCter as there's no flipflop there
	
	return

#=========================================================================================
# outputs
#=========================================================================================

def graph_xvg():

	#open files
	filename_svg = os.getcwd() + '/' + str(args.output_file) + '.svg'
	
	#create figure
	fig = plt.figure(figsize=(8, 6.2))
	if args.density:
		fig.suptitle("POPS density versus charge density at z = 0")
	else:
		fig.suptitle("POPS free energy versus charge density at z = 0")

	#plot data
	ax = fig.add_subplot(111)
	#plot AM_zCter
	plt.plot(data_density_charge[0], data_density_pops[0], color = colours[0], label = labels[0], linewidth = 0, marker = 'v', markersize=9)
	for f_index in range(1,6):
		plt.plot(data_density_charge[f_index], data_density_pops[f_index], color = colours[f_index], linewidth = 3, marker = 'v', markersize=9)
	#plot AM_zNter
	plt.plot(data_density_charge[6], data_density_pops[6], color = colours[0], label = labels[1], linewidth = 0, marker = '^', markersize=9)
	for f_index in range(7,12):
		plt.plot(data_density_charge[f_index], data_density_pops[f_index], color = colours[f_index-6], linewidth = 3, marker = '^', markersize=9)
	#plot SMa
	plt.plot(data_density_charge[12], data_density_pops[12], color = colours[0], label = labels[2], linewidth = 0, marker = 'o', markersize=8)
	for f_index in range(13,18):
		plt.plot(data_density_charge[f_index], data_density_pops[f_index], color = colours[f_index-12], linewidth = 3, marker = 'o', markersize=8)

	#plot boundaries
	ax.axhspan(boundaries_pops_min, boundaries_pops_max, xmin = 0, xmax = boundaries_charge_min/float(0.001), alpha=0.2, color="#262626", linewidth = 0) #, edgecolor = "#262626", linewidth = 0)
	ax.axhspan(boundaries_pops_min, boundaries_pops_max, xmin = boundaries_charge_max/float(0.001), xmax = 1, alpha=0.2, color="#262626", linewidth = 0) #, edgecolor = "#262626", linewidth = 0)
	ax.axvspan(boundaries_charge_min, boundaries_charge_max, alpha=0.2, color="#262626", linewidth = 0)
		
	fontP.set_size("small")
	ax.legend(prop = fontP, numpoints = 1)
	plt.xlabel('charge density [$e.\AA^{-3}$]')
	if args.density:
		plt.ylabel('POPS density [$\AA^{-3}$]')
	else:
		plt.ylabel('free energy [$kT$]')
	
	#save figure
	ax.set_xlim(0, 0.0010)
	ax.set_ylim(0, y_max)
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.xaxis.set_ticks_position('bottom')
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_major_locator(MaxNLocator(nbins=11))
	ax.yaxis.set_major_locator(MaxNLocator(nbins=12))
	ax.xaxis.labelpad = 20
	ax.yaxis.labelpad = 20
	ax.get_xaxis().set_tick_params(direction='out')
	ax.get_yaxis().set_tick_params(direction='out')
	plt.setp(ax.xaxis.get_majorticklabels(), fontsize = "small")
	plt.setp(ax.yaxis.get_majorticklabels(), fontsize = "small")
	plt.subplots_adjust(top = 0.9, bottom = 0.15, left = 0.15, right = 0.85)
	fig.savefig(filename_svg, transparent = True)
	plt.close()

	return
def write_xvg():
	#open files
	filename_xvg = os.getcwd() + '/' + str(args.output_file) + '.xvg'
	output_xvg = open(filename_xvg, 'w')
	
	#general header
	output_xvg.write("# [average xvg - written by xvg_plot_density_particle_vs_charge v" + str(version_nb) + "]\n")
	tmp_files = ""
	for f in args.particles_xvgfilenames:
		tmp_files += "," + str(f)
	output_xvg.write("# - particles files: " + str(tmp_files[1:]) + "\n")
	tmp_files = ""
	for f in args.charges_xvgfilenames:
		tmp_files += "," + str(f)
	output_xvg.write("# - particles files: " + str(tmp_files[1:]) + "\n")
	output_xvg.write("# - boundaries charge min: " + str(boundaries_charge_min) + "\n")
	output_xvg.write("# - boundaries charge max: " + str(boundaries_charge_max) + "\n")
	output_xvg.write("# - boundaries pops min: " + str(boundaries_pops_min) + "\n")
	output_xvg.write("# - boundaries pops max: " + str(boundaries_pops_max) + "\n")
	
	#xvg metadata
	output_xvg.write("@ title \"particles density vs charges density at z = 0 xvg\"\n")
	output_xvg.write("@ xaxis label \"charges density (e.Angstrom-3)\"\n")
	output_xvg.write("@ yaxis label \"particles density (Angstrom-3)\"\n")
	output_xvg.write("@ autoscale ONREAD xaxes\n")
	output_xvg.write("@ TYPE XY\n")
	output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
	output_xvg.write("@ legend on\n")
	output_xvg.write("@ legend box on\n")
	output_xvg.write("@ legend loctype view\n")
	output_xvg.write("@ legend 0.98, 0.8\n")
	output_xvg.write("@ legend length " + str(2*len(args.particles_xvgfilenames)) + "\n")
	for f_index in range(0,len(args.particles_xvgfilenames)):
		output_xvg.write("@ s" + str(f_index) + " legend \""+ str(args.particles_xvgfilenames[f_index]) +"\"\n")
	for f_index in range(0,len(args.charges_xvgfilenames)):
		output_xvg.write("@ s" + str(f_index) + " legend \""+ str(args.charges_xvgfilenames[f_index]) +"\"\n")
	
	#replace inf by nan (since they're not plotted)
	data_density_pops[np.isinf(data_density_pops)] = np.nan
	data_density_charge[np.isinf(data_density_charge)] = np.nan
	
	#data
	for f_index in range(0, len(args.particles_xvgfilenames)):
		results = str(data_density_pops[f_index]) + "	" + str(data_density_charge[f_index])
		output_xvg.write(results + "\n")		
	output_xvg.close()	

	return

##########################################################################################
# MAIN
##########################################################################################

load_xvg()
graph_xvg()
write_xvg()

#=========================================================================================
# exit
#=========================================================================================
print "\nFinished successfully! Check result in file '" + str(args.output_file) + "'."
print ""
sys.exit(0)
