#
# main.py
#
# Python program to execute read and analysis of simulation files.
# See script headers for more details on code I/O.
#
# Author: Benjamin MacFarlane
# Date: 27/07/2016
# Contact: bmacfarlane@uclan.ac.uk
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - VARIABLE DEFINITIONS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
arch_dir = "/Users/bmacfarlane1992/Documents/ACADEMIA_PHD/PROJECTS/FEEDBACK-STRUCTURE/"	# Location of project directory architecture
dat_dir = "/Volumes/SD_bm1992/DATA/"			# Location of data from analyse_disc.f90
ic_dir = "/Volumes/SD_bm1992/ICs/"				# Location of ICs to read in sink data
plot_form = "eps"								# Format of plots (either "eps" or "png")
col_arr = ["b", "g", "r", "c", "m", "k"]		# Array of color pointers for plots throughout
#
ea_run = [0, 1, 3, 4, 5]		# Select EA runs to process
#
r_start = 10						# Starting radius for analyses of pdisc and ea_comp
r_limit = 150					# Limit of radial plots in pdisc analyses
smooth = 3						# Value of AU that radial profiles are smoothed over
#
#
EA_timeref = ["BEFORE","POST-ONSET","PRE-OFFSET","AFTER"]		# Define names of EA snapshots for EA length reference
EA_lenref = ["SHORT","MEDIUM","LONG"]				# and time reference to EA outburst event
#
d = {}
snap0 = 'snap0' ; snap0_IA = 'snap0_IA'
snap1 = 'snap1' ; snap1_IA = 'snap1_IA'				# snaparr formatted as [before,during,after]  for [short,medium,long]
snap3 = 'snap3' ; snap3_IA = 'snap3_IA'             # accretion events. Uses dictionary to ensure that snaparr can be
snap4 = 'snap4' ; snap4_IA = 'snap4_IA'             # manipulated into snaparr{i,j} array. DE05 incdices also listed in ###
snap5 = 'snap5'	; snap5_IA = 'snap5_IA'
d[snap0] = [410, 460, 510, 560, 610, 660, 710, 760, 780, 800] ; d[snap0_IA] = []
	### [479, 529, 579, 629, 679, 729, 779, 829, 849, 869] ###
d[snap1] = [150, 300, 450, 600, 750, 900, 1050, 1200, 1350, 1500, 1650, 1800] ; d[snap1_IA] = []
	### [219, 369, 519, 669, 819, 969, 1119, 1269, 1419, 1569, 1719, 1869] ###
d[snap3] = [[100,150,200,250],[400,450,550,600],[1050,1150,1350,1450]]
	### [[00169,00219,00269,00319],[00469,00519,00619,00669],[01119,1219,01419,01519]] ###
d[snap3_IA] = [[290, 330, 370],[700,800,900]]
    ### [[359, 399, 439],[769,869,969]] ###
d[snap4] = [[200,209,212,220],[860,875,880,900]]	# For post-GF EA analysis add: ,[1525,1585,1595,1625]]
	### [[00269,00278,00281,00289],[00929,00944,00949,00969],[01594,01654,01664,01694]] ###
d[snap4_IA] = [370, 520, 670]
    ### [439, 589, 739] ###
d[snap5] = [[150,183,191,224],[719,755,765,801]]	# For post-GF EA analysis add: ,[1342,1377,1399,1434]]
	### [[00219,00252,00260,00293],[00788,00824,00834,00870],[01411,01446,01468,01503]] ###
d[snap5_IA] = [350, 475, 600]
    ### [419, 544, 669] ###
#
core0 = 'core0' ; core1 = 'core1'
core3 = 'core3' ; core4 = 'core4' ; core5 = 'core5'
d[core0] = [410, 610, 800]					# Choose core snapshots to be analysed for pdisc, pv and
d[core1] = [300, 1050, 1800]					# mass plots in non-episodic runs. Provide dummys for episodic
d[core3] = [] ; d[core4] = [] ; d[core5] = []			# runs to prevent I/O issues
#
logq = False         # As per referee comments - decide whether to plot log Q
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - CONSTANTS AND CONVERSIONS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
pcAU = 206265.		# Conversion from parsec -> AU
AUm = 1.496e11		# Conversion from AU -> m
G = 6.67e-11		# Gravitational constant
Msol_kg = 1.998e30	# Conversion from Solar masses -> kg
kb = 1.38e-23         # Bolztmann constant
mH = 1.67e-27         # Mass of Hydrogen atom
alpha_mm = 1.        # Artificial viscosity parameters from Morris & Monaghan (1997)
beta_mm = 2.*alpha_mm
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - MODULE IMPORTS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
import os			# Standard Python modules
import numpy as np
#
import r1			# Local modules
import pdisc
import exp_comp
import ea_comp
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
	# Before beginning loops, remove old .dat files with radial exponents for SD and T
	# Add dummy exponent-time points to ensure that exp_comp.py plots appropriate legend symbols
#
os.system('rm -r '+arch_dir+'SD_exps.dat '+arch_dir+'T_exps.dat')
for i in range(0, len(ea_run)):
	f = open(arch_dir+'SD_exps.dat','a')
	f.write(str(ea_run[i])+' 0 0 0 0 0'+'\n' )
	if (i+1 > 1):
		f.write(str(ea_run[i])+' 0 78.70 -1 0 0'+'\n' )
	f.close()
for i in range(0, len(ea_run)):
	f = open(arch_dir+'T_exps.dat','a')
	f.write(str(ea_run[i])+' 0 0 0 0'+'\n' )
	if (i+1 > 1):
		f.write(str(ea_run[i])+' 0 78.70 -1 0'+'\n' )
	f.close()
	f.close()
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# Loop over ea runs (L22) to analyse disc structure
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
for i in range(0, len(ea_run)):
#
	print "\n Run "+str(ea_run[i])+" is now being analysed \n"
#
	plotdir = arch_dir+"PLOTS/"+str(ea_run[i])+"/"
	os.system('mkdir '+str(plotdir))
#
	dat_dir_tmp = dat_dir+str(ea_run[i])+"/vK90_0i/"
#
	snaparr = d['snap'+str(ea_run[i])] ; snaparr = np.array(snaparr)
	snapcore = d['core'+str(ea_run[i])] ; snapcore = np.array(snapcore)
	snaparr_IA = d['snap'+str(ea_run[i])+'_IA'] ; snaparr_IA = np.array(snaparr_IA)
#
	# Read in rdisc.1 file
#
	hasharr_app, n_accr = r1.read(dat_dir_tmp, plotdir, ea_run[i], snaparr, plot_form)
#
	# pdisc read (Disc parameters vs. radius for defined time)
#
	pdisc.read(arch_dir, dat_dir_tmp, plotdir, col_arr, ea_run[i], \
	   snaparr, snaparr_IA, snapcore, EA_timeref, EA_lenref, r_start, r_limit, \
	   smooth, AUm, kb, mH, plot_form, logq, alpha_mm, beta_mm)
#
	# Convert all images to .eps format, suitable for publication
#
#	os.system('for f in '+plotdir+'*.pdf; do pdftops -eps $f; done')
#	os.system('rm -r '+plotdir+'*.pdf')
#
	# After looping over all selected ea runs, compare radial profile exponents for T and SD
#
plotdir = arch_dir + "/PLOTS/"
exp_comp.comp(arch_dir, dat_dir, plotdir, col_arr, plot_form)
#
	# Generate plots of ea comparative snapshots for Q, SD and T profiles
#
ea_comp.comp(arch_dir, dat_dir, plotdir, col_arr, plot_form, logq, r_start, \
   r_limit, smooth, alpha_mm, beta_mm, kb, mH)
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# Remove garbage from analysis directory, then end program
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
os.system("rm -r *.pyc")
os.system("rm -r *~")
exit()
