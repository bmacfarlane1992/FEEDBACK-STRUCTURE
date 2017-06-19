#
# r1.py
#
# Python program to read in full rdisc.1 data structures from vK90_0i runs (/rdisc_DS), to evaluate
# evolution of simulation disc mass and radius based on set criteria.
# Script also evaluates disc mass and radius based on imposed criteria, for chosen snapshots
#
#
# Author: Benjamin MacFarlane
# Date: 26/09/2016
# Contact: bmacfarlane@uclan.ac.uk
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - VARIABLE DEFINITIONS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
time_check = "FALSE"		# Choose whether ("TRUE") or not ("FALSE") to output time analysis of run
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - MODULE IMPORTS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['ps.fonttype'] = 42
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
def read(dat_dir, plotdir, ea_run, snaparr, plot_form):
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# Read rdisc.1 file if original DS runs are being analysed
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
	print "rdisc.1 file being read // accretion events being noted // plots being generated"
#
	arch_ext = 'rdisc_DS/DE05.rdisc.1'
	filename = dat_dir+arch_ext
#
	# Define arrays to fill
#
	time = [] ; m_s = [] ; m_mri_d = [] ; r_d_kep = [] ; m_d_kep = [] ; r_d_piv = []
	m_d_piv = [] ; r_d_sigALMA = [] ; m_d_sigALMA = [] ; mcloud1 = [] ; acctag_num = []
#
	# Read in ea runs with episodic accretion
#
	if (float(ea_run) > 1):
		f = open(filename, 'r')
		for line in f:
			line = line.strip()
			columns = line.split()
			time.append(float(columns[0])/1000.) ; m_s.append(float(columns[1]))
			m_mri_d.append(float(columns[2])) ; r_d_kep.append(float(columns[3]))
			m_d_kep.append(float(columns[4])) ; r_d_piv.append(float(columns[5]))
			m_d_piv.append(float(columns[6])) ; r_d_sigALMA.append(float(columns[11]))
			m_d_sigALMA.append(float(columns[12])) ; acctag_num.append(float(columns[14]))
			mcloud1.append(float(columns[17]))
		f.close()
#
	# Read in ea runs without episodic accretion
#
	if (float(ea_run) <= 1):
		f = open(filename, 'r')
		for line in f:
			line = line.strip()
			columns = line.split()
			time.append(float(columns[0])/1000.) ; m_s.append(float(columns[1]))
			m_mri_d.append(float(columns[2])) ; r_d_kep.append(float(columns[3]))
			m_d_kep.append(float(columns[4])) ; r_d_piv.append(float(columns[5]))
			m_d_piv.append(float(columns[6])) ; r_d_sigALMA.append(float(columns[11]))
			m_d_sigALMA.append(float(columns[12])) ; mcloud1.append(float(columns[15]))
			acctag_num.append(0)
		f.close()
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# Carry out analysis of file in preparation for plotting
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
	# Define time indices for accretion times, in order to generate continuous fill
#
	n_snaps = int(len(time))
#
	hasharr = []
	jstart = 0
	for i in range(0, n_snaps-1):
		for j in range(jstart, n_snaps-1):
			if (int(acctag_num[j]) == 1):
				hasharr.append(j)
				for k in range(j, n_snaps):
					if(int(acctag_num[k]) != int(acctag_num[j])):
						hasharr.append(k) ; jstart = k ; break
				break
	hasharr_app=[]
	for i in hasharr:
		if i not in hasharr_app:
			hasharr_app.append(i)
	if (acctag_num[n_snaps-1] == "1"):
		hasharr_app.append(n_snaps-1)
	n_accr = len(hasharr_app)/2
	hasharr_app = np.reshape(hasharr_app, (n_accr, 2))
#
	hasharr_tmp = hasharr_app*0.01 + min(time)
#
	# Print time limits for individual runs for refining of times plotted
#
	if (time_check == "TRUE"):
		print "Minimum time of run is: ", min(time), " kyr"
		print "Maximum time of run is: ", max(time), " kyr"
#
	# Check the times of snapshots entered for before, during and after accretion event chosen
#
		print "\n"
		for i in range(0,n_accr):
			print "Accretion event "+str(i+1)+" begins at snapshot "+ \
			   str(hasharr_app[i][0])+" for "+ \
			   str((hasharr_app[i][1]-hasharr_app[i][0])*10)+" years" \
			   " until snapshot "+str(hasharr_app[i][1])
		print "\n The final snapshot of ea run "+str(ea_run)+" is "+str(len(time))+"\n"
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# Plot rdisc.1 temporally evolved parameters
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
	# Plot mass of star and MRI disc vs. time
#
	plt.figure(1)
	ax1 = plt.subplot(111)
	line1 = plt.plot(time, m_s, color = 'b', label="Stellar Mass")
	line2 = plt.plot(time, m_mri_d, color = 'r', label="MRI disc mass")
	plt.ylabel("Mass "+(r'(M$_{\odot}$)'))
	plt.xlabel('Time (kyr)') ; plt.xlim([int(min(time)),int(max(time))])
	legend = plt.legend(loc = 'upper left', fontsize=11)
	plt.savefig(plotdir+str(ea_run)+'_iad_star_mass.'+str(plot_form), \
	   format=str(plot_form), dpi=150)
	plt.clf()
#
	# Plot mass of disc with Keplerian and azimuthal velocity criteria
#
	plt.figure(1)
	ax1 = plt.subplot(111)
	line11 = plt.plot(time, m_d_kep, label="Keplerian")
	line22 = plt.plot(time, m_d_piv, label="Azimuthal")
	line33 = plt.plot(time, m_d_sigALMA, label="ALMA SD")
	for i in range(0,n_accr):
		plt.fill_between([hasharr_tmp[i][0],hasharr_tmp[i][1]], \
		0, ax1.get_ylim()[1], color='k', alpha = 0.25)
	plt.ylabel("Criterion Mass "+(r'(M$_{\odot}$)'))
	plt.xlabel('Time (kyr)') ; plt.xlim([int(min(time)),math.ceil(max(time))])
	legend = plt.legend(loc = 'upper left', fontsize=11)
	plt.savefig(plotdir+str(ea_run)+'_disc_mass.'+str(plot_form), \
	   format=str(plot_form), dpi=150)
	plt.clf()
#
	# Plot radius of disc with Keplerian and azimuthal velcity criteria
#
	plt.figure(1)
	ax1 = plt.subplot(111)
	line11 = plt.plot(time, r_d_kep, label="Keplerian")
	line22 = plt.plot(time, r_d_piv, label="Azimuthal")
	line33 = plt.plot(time, r_d_sigALMA, label="ALMA SD")
	for i in range(0,n_accr):
		plt.fill_between([hasharr_tmp[i][0],hasharr_tmp[i][1]], \
		   0, ax1.get_ylim()[1], color='k', alpha = 0.25)
	legend = plt.legend(loc = 'upper left', fontsize=11)
	plt.xlabel('Time (kyr)') ; plt.xlim([int(min(time)),math.ceil(max(time))])
	plt.ylim(0, ax1.get_ylim()[1]) ; plt.ylabel("Criterion Radius (AU)")
	plt.savefig(plotdir+str(ea_run)+'_disc_radius.'+str(plot_form), \
	   format=str(plot_form), dpi=150)
	plt.clf()
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# Write hasharr_app data to ensure mass_comp.py successfully traces accretion events
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
	if (n_accr != 0):
		f = open(dat_dir+'../acc_params.dat','w')
		for i in range(0, n_accr):
			f.write( str(hasharr_tmp[i][0])+' '+str(hasharr_tmp[i][1])+'\n' )
		f.close()
#
#
	return hasharr_app, n_accr
