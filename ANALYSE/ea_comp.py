#
# ea_comp.py
#
# Analysis script to compare Q, T and SD vs. r profiles between different ea runs
#
# Author: Benjamin MacFarlane
# Date: 22/09/2016
# Contact: bmacfarlane@uclan.ac.uk
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - VARIABLE DEFINITIONS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
snaps = [800, 900, 1050, 1150]	# Snapshot and EA run association for radial profiles being compared
ea_snaps = [0, 1, 3, 3]
#snaps = [1350, 1500, 1650, 1800]
#ea_snaps = [1, 1, 1, 1]
ea_leg = ["NRF", "CRF", "ERF", "ERF-O"]
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - MODULE IMPORTS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
import os			# Standard Python modules
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['ps.fonttype'] = 42
from scipy import interpolate
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
def comp(arch_dir, dat_dir, plotdir, col_arr, plot_form, logq, r_start, r_limit, \
   smooth, alpha_mm, beta_mm, kb, mH):
#
	print "EA-dependent Q, Sigma and T profiles are being plotted"
#
	file_n = len(snaps)
#
	# Arrays to fill
#
	time = [[] for i in range(file_n)] ; r = [[] for i in range(file_n)]
	Q = [[] for i in range(file_n)] ; T = [[] for i in range(file_n)]
	sig = [[] for i in range(file_n)] ; H = [[] for i in range(file_n)]
	h_ave = [[] for i in range(file_n)] ; alpha_ss_lp = [[] for i in range(file_n)]
#
	# Define respective DE05 pointers
#
	file_list = []
	for a in range(0, len(snaps)):
		if (snaps[a] < (1000-70)):
			file_list.append('pdisc/DE05.du.00'+str(snaps[a]+69)+'.pdisc.1')
		elif (snaps[a] > (1000-70)):
			file_list.append('pdisc/DE05.du.0'+str(snaps[a]+69)+'.pdisc.1')
#
	for i in range(0, file_n):
#
	# Extract variables from file loop and convert to numpy arrays
#
		file_list[i] = dat_dir+str(ea_snaps[i])+"/"+"vK90_90i/"+file_list[i]
		f = open(file_list[i], 'r')
#
		header = f.readline()
#
		j = 0
		for line in f:
			line = line.strip() ; columns = line.split()
			time[i].append(float(columns[0])/1000.) ; r[i].append(float(columns[1]))
			Q[i].append(float(columns[2])) ; T[i].append(float(columns[3]))
			sig[i].append(float(columns[4]))
			H[i].append( np.sqrt( (kb * T[i][j]) / (2.3 * mH) ) / (float(columns[13])*1000.) )
			h_ave[i].append(float(columns[14]))
			alpha_ss_lp[i].append( 1./10. * alpha_mm * (h_ave[i][j] / ( H[i][j] * r[i][j] ) ) )
			j = j + 1
		f.close()
#
	time = np.array(time) ; r = np.array(r) ; Q = np.array(Q)
	T = np.array(T) ; sig = np.array(sig) ; H = np.array(H)
	h_ave = np.array(h_ave) ; alpha_ss_lp = np.array(alpha_ss_lp)
#
# Begin plotting comparison of Q, T, SD and alpha_ss
#
	fig = plt.figure(1)
	fig.set_size_inches(6.0, 6.0)
#
	r_int = np.arange(r_start,r_limit,r_limit/(r_limit/smooth))
#
	ax1 = plt.subplot(411)
	for i in range(0, len(snaps)):
#
		tck = interpolate.splrep(r[i,r_start:r_limit], \
		   T[i,r_start:r_limit])
		y_int = interpolate.splev(r_int, tck, der = 0)
#
		line1 = plt.plot(r_int, y_int, color = col_arr[i])
		ax1.set_xscale('log') ; ax1.set_xlim(r_start, 150) ; plt.xticks(fontsize = 9)
		plt.ylabel("log T (K)", fontsize = 9, labelpad=0.5)
		ax1.set_yscale('log') ; ax1.get_yaxis().set_label_coords(-0.1,0.5)
		ax1.set_ylim(0, 1500) ; plt.yticks(fontsize = 9)
		ax1.xaxis.set_major_formatter(plt.NullFormatter())
#
	ax2 = plt.subplot(412)
	for i in range(0, len(snaps)):
#
		tck = interpolate.splrep(r[i,r_start:r_limit], \
		   sig[i,r_start:r_limit])
		y_int = interpolate.splev(r_int, tck, der = 0)
#
		GI1_x = 35. ; GI2_x = 55. ; GI3_x = 90.
		GI_y0 = 1000. ; GI1_y1 = 300. ; GI2_y1 = 150. ; GI3_y1 = 50.
#
		line2 = plt.plot(r_int, y_int, color = col_arr[i])
		ax2.arrow(GI1_x, GI_y0, 0., GI1_y1-GI_y0, head_width=3., head_length=50., fc='k', ec='k')
		ax2.arrow(GI2_x, GI_y0, 0., GI2_y1-GI_y0, head_width=3.*(GI2_x/GI1_x), head_length=50./(GI1_y1/GI2_y1), fc='k', ec='k')
		ax2.arrow(GI3_x, GI_y0, 0., GI3_y1-GI_y0, head_width=3.*(GI3_x/GI1_x), head_length=50./(GI1_y1/GI3_y1), fc='k', ec='k')
		ax2.set_xscale('log') ; ax2.set_xlim(r_start, 150) ; plt.xticks(fontsize = 9)
		plt.ylabel("log "+r'$\Sigma$ (g cm$^{-2}$)', fontsize = 9, labelpad=0.5)
		ax2.set_yscale('log') ; ax2.get_yaxis().set_label_coords(-0.1,0.5) ;
		ax2.set_ylim(0, 1200) ; plt.yticks(fontsize = 9)
		ax2.xaxis.set_major_formatter(plt.NullFormatter())
#
	ax3 = plt.subplot(413)
	for i in range(0, len(snaps)):
#
		tck = interpolate.splrep(r[i,r_start:r_limit], \
		   Q[i,r_start:r_limit])
		y_int = interpolate.splev(r_int, tck, der = 0)
#
		line3 = plt.plot(r_int, y_int, color = col_arr[i], \
		   label = ea_leg[i])
		ax3.axhline(y=1.,color='k',ls='dashed')
		ax3.set_xscale('log')
		ax3.set_xlim(r_start, 150) ; plt.xticks(fontsize = 9)
		plt.ylabel("Q", fontsize = 9)
		ax3.get_yaxis().set_label_coords(-0.1,0.5)
		ax3.xaxis.set_major_formatter(plt.NullFormatter())
#
		if logq is True:
			ax3.set_yscale('log')
			ax3.set_ylim(0.25,4)
		elif logq is False:
			ytck = [0.0, 1.0, 2.0, 3.0, 4.0]
			ax3.set_ylim(0, 4)
			ax3.set_yticks(ytck) ; plt.yticks(fontsize = 9)
#
#
	ax4 = plt.subplot(414)
#
	for i in range(0, len(snaps)):
#
		tck = interpolate.splrep(r[i,r_start:r_limit], \
		   alpha_ss_lp[i,r_start:r_limit])
		y_int = interpolate.splev(r_int, tck, der = 0)
#
		line4 = plt.plot(r_int, y_int, color = col_arr[i], \
		   label = ea_leg[i])
		ax4.set_xlim(r_start, 150) ; plt.xticks(fontsize = 9)
		ax4.set_yscale('log')
#
		plt.xlabel("Radius (AU)", fontsize = 9) ; ax4.set_xscale('log')
		plt.ylabel(r'$\alpha_{SS}$', fontsize = 9, labelpad=1.5)
		ax4.get_yaxis().set_label_coords(-0.1,0.5)
		ax4.set_ylim(0.01, 1.)
		plt.yticks(fontsize = 9)
	plt.legend(loc='upper left', fontsize = 7.5)
#
	plt.gcf().subplots_adjust(bottom=0.075)
	plt.savefig(str(plotdir)+'AziProf_comp.'+str(plot_form), format=str(plot_form), dpi=150)
	plt.clf()
#
#
