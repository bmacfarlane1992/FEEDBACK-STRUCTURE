"""
 exp_comp.py

 Programme to plot comparison of radial T and SD exponents for no, continuous and episodic
 feedback regimes

 Author: Benjamin MacFarlane
 Date: 26/09/2016
 Contact: bmacfarlane@uclan.ac.uk

"""
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - VARIABLE DEFINITIONS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
	# Select snapshots (from respective snaparr indices) from which to compare exponents from.
	# Also select which EA run to draw exponents from (see main.py for selection available)
	# Include relevant IA snapshots in EF runs too

snaparr = [[0, 310, 360, 410, 460, 510, 560, 610, 660, 710, 760, 780, 800], \
	[0, 150, 300, 450, 600, 750, 900, 1050, 1200, 1350, 1500, 1650, 1800], \
	[0, 100, 150, 200, 250, 290, 330, 370, 400, 450, 550, 600, 700, 800, 900, 1050, 1150, 1350, 1450], \
	[0, 200, 209, 212, 220, 370, 520, 670, 860, 875, 880, 1525, 1585, 1595, 1625], \
	[0, 150, 183, 191, 224, 350, 475, 600, 719, 755, 765, 1342, 1377, 1399, 1434] ]
ea_run = ["NRF", "CRF", "ERF-A","ERF-B","ERF-C"]
ea_run_n = [0, 1, 3, 4, 5]
#
	# Arrays of snapshots that have GI/Smooth morphology to compare average exponent values between
	# feedback regimes. Also identify the EF snapshots in outburst for each morphology.
	# Nested array format: [[NF],[CF],[EF-A],[EF-B],[EF-C]], and [[EF-A],[EF-B],[EF-C]] (for _O arrays)
#
terminal_expcomp = "TRUE"
GI = [[410, 460, 510, 560, 610, 660, 710, 760, 780, 800],[1200,1350,1500,1750], \
	[370,400,450,800,900,1050,1150],[370,520,670,860,875,880],[475,600,719,755,765]]
SMOOTH = [[310,360],[150,300,450,600,750,900,1050,1650],[100,150,200,250,290,330,550,600,700,1350,1450], \
	[200,209,212,220],[150,183,191,224,350]]
#
GI_O = [[450,1150],[875,880],[755,765]]
SMOOTH_O = [[150,200,550,1350],[209,212],[183,191]]
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - MODULE IMPORTS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import rcParams
rcParams['ps.fonttype'] = 42
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
		# # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
def comp(arch_dir, dat_dir, plotdir, col_arr, plot_form):
#
	print "\nExponent comparisons ([Sigma, T] profiles) between feedback regimes are now being plotted"
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# Read in files with exponent data
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
	SD_ea = [] ; SD_snap = [] ; SD_time = [] ; SD_exp = [] ; p_err = [] ; exp_err = []
	f = open(arch_dir+'SD_exps.dat','r')
	for line in f:
		line = line.strip() ; columns = line.split()
		SD_ea.append(float(columns[0])) ; SD_snap.append(float(columns[1]))
		SD_time.append(float(columns[2])) ; SD_exp.append(float(columns[3]))
		p_err.append(float(columns[4])) ; exp_err.append(float(columns[5]))
	f.close()
	SD_ea = np.array(SD_ea) ; SD_snap = np.array(SD_snap) ; SD_time = np.array(SD_time)
	SD_exp = np.array(SD_exp) ; p_err = np.array(p_err) ; exp_err = np.array(exp_err)
#
	T_ea = [] ; T_snap = [] ; T_time = [] ; T_exp = [] ; q_err = []
	f = open(arch_dir+'T_exps.dat','r')
	for line in f:
		line = line.strip() ; columns = line.split()
		T_ea.append(float(columns[0])) ; T_snap.append(float(columns[1]))
		T_time.append(float(columns[2])) ; T_exp.append(float(columns[3]))
		q_err.append(float(columns[4]))
	f.close()
	T_ea = np.array(T_ea) ; T_snap = np.array(T_snap)
	T_time = np.array(T_time) ; T_exp = np.array(T_exp)
	q_err = np.array(q_err)
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# Plotting of exponent comparisons for different feedback regimes
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
	# SD exponent
#
	fig = plt.figure(1)
	fig.set_figheight(15)  ; fig.set_figwidth(5)
	axcount = 0
#
	for i in range(0,len(ea_run)):
#
		axcount = axcount+1
		axpoint = 510 + axcount
		ax1 = plt.subplot(axpoint)
		if (i != len(ea_run)-1):
			ax1.xaxis.set_major_formatter(plt.NullFormatter())
		plt.xlim(79., 100.)
		plt.ylim(0, 4) ; plt.ylabel('p')
#
	# Read in accretion parameters and specify T_snap and SD_snap indices are in outburst
#
		n_accr = 0
		if (ea_run_n[i] >  1):
			t_s = [] ; t_e = []
			f = open(dat_dir+str(ea_run_n[i])+'/acc_params.dat','r')
			for line in f:
				line = line.strip() ; columns = line.split()
				t_s.append(float(columns[0])) ; t_e.append(float(columns[1]))
			f.close()
			t_s = np.array(t_s) ; t_e = np.array(t_e) ; n_accr = len(t_s)
#
	# Define which EA snapshots occur during outburst phase, then plot with relevant symbol
#
		n_out = 0 ; y_out = 0
		for j in range(0,len(snaparr[i])):
			for k in range(0,len(SD_ea)):
#
				outburst = 0
				for a in range(0, n_accr):
					if ((SD_ea[k] > 1) and (SD_time[k] > t_s[a]) and (SD_time[k] < t_e[a])):
						outburst = 1
#
				if ((snaparr[i][j] == SD_snap[k]) and (SD_ea[k] == ea_run_n[i]) and (exp_err[k] == 0) and (outburst == 0)):
					plt.scatter(SD_time[k], abs(SD_exp[k]), s=80, \
					   facecolors='none', edgecolors=col_arr[i], \
					   label = ea_run[i] if n_out == 0 else "")
#					plt.errorbar(SD_time[k], abs(SD_exp[k]), yerr = p_err[k], color='k')
					n_out = n_out + 1
				if ((snaparr[i][j] == SD_snap[k]) and (SD_ea[k] == ea_run_n[i]) and (exp_err[k] == 1) and (outburst == 0)):
					plt.scatter(SD_time[k], abs(SD_exp[k]), s=80, marker="^", \
					   facecolors='none', edgecolors=col_arr[i], \
					   label = ea_run[i] if n_out == 0 else "")
#					plt.errorbar(SD_time[k], abs(SD_exp[k]), yerr = p_err[k], color='k')
					n_out = n_out + 1
				if ((snaparr[i][j] == SD_snap[k]) and (SD_ea[k] == ea_run_n[i]) and (exp_err[k] == 0) and (outburst == 1)):
					plt.scatter(SD_time[k], abs(SD_exp[k]), s=80, \
					   facecolors=col_arr[i], edgecolors=col_arr[i], \
					   label = ea_run[i]+' Outburst' if y_out == 0 else "")
#					plt.errorbar(SD_time[k], abs(SD_exp[k]), yerr = p_err[k], color='k')
					y_out = y_out + 1
				if ((snaparr[i][j] == SD_snap[k]) and (SD_ea[k] == ea_run_n[i]) and (exp_err[k] == 1) and (outburst == 1)):
					plt.scatter(SD_time[k], abs(SD_exp[k]), s=80, marker="^", \
					   facecolors=col_arr[i], edgecolors=col_arr[i], \
					   label = ea_run[i]+' Outburst' if y_out == 0 else "")
#					plt.errorbar(SD_time[k], abs(SD_exp[k]), yerr = p_err[k], color='k')
					y_out = y_out + 1
			plt.legend(loc='upper right', fontsize = 8, scatterpoints = 1)
	plt.xlabel('Time (kyr)' )
	fig.tight_layout()
	plt.savefig(plotdir+'p_time.'+str(plot_form), format=str(plot_form), dpi=150)
	plt.clf()
#
	# T exponent comparison
#
	plt.figure(1)
	fig.set_figheight(15)  ; fig.set_figwidth(5)
#
	axcount = 0
#
	for i in range(0,len(ea_run)):
#
		axcount = axcount+1
		axpoint = 510 + axcount
		ax1 = plt.subplot(axpoint)
		if (i != len(ea_run)-1 ):
			ax1.xaxis.set_major_formatter(plt.NullFormatter())
		plt.xlim(79., 100.) ;
		plt.ylim(0, 2.5) ; plt.ylabel('q')
#
		n_accr = 0
		if (ea_run_n[i] >  1):
			t_s = [] ; t_e = []
			f = open(dat_dir+str(ea_run_n[i])+'/acc_params.dat','r')
			for line in f:
				line = line.strip() ; columns = line.split()
				t_s.append(float(columns[0])) ; t_e.append(float(columns[1]))
			f.close()
			t_s = np.array(t_s) ; t_e = np.array(t_e) ; n_accr = len(t_s)
#
		n_out = 0 ; y_out = 0
		for j in range(0,len(snaparr[i])):
			for k in range(0,len(T_ea)):
#
				outburst = 0
				for a in range(0, n_accr):
					if ((T_ea[k] > 1) and (T_time[k] > t_s[a]) and (T_time[k] < t_e[a])):
						outburst = 1
#
				if ((snaparr[i][j] == T_snap[k]) and (SD_ea[k] == ea_run_n[i]) and (exp_err[k] == 0) and (outburst == 0)):
					plt.scatter(T_time[k], abs(T_exp[k]), s=80, \
					   facecolors='none', edgecolors=col_arr[i], \
					   label = ea_run[i] if n_out == 0 else "")
#					plt.errorbar(T_time[k], abs(T_exp[k]), yerr = q_err[k], color='k')
					n_out = n_out + 1
				if ((snaparr[i][j] == T_snap[k]) and (SD_ea[k] == ea_run_n[i]) and (exp_err[k] == 1) and (outburst == 0)):
					plt.scatter(T_time[k], abs(T_exp[k]), s=80, marker="^", \
					   facecolors='none', edgecolors=col_arr[i], \
					   label = ea_run[i] if n_out == 0 else "")
#					plt.errorbar(T_time[k], abs(T_exp[k]), yerr = q_err[k], color='k')
					n_out = n_out + 1
				if ((snaparr[i][j] == T_snap[k]) and (T_ea[k] == ea_run_n[i]) and (exp_err[k] == 0) and (outburst == 1)):
					plt.scatter(T_time[k], abs(T_exp[k]), s=80, \
					   facecolors=col_arr[i], edgecolors=col_arr[i], \
					   label = ea_run[i]+' Outburst' if y_out == 0 else "")
#					plt.errorbar(T_time[k], abs(T_exp[k]), yerr = q_err[k], color='k')
					y_out = y_out + 1
				if ((snaparr[i][j] == T_snap[k]) and (T_ea[k] == ea_run_n[i]) and (exp_err[k] == 1) and (outburst == 1)):
					plt.scatter(T_time[k], abs(T_exp[k]), s=80, marker="^", \
					   facecolors=col_arr[i], edgecolors=col_arr[i], \
					   label = ea_run[i]+' Outburst' if y_out == 0 else "")
#					plt.errorbar(T_time[k], abs(T_exp[k]), yerr = q_err[k], color='k')
					y_out = y_out + 1
			plt.legend(loc='upper right', fontsize = 8, scatterpoints = 1)
	plt.xlabel('Time (kyr)' )
	fig.tight_layout()
	plt.savefig(plotdir+'q_time.'+str(plot_form), format=str(plot_form), dpi=150)
	plt.clf()
#
#
#
	# p + q exponent comparison
#
	plt.figure(1)
	fig.set_figheight(15)  ; fig.set_figwidth(5)
#
	axcount = 0
#
	for i in range(0,len(ea_run)):
#
		axcount = axcount+1
		axpoint = 510 + axcount
		ax1 = plt.subplot(axpoint)
		if (i != len(ea_run)-1 ):
			ax1.xaxis.set_major_formatter(plt.NullFormatter())
		plt.xlim(79., 100.) ;
		plt.ylim(0, 7) ; plt.ylabel('p + q')
#
		n_accr = 0
		if (ea_run_n[i] >  1):
			t_s = [] ; t_e = []
			f = open(dat_dir+str(ea_run_n[i])+'/acc_params.dat','r')
			for line in f:
				line = line.strip() ; columns = line.split()
				t_s.append(float(columns[0])) ; t_e.append(float(columns[1]))
			f.close()
			t_s = np.array(t_s) ; t_e = np.array(t_e) ; n_accr = len(t_s)
#
		n_out = 0 ; y_out = 0
		for j in range(0,len(snaparr[i])):
			for k in range(0,len(T_ea)):
#
				outburst = 0
				for a in range(0, n_accr):
					if ((T_ea[k] > 1) and (T_time[k] > t_s[a]) and (T_time[k] < t_e[a])):
						outburst = 1
#
				if ((snaparr[i][j] == T_snap[k]) and (SD_ea[k] == ea_run_n[i]) \
				   and (exp_err[k] == 0) and (outburst == 0)):
					plt.scatter(T_time[k], abs(T_exp[k] + SD_exp[k]), s=80, \
					   facecolors='none', edgecolors=col_arr[i], \
					   label = ea_run[i] if n_out == 0 else "")
#					plt.errorbar(T_time[k], abs(T_exp[k]), yerr = q_err[k], color='k')
					n_out = n_out + 1
				if ((snaparr[i][j] == T_snap[k]) and (SD_ea[k] == ea_run_n[i]) \
				   and (exp_err[k] == 1) and (outburst == 0)):
					plt.scatter(T_time[k], abs(T_exp[k] + SD_exp[k]), s=80, marker="^", \
					   facecolors='none', edgecolors=col_arr[i], \
					   label = ea_run[i] if n_out == 0 else "")
#					plt.errorbar(T_time[k], abs(T_exp[k]), yerr = q_err[k], color='k')
					n_out = n_out + 1
				if ((snaparr[i][j] == T_snap[k]) and (T_ea[k] == ea_run_n[i]) \
				   and (exp_err[k] == 0) and (outburst == 1)):
					plt.scatter(T_time[k], abs(T_exp[k] + SD_exp[k]), s=80, \
					   facecolors=col_arr[i], edgecolors=col_arr[i], \
					   label = ea_run[i]+' Outburst' if y_out == 0 else "")
#					plt.errorbar(T_time[k], abs(T_exp[k]), yerr = q_err[k], color='k')
					y_out = y_out + 1
				if ((snaparr[i][j] == T_snap[k]) and (T_ea[k] == ea_run_n[i]) \
				   and (exp_err[k] == 1) and (outburst == 1)):
					plt.scatter(T_time[k], abs(T_exp[k] + SD_exp[k]), s=80, marker="^", \
					   facecolors=col_arr[i], edgecolors=col_arr[i], \
					   label = ea_run[i]+' Outburst' if y_out == 0 else "")
#					plt.errorbar(T_time[k], abs(T_exp[k]), yerr = q_err[k], color='k')
					y_out = y_out + 1
			plt.legend(loc='upper right', fontsize = 8, scatterpoints = 1)
	plt.xlabel('Time (kyr)' )
	fig.tight_layout()
	plt.savefig(plotdir+'exp_analysis1.'+str(plot_form), format=str(plot_form), dpi=150)
	plt.clf()
#
#
	# p + q exponent comparison
#
	plt.figure(1)
	fig.set_figheight(15)  ; fig.set_figwidth(5)
#
	axcount = 0
#
	for i in range(0,len(ea_run)):
#
		axcount = axcount+1
		axpoint = 510 + axcount
		ax1 = plt.subplot(axpoint)
		if (i != len(ea_run)-1 ):
			ax1.xaxis.set_major_formatter(plt.NullFormatter())
		plt.xlim(79., 100.) ;
		plt.ylim(0, 7) ; plt.ylabel('p - (q/2) - 3/2')
#
		n_accr = 0
		if (ea_run_n[i] >  1):
			t_s = [] ; t_e = []
			f = open(dat_dir+str(ea_run_n[i])+'/acc_params.dat','r')
			for line in f:
				line = line.strip() ; columns = line.split()
				t_s.append(float(columns[0])) ; t_e.append(float(columns[1]))
			f.close()
			t_s = np.array(t_s) ; t_e = np.array(t_e) ; n_accr = len(t_s)
#
		n_out = 0 ; y_out = 0
		for j in range(0,len(snaparr[i])):
			for k in range(0,len(T_ea)):
#
				outburst = 0
				for a in range(0, n_accr):
					if ((T_ea[k] > 1) and (T_time[k] > t_s[a]) and (T_time[k] < t_e[a])):
						outburst = 1
#
				if ((snaparr[i][j] == T_snap[k]) and (SD_ea[k] == ea_run_n[i]) \
				   and (exp_err[k] == 0) and (outburst == 0)):
					plt.scatter(T_time[k], abs(SD_exp[k] + (T_exp[k]/2.) - (3./2.) ), \
					   s=80, facecolors='none', edgecolors=col_arr[i], \
					   label = ea_run[i] if n_out == 0 else "")
#					plt.errorbar(T_time[k], abs(T_exp[k]), yerr = q_err[k], color='k')
					n_out = n_out + 1
				if ((snaparr[i][j] == T_snap[k]) and (SD_ea[k] == ea_run_n[i]) \
				   and (exp_err[k] == 1) and (outburst == 0)):
					plt.scatter(T_time[k], abs(SD_exp[k] + (T_exp[k]/2.) - (3./2.)), \
					   s=80, marker="^", facecolors='none', edgecolors=col_arr[i], \
					   label = ea_run[i] if n_out == 0 else "")
#					plt.errorbar(T_time[k], abs(T_exp[k]), yerr = q_err[k], color='k')
					n_out = n_out + 1
				if ((snaparr[i][j] == T_snap[k]) and (T_ea[k] == ea_run_n[i]) \
				   and (exp_err[k] == 0) and (outburst == 1)):
					plt.scatter(T_time[k], abs(SD_exp[k] + (T_exp[k]/2.) - (3./2.)), \
					   s=80, facecolors=col_arr[i], edgecolors=col_arr[i], \
					   label = ea_run[i]+' Outburst' if y_out == 0 else "")
#					plt.errorbar(T_time[k], abs(T_exp[k]), yerr = q_err[k], color='k')
					y_out = y_out + 1
				if ((snaparr[i][j] == T_snap[k]) and (T_ea[k] == ea_run_n[i]) \
				   and (exp_err[k] == 1) and (outburst == 1)):
					plt.scatter(T_time[k], abs(SD_exp[k] + (T_exp[k]/2.) - (3./2.)), \
					   s=80, marker="^", facecolors=col_arr[i], edgecolors=col_arr[i], \
					   label = ea_run[i]+' Outburst' if y_out == 0 else "")
#					plt.errorbar(T_time[k], abs(T_exp[k]), yerr = q_err[k], color='k')
					y_out = y_out + 1
			plt.legend(loc='upper right', fontsize = 8, scatterpoints = 1)
	plt.xlabel('Time (kyr)' )
	fig.tight_layout()
	plt.savefig(plotdir+'exp_analysis2.'+str(plot_form), format=str(plot_form), dpi=150)
	plt.clf()
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	# Comparison of p and q exponents as a function of both feedback regime and disc morphology
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
	if (terminal_expcomp == "TRUE"):
#
	# GI morphology
#
# p exponent averages
#
		print "\nGI Morphology \n"
		print "p exponent \n"
		p_EFtmp = []
		for i in range(0, len(GI)):
			p_tmp = []
			if (i <= 1):
				GI_ea_n = i
			elif (i > 1):
				GI_ea_n = i+1
			for j in range(0, len(GI[i])):
				for k in range(0, len(SD_exp)):
					if ((GI[i][j] == SD_snap[k]) and (GI_ea_n == SD_ea[k])):
						p_tmp.append(float(SD_exp[k]))
			if ((i <=1) and (len(p_tmp) == 0)):
				print ea_run[i],": N/A"
			elif ((i <= 1) and (len(p_tmp) > 0)):
				print ea_run[i],": ", np.round(abs(np.mean(p_tmp)),1)
			elif (i > 1):
				p_EFtmp = p_EFtmp + p_tmp
		print "ERF runs: ", np.round(abs(np.mean(p_EFtmp)),1), "\n"
#
# q exponent averages
#
		print "q exponent \n"
#
		q_EF = [] ; q_EFO = []
		for i in range(0, len(GI)):
			q_EFtmp = [] ; q_EFOtmp = []
#
			if (i <= 1):
				q_tmp = []
				GI_ea_n = i
				for j in range(0, len(GI[i])):
					for k in range(0, len(T_exp)):
						if ((GI[i][j] == T_snap[k]) and (GI_ea_n == T_ea[k])):
							q_tmp.append(float(T_exp[k]))
				if (len(q_tmp) == 0):
					print ea_run[i],": N/A"
				elif (len(q_tmp) > 0):
					print ea_run[i],": ", np.round(abs(np.mean(q_tmp)),1)
#
			elif (i > 1):
				GI_ea_n = i+1
				for j in GI[i]:
					if j in GI_O[i-2]:
						for k in range(0, len(T_exp)):
							if ((j == T_snap[k]) and (GI_ea_n == T_ea[k])):
								q_EFOtmp.append(T_exp[k])
					else:
						for k in range(0, len(T_exp)):
							if ((j == T_snap[k]) and (GI_ea_n == T_ea[k])):
								q_EFtmp.append(T_exp[k])
				q_EF = q_EF + q_EFtmp
				q_EFO = q_EFO + q_EFOtmp
#
		print "ERF (quiescent) snapshots: ", np.round(abs(np.mean(q_EF)),1)
		print "ERF-O snapshots: ", np.round(abs(np.mean(q_EFO)),1), "\n"
#
#
	# SMOOTH morphology
#
# p exponent averages
#
		print "SMOOTH Morphology \n"
		print "p exponent \n"
		p_EFtmp = []
		for i in range(0, len(SMOOTH)):
			p_tmp = []
			if (i <= 1):
				SMOOTH_ea_n = i
			elif (i > 1):
				SMOOTH_ea_n = i+1
			for j in range(0, len(SMOOTH[i])):
				for k in range(0, len(SD_exp)):
					if ((SMOOTH[i][j] == SD_snap[k]) and (SMOOTH_ea_n == SD_ea[k])):
							p_tmp.append(float(SD_exp[k]))
			if ((i <=1) and (len(p_tmp) == 0)):
				print ea_run[i],": N/A"
			elif ((i <= 1) and (len(p_tmp) > 0)):
				print ea_run[i],": ", np.round(abs(np.mean(p_tmp)),1)
			elif (i > 1):
				p_EFtmp = p_EFtmp + p_tmp
		print "ERF runs: ", np.round(abs(np.mean(p_EFtmp)),1), "\n"
#
# q exponent averages
#
		print "q exponent \n"
#
		q_EF = [] ; q_EFO = []
		for i in range(0, len(SMOOTH)):
			q_EFtmp = [] ; q_EFOtmp = []
#
			if (i <= 1):
				q_tmp = []
				SMOOTH_ea_n = i
				for j in range(0, len(SMOOTH[i])):
					for k in range(0, len(T_exp)):
						if ((SMOOTH[i][j] == T_snap[k]) and (SMOOTH_ea_n == T_ea[k])):
							q_tmp.append(float(T_exp[k]))
				if (len(q_tmp) == 0):
					print ea_run[i],": N/A"
				elif (len(q_tmp) > 0):
					print ea_run[i],": ", np.round(abs(np.mean(q_tmp)),1)
#
			elif (i > 1):
				SMOOTH_ea_n = i+1
				for j in SMOOTH[i]:
					if j in SMOOTH_O[i-2]:
						for k in range(0, len(T_exp)):
							if ((j == T_snap[k]) and (SMOOTH_ea_n == T_ea[k])):
								q_EFOtmp.append(T_exp[k])
					else:
						for k in range(0, len(T_exp)):
							if ((j == T_snap[k]) and (SMOOTH_ea_n == T_ea[k])):
								q_EFtmp.append(T_exp[k])
				q_EF = q_EF + q_EFtmp
				q_EFO = q_EFO + q_EFOtmp
#
		print "ERF (quiescent) snapshots: ", np.round(abs(np.mean(q_EF)),1)
		print "ERF-O snapshots: ", np.round(abs(np.mean(q_EFO)),1), "\n"
