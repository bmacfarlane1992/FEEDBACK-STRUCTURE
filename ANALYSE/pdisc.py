#
# pdisc.py
#
# Python program to read in pdisc (radially varied, temporally static) data structures
# and compute surface density (SD) and temperature (T)  power exponents
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
def func(x, a, b):                          # Linear function used to fit power indices to log-log Sig/T vs. radius distributions
  return a*x + b
#
loglin = "FALSE"                            # Choose either log-linear ("TRUE") or log-log ("FALSE") distributions of SD and T
#
d = {}                                          # Radial value (AU) to start of p and q fits for each snaparr (see main.py for formatting)
r_fitS_0 = 'r_fitS_0' ; r_fitS_1 = 'r_fitS_1'    # using dictionaries
r_fitS_3 = 'r_fitS_3' ; r_fitS_3_IA = 'r_fitS_3_IA'         # Set dictionary and relevent start/stop radii for snapshots inbetween accretion (IA)
r_fitS_4 = 'r_fitS_4' ; r_fitS_4_IA = 'r_fitS_4_IA'         # bursts too.
r_fitS_5 = 'r_fitS_5' ; r_fitS_5_IA = 'r_fitS_5_IA'
d[r_fitS_0] = [20, 15, 15, 20, 20, 20, 20, 20, 20, 20]
d[r_fitS_1] = [5, 10, 20, 20, 20, 20, 20, 20, 20, 20, 15, 15]
d[r_fitS_3] = [[5,5,10,10],[20,20,20,15],[20,20,20,20]]
d[r_fitS_3_IA] = [[15, 15, 15],[20,20,20]]
d[r_fitS_4] = [[15,15,22,12],[20,12,30,20]]
d[r_fitS_4_IA] = [20, 20, 20]
d[r_fitS_5] = [[10,12,12,10],[15,15,15,25]]
d[r_fitS_5_IA] = [15, 25, 25]
#
r_fitE_0 = 'r_fitE_0' ; r_fitE_1 = 'r_fitE_1'    # As above, for end of fits for each snaparr.
r_fitE_3 = 'r_fitE_3' ; r_fitE_3_IA = 'r_fitE_3_IA'
r_fitE_4 = 'r_fitE_4' ; r_fitE_4_IA = 'r_fitE_4_IA'
r_fitE_5 = 'r_fitE_5' ; r_fitE_5_IA = 'r_fitE_5_IA'
d[r_fitE_0] = [60, 60, 60, 70, 70, 70, 70, 70, 70, 70]
d[r_fitE_1] = [20, 50, 60, 80, 80, 80, 80, 80, 70, 70, 65, 80]
d[r_fitE_3] = [[15,20,40,40],[60,60,60,50],[60,60,60,60]]
d[r_fitE_3_IA] = [[50, 50, 50],[70,70,70]]
d[r_fitE_4] = [[40,40,40,40],[70,45,75,30]]
d[r_fitE_4_IA] = [55, 60, 70]
d[r_fitE_5] = [[18,40,50,30],[60,60,60,40]]
d[r_fitE_5_IA] = [45, 60, 70]
#
exp_err_0 = 'exp_err_0' ; exp_err_1 = 'exp_err_1'    # Error tag on fitting to SD profile, to indicate "bumps" in fitting routine
exp_err_3 = 'exp_err_3' ; exp_err_3_IA = 'exp_err_3_IA'
exp_err_4 = 'exp_err_4' ; exp_err_4_IA = 'exp_err_4_IA'
exp_err_5 = 'exp_err_5' ; exp_err_5_IA = 'exp_err_5_IA'
d[exp_err_0] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
d[exp_err_1] = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
d[exp_err_3] = [[1, 1, 0, 0],[0, 0, 0, 1],[1, 0, 0, 0]]
d[exp_err_3_IA] = [[0, 0, 0],[0, 0, 0]]
d[exp_err_4] = [[1, 1, 1, 1],[0, 0, 0, 1]]
d[exp_err_4_IA] = [0, 0, 0]
d[exp_err_5] = [[1, 1, 0, 1],[0, 0, 0, 1]]
d[exp_err_5_IA] = [0, 0, 0]
#
logq = True         # As per referee comments - decide whether to plot log Q
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
from scipy.optimize import curve_fit
from scipy import interpolate
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
def read(arch_dir, dat_dir, plotdir, col_arr, ea_run, snaparr, snaparr_IA, snapcore, \
    EA_timeref, EA_lenref, r_start, r_limit, smooth, AUm, kb, mH, plot_form, \
    logq, alpha_mm, beta_mm):
#
    print("pdisc files being read // exponents being computed // plots being generated")
#
    # Define number of files to be read dependent on snaparr dimensions
#
    if (snaparr.ndim == 1):
        file_n = len(snaparr)
    elif (snaparr.ndim == 2):
        file_n = len(snaparr)*len(snaparr[0])
#
    # Arrays to fill
#
    time = [[] for i in range(file_n)] ; r = [[] for i in range(file_n)]
    Q = [[] for i in range(file_n)] ; T = [[] for i in range(file_n)]
    sig = [[] for i in range(file_n)] ; omega = [[] for i in range(file_n)]
    H = [[] for i in range(file_n)]
    vkep = [[] for i in range(file_n)] ; h_ave = [[] for i in range(file_n)]
#
    # Create string array as pointer to pdisc files based on snaparr
#
    file_list = []
    if (snaparr.ndim == 1):
        snaparr_tmp = np.array([0]*len(snaparr))
        snapcore_tmp = np.array([0]*len(snapcore))
        fcount = 0
        for a in range(0, len(snaparr)):
            if (snaparr[a] < (1000-70)):
                file_list.append(dat_dir+'pdisc/DE05.du.00'+ \
                   str(snaparr[a]+69)+'.pdisc.1')
            elif (snaparr[a] > (1000-70)):
                file_list.append(dat_dir+'pdisc/DE05.du.0'+ \
                   str(snaparr[a]+69)+'.pdisc.1')
            for b in range(0, len(snapcore)):
                if (snaparr[a] == snapcore[b]):
                    snapcore_tmp[b] = fcount
            snaparr_tmp[a] = fcount
            fcount = fcount + 1
#
    elif (snaparr.ndim == 2):
        snaparr_tmp = np.array([[0]*len(snaparr[0])]*len(snaparr))
        fcount = 0
        for a in range(0, len(snaparr)):
            for b in range(0, len(snaparr[0])):
                if (snaparr[a][b] < (1000-70)):
                    file_list.append(dat_dir+'pdisc/DE05.du.00'+ \
                       str(snaparr[a][b]+69)+'.pdisc.1')
                elif (snaparr[a][b] > (1000-70)):
                    file_list.append(dat_dir+'pdisc/DE05.du.0'+ \
                       str(snaparr[a][b]+69)+'.pdisc.1')
                snaparr_tmp[a][b] = fcount
                fcount = fcount + 1
#
    # Loop over files in pdisc directory to read all temporally evolved values
#
    alpha_ss_lp = [[] for i in range(file_n)]
    alpha_ss_meru1 = [[] for i in range(file_n)]
    alpha_ss_meru2 = [[] for i in range(file_n)]
#
    for i in range(0, file_n):
#
# "#t(yr)/r(AU)/Q/T/sig(gcm-2)/Omeg/DiscM(Msun)/rInM/M*/Mdis/vr/vthe/vz/vkep/h/h_mid/rpart"
#
        f = open(file_list[i], 'r')
        header = f.readline()
        j = 0
        for line in f:
            line = line.strip() ; columns = line.split()
            time[i].append(float(columns[0])/1000.) ; r[i].append(float(columns[1]))
            Q[i].append(float(columns[2])) ; T[i].append(float(columns[3]))
            sig[i].append(float(columns[4])) ; omega[i].append(float(columns[5]))
            H[i].append( np.sqrt( (kb * T[i][j]) / (2.3 * mH) ) / (float(columns[13])*1000.) )
            vkep[i].append(float(columns[13])) ; h_ave[i].append(float(columns[14]))
#
    # In file read, compute effective Shakura-Sunyaev alpha from Lodato & Rice (2010)
    # and Meru et al. (2012) formalisms
#
            alpha_ss_lp[i].append( 1./10. * alpha_mm * (h_ave[i][j] / ( H[i][j] * r[i][j] )) )
            alpha_ss_meru1[i].append( 31./225. * alpha_mm * (h_ave[i][j] / ( H[i][j] * r[i][j] ) ) )
            alpha_ss_meru2[i].append( 9./(7.*math.pi) * beta_mm * (h_ave[i][j] / ( H[i][j] * r[i][j] ) ) )
            j = j + 1
        f.close()
#
        print np.mean(T[i],r_start:r_limit])
#
    # Set run specific fit start/end locations as per dictionary definitions
#
        r_fit_S = d['r_fitS_'+str(ea_run)] ; r_fit_S = np.array(r_fit_S)
        r_fit_E = d['r_fitE_'+str(ea_run)] ; r_fit_E = np.array(r_fit_E)
        exp_err = d['exp_err_'+str(ea_run)]
#
    r = np.array(r) ; Q = np.array(Q) ; T = np.array(T) ; sig = np.array(sig)
    alpha_ss_lp = np.array(alpha_ss_lp)
    alpha_ss_meru1 = np.array(alpha_ss_meru1)
    alpha_ss_meru2 = np.array(alpha_ss_meru2)
#
    # Define timearr array as list of physical time of snapshot
#
    timearr = []
#
    if (snaparr.ndim == 1):
        for i in range(len(snaparr)):
            timearr.append( (snaparr[i] * 0.01) + 78.68 )
        timearr = np.array(timearr)
    elif (snaparr.ndim == 2):
        for i in range(0,len(snaparr)):
            for j in range(0,len(snaparr[i])):
                timearr.append( (snaparr[i][j] * .01) + 78.68  )
        timearr = np.array(timearr)
        timearr = np.reshape(timearr, (len(snaparr),len(snaparr[i])) )
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # - - - NON-EA RUNS
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
# Begin plotting comparison of simulation time references
#
    if (snaparr_tmp.ndim == 1):
#
    # Plot both Lodato & Price (2010) and Meru+ (2012) formalisms for alpha
#
        fig = plt.figure(1)
        ax1 = plt.subplot(111)
        for i in range(0, len(snaparr_tmp)):
            for j in range(0, len(snapcore_tmp)):
                if (snaparr_tmp[i] == snapcore_tmp[j]):
#
                    rnew = np.arange(1,r_limit,r_limit/(r_limit/smooth))
                    tck1 = interpolate.splrep(r[snaparr_tmp[i],1:r_limit], \
                       alpha_ss_lp[snaparr_tmp[i],1:r_limit])
                    alph1new = interpolate.splev(rnew, tck1, der = 0)
                    line1 = plt.plot(rnew, alph1new, color = col_arr[j], \
                       label = 'Lodato & Price (2010)')
#
#                    tck2 = interpolate.splrep(r[snaparr_tmp[i],1:r_limit], \
#                       alpha_ss_meru1[snaparr_tmp[i],1:r_limit])
#                    alph2new = interpolate.splev(rnew, tck2, der = 0)
#                    line2 = plt.plot(rnew, alph2new, color = col_arr[j], \
#                       label = 'Meru et al. (2012) - Linear', \
#                       linestyle = 'dashed')
#
#                    tck3 = interpolate.splrep(r[snaparr_tmp[i],1:r_limit], \
#                       alpha_ss_meru2[snaparr_tmp[i],1:r_limit])
#                    alph3new = interpolate.splev(rnew, tck3, der = 0)
#                    line1 = plt.plot(rnew, alph3new, color = col_arr[j], \
#                       label = 'Meru et al. (2012) - Quadratic',
#                       linestyle = 'dotted')
#
        plt.xlabel("Radius (AU)", fontsize = 18, labelpad=0.5)
        plt.xticks(fontsize = 15) ; ax1.set_xlim(r_start, 120) ; ax1.set_xscale('log')
        ax1.set_ylim(0.01, 1.)
        plt.ylabel(r'$\alpha_{SS}$', fontsize = 18, labelpad=0.5)
        plt.yscale('log')
#        plt.legend(loc = 'upper left', fontsize = 6)
#
        plt.yticks(fontsize = 15)
        plt.savefig(str(plotdir)+'alpha_r_'+str(r_limit)+'AU.'+str(plot_form), \
           format=str(plot_form), dpi=150)
        plt.clf()
#
#
    # Plot individual r vs. Q distribution
#
        fig = plt.figure(1)
        ax1 = plt.subplot(111)
        for i in range(0, len(snaparr_tmp)):
            for j in range(0, len(snapcore_tmp)):
                if (snaparr_tmp[i] == snapcore_tmp[j]):
#
                    rnew = np.arange(1,r_limit,r_limit/(r_limit/smooth))
                    tck = interpolate.splrep(r[snaparr_tmp[i],1:r_limit], \
                       Q[snaparr_tmp[i],1:r_limit])
                    Qnew = interpolate.splev(rnew, tck, der = 0)
#
                    line1 = plt.plot(rnew, Qnew, color = col_arr[j])
#
        plt.xlabel("Radius (AU)", fontsize = 18, labelpad=0.5)
        plt.xticks(fontsize = 15) ; ax1.set_xlim(r_start, 120) ; ax1.set_xscale('log')
        plt.ylabel('Toomre Q value', fontsize = 18, labelpad=0.5)
#        ax1.axhline(y=1.,color='k',ls='dashed')
#
        if logq is True:
            ax1.set_yscale('log')
            ax1.set_ylim(0.25,4)
        elif logq is False:
            ax1.set_ylim(0, 4)
#
        plt.yticks(fontsize = 15)
        plt.savefig(str(plotdir)+'Q_r_'+str(r_limit)+'AU.'+str(plot_form), \
           format=str(plot_form), dpi=150)
        plt.clf()
#
# Plot surface density vs. radius for non-ea runs with fit to log-log data to find power index, p
#
    # Set y_fitted array sizes dependent on disc radial extents defined in r_fitS and r_fitE
#
        coeffs = [0]*len(snaparr_tmp) ; matcov = [0]*len(snaparr_tmp) ; p_err = [0]*len(snaparr_tmp)
        y_fitted = [[] for i in range(len(snaparr_tmp))]
#
        for i in range(0, len(snaparr_tmp)):
            for j in range(0, len(sig[snaparr_tmp[i],r_fit_S[i]:r_fit_E[i]] ) ):
                y_fitted[i].append(0)
#
        plt.figure(1)
        ax1 = plt.subplot(111)
#
        leg_time = []
        handles = []
#
        for i in range(0, len(snaparr_tmp)):
#
            coeffs[i], matcov[i] = curve_fit(func, np.log10(r[snaparr_tmp[i],r_fit_S[i]:r_fit_E[i]]), \
               np.log10(sig[snaparr_tmp[i],r_fit_S[i]:r_fit_E[i]]), [1, 1])
            p_err[i] = np.sqrt(np.diag(matcov[i]))
#
    # Write exponents to file
#
            f = open(arch_dir+'SD_exps.dat','a')
            f.write(str(ea_run)+' '+str(snaparr[i])+' '+str(round(timearr[i],4))+' ' \
               +str( round(coeffs[i][0], 2) )+' '+str(p_err[i][0])+' '+str(exp_err[i]) \
               +' '+str(np.round(np.mean(H[i][5:25]), 2))+'\n' )
            f.close()
#
            for j in range(0, len(snapcore_tmp)):
                if (snaparr_tmp[i] == snapcore_tmp[j]):\
#
                    leg_time.append(str(round(timearr[i], 3) ) )
#
                    y_fitted[i] = func(np.log10(r[snaparr_tmp[i],r_fit_S[i]:r_fit_E[i]]), \
                        coeffs[i][0], coeffs[i][1])
#
                    rnew = np.arange(r_start,r_limit,r_limit/(r_limit/smooth))
                    tck = interpolate.splrep(r[snaparr_tmp[i],r_start:r_limit], \
                       sig[snaparr_tmp[i],r_start:r_limit])
                    signew = interpolate.splev(rnew, tck, der = 0)
#
                    if (loglin == "TRUE"):
                        line1 = plt.plot(rnew, signew, color = col_arr[j])
                        line2 = plt.plot(r[snaparr_tmp[i],r_fit_S[i]:r_fit_E[i]], \
                           np.power(10.,y_fitted[i]), color = col_arr[j], linewidth = 2, \
                           linestyle = 'dashed', label = str(round(abs(coeffs[i][0]), 1)) )
                        h, l = ax1.get_legend_handles_labels()
                        ax1.set_yscale('log')

                    elif (loglin == "FALSE"):
                        line1 = plt.plot(rnew, signew, color = col_arr[j])
                        line2 = plt.plot(r[snaparr_tmp[i],r_fit_S[i]:r_fit_E[i]], \
                           np.power(10.,y_fitted[i]), color = col_arr[j], linewidth = 2, \
                           linestyle = 'dashed', label = str(round(abs(coeffs[i][0]), 1)) )
                        h, l = ax1.get_legend_handles_labels()
                        ax1.set_xscale('log') ; ax1.set_yscale('log')

        plt.xlabel("Radius (AU)", fontsize = 18, labelpad=0.5); plt.xticks(fontsize = 15)
        plt.ylabel('log '+r'$\Sigma$ (cm$^{-2}$)', fontsize = 18, labelpad=0.5)
        ax1.set_ylim(0, 1550) ; ax1.set_xlim(r_start, 120) ; plt.yticks(fontsize = 15)
        legend1 = plt.legend(h, l, loc='upper right', title = "p values", fontsize=14)
        ax1 = plt.gca().add_artist(legend1)
        if (len(leg_time) == 3):
            plt.legend(h, [leg_time[0], leg_time[1], leg_time[2]], \
               loc = 'lower left', title = "Time (kyr)", fontsize = 14)
        plt.savefig(str(plotdir)+'SD_r_'+str(r_limit)+'AU.'+str(plot_form), \
           format=str(plot_form), dpi=150)
        plt.clf()
#
# Plot temperature vs. radius for non-ea runs with fit to log-log data to find power index, q
#
        coeffs = [0]*len(snaparr_tmp) ; matcov = [0]*len(snaparr_tmp) ; q_err = [0]*len(snaparr_tmp)
        y_fitted = [[] for i in range(len(snaparr_tmp))]
#
        for i in range(0, len(snaparr_tmp)):
            for j in range(0, len(T[snaparr_tmp[i],r_fit_S[i]:r_fit_E[i]] ) ):
                y_fitted[i].append(0)
#
        plt.figure(1)
        ax1 = plt.subplot(111)
        for i in range(0, len(snaparr_tmp)):
#
            coeffs[i], matcov[i] = curve_fit(func, np.log10(r[snaparr_tmp[i],r_fit_S[i]:r_fit_E[i]]), \
               np.log10(T[snaparr_tmp[i],r_fit_S[i]:r_fit_E[i]]), [1, 1])
            q_err[i] = np.sqrt(np.diag(matcov[i]))
            y_fitted[i] = func(np.log10(r[snaparr_tmp[i],r_fit_S[i]:r_fit_E[i]]), \
               coeffs[i][0], coeffs[i][1])
#
    # Write exponents to file
#
            f = open(arch_dir+'T_exps.dat','a')
            f.write(str(ea_run)+' '+str(snaparr[i])+' '+str(round(timearr[i],4))+ \
               ' '+str( round(coeffs[i][0], 2) )+' '+str(q_err[i][0])+'\n' )
            f.close()
#
            for j in range(0, len(snapcore_tmp)):
                if (snaparr_tmp[i] == snapcore_tmp[j]):
#
                    rnew = np.arange(r_start,r_limit,r_limit/(r_limit/smooth))
                    tck = interpolate.splrep(r[snaparr_tmp[i],r_start:r_limit], \
                       T[snaparr_tmp[i],r_start:r_limit])
                    Tnew = interpolate.splev(rnew, tck, der = 0)
#
                    if (loglin == "TRUE"):
                        line1 = plt.plot(rnew, Tnew, color = col_arr[j])
                        line2 = plt.plot(r[snaparr_tmp[i],r_fit_S[i]:r_fit_E[i]], \
                           np.power(10.,y_fitted[i]), color = col_arr[j], linewidth = 2, \
                           linestyle = 'dashed', label = str(round(abs(coeffs[i][0]), 1)) )
                        ax1.set_yscale('log')
                    elif (loglin == "FALSE"):
                        line1 = plt.plot(rnew, Tnew, color = col_arr[j])
                        line2 = plt.plot(r[snaparr_tmp[i],r_fit_S[i]:r_fit_E[i]], \
                           np.power(10,y_fitted[i]), color = col_arr[j], linewidth = 2, \
                           linestyle = 'dashed', label = str(round(abs(coeffs[i][0]), 1)) )
                        ax1.set_xscale('log') ; ax1.set_yscale('log')
#
        plt.xlabel("Radius (AU)", fontsize = 18, labelpad=0.5) ; ax1.set_xlim(r_start, 120)
        plt.xticks(fontsize = 15)
        plt.ylabel("log Temperature (K)", fontsize = 18, labelpad=0.5)
        ax1.set_ylim(0, 550) ; plt.yticks(fontsize = 15)
        plt.legend(loc='upper right', title = "q values", fontsize=14)
        plt.savefig(str(plotdir)+'T_r_'+str(r_limit)+'AU.'+str(plot_form), \
           format=str(plot_form), dpi=150)
        plt.clf()
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # - - - EA RUNS [ VARYING OUTBURST DURATION ]
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
    if (snaparr_tmp.ndim == 2):
#
    # Invoke a switch statement to loop over pdisc analysis of runs with ea to focus on (a) single
    # accretion event and (b) focus on snapshot relative to accretion state for different events
#
        for lendur in range(0, 2):
            title_point = [] ; leg_point = []
#
            if (lendur == 0):
                for i in range(0, len(snaparr_tmp)):
                    title_point.append(EA_lenref[i])
                for i in range(0, len(snaparr_tmp[0])):
                    leg_point.append(EA_timeref[i])
            elif (lendur == 1):
                for i in range(0, len(snaparr_tmp[0])):
                    title_point.append(EA_timeref[i])
                for i in range(0, len(snaparr_tmp)):
                    leg_point.append(EA_lenref[i])
#
    # Rotate snaparr, timearr and planet mass/radius data dependent on switch
#
                snaparr_tmp = np.rot90(snaparr_tmp, 1) ; snaparr_tmp = snaparr_tmp[::-1]
                timearr = np.rot90(timearr, 1) ; timearr = timearr[::-1]
                r_fit_S = np.rot90(r_fit_S, 1) ; r_fit_S = r_fit_S[::-1]
                r_fit_E = np.rot90(r_fit_E, 1) ; r_fit_E = r_fit_E[::-1]
#
    # Now loop over files dependent on switch condition
#
            for i in range(0,len(title_point)):

                leg_time = []
#
#
    # Plot both Lodato & Price (2010) and Meru+ (2012) formalisms for alpha
#
                fig = plt.figure(1)
                ax1 = plt.subplot(111)
                for j in range(0, len(snaparr_tmp[0])):
#
                    rnew = np.arange(1,r_limit,r_limit/(r_limit/smooth))
                    tck1 = interpolate.splrep(r[snaparr_tmp[i][j],1:r_limit], \
                       alpha_ss_lp[snaparr_tmp[i][j],1:r_limit])
                    alph1new = interpolate.splev(rnew, tck1, der = 0)
                    line1 = plt.plot(rnew, alph1new, color = col_arr[j], \
                       label = 'Lodato & Rice (2010)')
#
#                    tck2 = interpolate.splrep(r[snaparr_tmp[i][j],1:r_limit], \
#                       alpha_ss_meru1[snaparr_tmp[i][j],1:r_limit])
#                    alph2new = interpolate.splev(rnew, tck2, der = 0)
#                    line2 = plt.plot(rnew, alph2new, color = col_arr[j], \
#                       label = 'Meru et al. (2012) - Linear', \
#                       linestyle = 'dashed')
#
#                    tck3 = interpolate.splrep(r[snaparr_tmp[i][j],1:r_limit], \
#                       alpha_ss_meru2[snaparr_tmp[i][j],1:r_limit])
#                    alph3new = interpolate.splev(rnew, tck3, der = 0)
#                    line1 = plt.plot(rnew, alph3new, color = col_arr[j], \
#                       label = 'Meru et al. (2012) - Quadratic', \
#                       linestyle = 'dotted')
#
                    plt.xlabel("Radius (AU)", fontsize = 18, labelpad=0.5)
                    plt.xticks(fontsize = 15) ; ax1.set_xlim(r_start, 120) ; ax1.set_xscale('log')
                    ax1.set_ylim(0.01, 1.)
                    plt.ylabel(r'$\alpha_{SS}$', fontsize = 18, labelpad=0.5)
                    plt.yscale('log')
                    plt.yticks(fontsize = 15)
#                    plt.legend(loc = 'upper left', fontsize = 6)
#
                plt.savefig(str(plotdir)+'alpha_r_'+str(r_limit)+'AU_'+title_point[i]+'.'+str(plot_form), \
                   format=str(plot_form), dpi=150)
                plt.clf()
#
#
#
    # Plot individual r vs. Q distribution
#
                fig = plt.figure(1)
                ax1 = plt.subplot(111)
                for j in range(0, len(snaparr_tmp[0])):
#
                    rnew = np.arange(1,r_limit,r_limit/(r_limit/smooth))
                    tck = interpolate.splrep(r[snaparr_tmp[i][j],1:r_limit], \
                       Q[snaparr_tmp[i][j],1:r_limit])
                    Qnew = interpolate.splev(rnew, tck, der = 0)
#
                    line1 = plt.plot(rnew, Qnew, color = col_arr[j])
#
                plt.xlabel("Radius (AU)", fontsize = 18, labelpad=0.5) ; plt.xticks(fontsize = 15)
                ax1.set_xlim(r_start, 120) ; ax1.set_xscale('log')
                plt.ylabel("Toomre Q parameter", fontsize = 18, labelpad=0.5)
#
                if logq is True:
                    ax1.set_yscale('log')
                    ax1.set_ylim(0.25,4)
                elif logq is False:
                    ax1.set_ylim(0, 4)
#
                plt.yticks(fontsize = 15)
#                ax1.axhline(y=1.,color='k',ls='dashed')
                plt.savefig(str(plotdir)+'Q_r_'+str(r_limit)+'AU_'+title_point[i]+'.'+str(plot_form), \
                   format=str(plot_form), dpi=150)
                plt.clf()
#
# Plot surface density vs. radius for ea runs with fit to log-log data to find power index, p
#
                coeffs = [0]*len(snaparr_tmp[0]) ; matcov = [0]*len(snaparr_tmp[0]) ; p_err = [0]*len(snaparr_tmp[0])
                y_fitted = [[] for j in range(len(snaparr_tmp[0]))]
#
                for j in range(0, len(snaparr_tmp[0])):
#
                    for k in range(0, len(sig[snaparr_tmp[i][j],r_fit_S[i][j]:r_fit_E[i][j]]) ):
                        y_fitted[j].append(0)
#
                plt.figure(1)
                ax1 = plt.subplot(111)
#
                for j in range(0, len(snaparr_tmp[0])):
#
                    leg_time.append(str(round(timearr[i][j], 3)) )
#
                    coeffs[j], matcov[j] = curve_fit(func, np.log10(r[snaparr_tmp[i][j],r_fit_S[i][j]:r_fit_E[i][j]]), \
                       np.log10(sig[snaparr_tmp[i][j],r_fit_S[i][j]:r_fit_E[i][j]]), [1, 1])
                    p_err[j] = np.sqrt(np.diag(matcov[j]))
                    y_fitted[j] = func(np.log10(r[snaparr_tmp[i][j],r_fit_S[i][j]:r_fit_E[i][j]]), \
                       coeffs[j][0], coeffs[j][1])
#
    # Write exponents to file
#
                    if( lendur == 0 ):
                        f = open(arch_dir+'SD_exps.dat','a')
                        f.write(str(ea_run)+' '+str(snaparr[i][j])+' '+str(round(timearr[i][j],4))+ \
                           ' '+str( round(coeffs[j][0], 2) )+' '+str(p_err[j][0])+' '+str(exp_err[i][j]) \
                           +' '+str(np.round(np.mean(H[i][2:25]), 2))+'\n' )
                        f.close()
#
                    rnew = np.arange(r_start,r_limit,r_limit/(r_limit/smooth))
                    tck = interpolate.splrep(r[snaparr_tmp[i][j],r_start:r_limit], \
                       sig[snaparr_tmp[i][j],r_start:r_limit])
                    signew = interpolate.splev(rnew, tck, der = 0)
#
                    if (loglin == "TRUE"):
                        line1 = plt.plot(rnew, signew, color = col_arr[j])
                        line2 = plt.plot(r[snaparr_tmp[i][j],r_fit_S[i][j]:r_fit_E[i][j]], \
                           np.power(10.,y_fitted[j]), color = col_arr[j], linewidth = 2, \
                           linestyle = 'dashed', \
                           label = str(leg_point[j])+': '+str(round(abs(coeffs[j][0]), 1))+ \
                              ' - '+str(round(abs(timearr[i][j]), 1))+' kyr' )
                        h, l = ax1.get_legend_handles_labels()
                        ax1.set_yscale('log')
                    elif (loglin == "FALSE"):
                        line1 = plt.plot(rnew, signew, color = col_arr[j])
                        line2 = plt.plot(r[snaparr_tmp[i][j],r_fit_S[i][j]:r_fit_E[i][j]], \
                           np.power(10.,y_fitted[j]), color = col_arr[j], linewidth = 2, \
                           linestyle = 'dashed', \
                           label = str(leg_point[j])+': '+str(round(abs(coeffs[j][0]), 1)) )
                        h, l = ax1.get_legend_handles_labels()
                        ax1.set_xscale('log') ; ax1.set_yscale('log')
#
                plt.xlabel("Radius (AU)", fontsize = 18, labelpad=0.5) ; ax1.set_xlim(r_start, 120)
                plt.xticks(fontsize = 15)
                plt.ylabel('log '+r'$\Sigma$ (cm$^{-2}$)', fontsize = 18, labelpad=0.5)
                ax1.set_ylim(0, 1250) ; plt.yticks(fontsize = 15)
                legend1 = plt.legend(loc='upper right', title = "p values", fontsize=14)
                ax1 = plt.gca().add_artist(legend1)
                if (len(leg_time) == 4):
                    plt.legend(h, [leg_time[0], leg_time[1], leg_time[2], leg_time[3]], \
                      loc = 'lower left', title = "Time (kyr)", fontsize = 14)
                plt.savefig(str(plotdir)+'SD_r_'+str(r_limit)+'AU_'+str(title_point[i])+'.'+str(plot_form), \
                   format=str(plot_form), dpi=150)
                plt.clf()
#
# Plot temperature vs. radius for ea runs with fit to log-log data to find power index, q
#
                coeffs = [0]*len(snaparr_tmp[0]) ; matcov = [0]*len(snaparr_tmp[0]) ; q_err = [0]*len(snaparr_tmp[0])
                y_fitted = [[] for j in range(len(snaparr_tmp[0]))]
#
                for j in range(0, len(snaparr_tmp[0])):
                    for k in range(0, len(T[snaparr_tmp[i][j],r_fit_S[i][j]:r_fit_E[i][j]]) ):
                        y_fitted[j].append(0)
#
                plt.figure(1)
                ax1 = plt.subplot(111)
#
                for j in range(0, len(snaparr_tmp[0])):
#
                    coeffs[j], matcov[j] = curve_fit(func, np.log10(r[snaparr_tmp[i][j],r_fit_S[i][j]:r_fit_E[i][j]]), \
                       np.log10(T[snaparr_tmp[i][j],r_fit_S[i][j]:r_fit_E[i][j]]), [1, 1])
                    q_err[j] = np.sqrt(np.diag(matcov[j]))
                    y_fitted[j] = func(np.log10(r[snaparr_tmp[i][j],r_fit_S[i][j]:r_fit_E[i][j]]),
                       coeffs[j][0], coeffs[j][1])
#
    # Write exponents to file
#
                    if( lendur == 0 ):
                        f = open(arch_dir+'T_exps.dat','a')
                        f.write(str(ea_run)+' '+str(snaparr[i][j])+' '+str(round(timearr[i][j],4))+ \
                           ' '+str( round(coeffs[j][0], 2) )+' '+str(q_err[j][0])+'\n' )
                        f.close()
#
                    rnew = np.arange(r_start,r_limit,r_limit/(r_limit/smooth))
                    tck = interpolate.splrep(r[snaparr_tmp[i][j],r_start:r_limit], \
                       T[snaparr_tmp[i][j],r_start:r_limit])
                    Tnew = interpolate.splev(rnew, tck, der = 0)
#
                    if (loglin == "TRUE"):
                        line1 = plt.plot(rnew, Tnew, color = col_arr[j])
                        line2 = plt.plot(r[snaparr_tmp[i][j],r_fit_S[i][j]:r_fit_E[i][j]], \
                           np.power(10.,y_fitted[j]), color = col_arr[j], linewidth = 2, \
                           linestyle = 'dashed', \
                           label = str(leg_point[j])+': '+str(round(abs(coeffs[j][0]), 1)) )
                        ax1.set_yscale('log')
                    elif (loglin == "FALSE"):
                        line1 = plt.plot(rnew, Tnew, color = col_arr[j])
                        line2 = plt.plot(r[snaparr_tmp[i][j],r_fit_S[i][j]:r_fit_E[i][j]], \
                           np.power(10.,y_fitted[j]), color = col_arr[j], linewidth = 2, \
                           linestyle = 'dashed', \
                           label = str(leg_point[j])+': '+str(round(abs(coeffs[j][0]), 1)) )
                        ax1.set_xscale('log') ; ax1.set_yscale('log')
#
                plt.xlabel("Radius (AU)", fontsize = 18, labelpad=0.5) ; ax1.set_xlim(r_start, 120)
                plt.xticks(fontsize = 15)
                plt.ylabel("log Temperature (K)", fontsize = 18, labelpad=0.5) ; ax1.set_ylim(0, 350)
                plt.yticks(fontsize = 15)
                plt.legend(loc='upper right', title = "q values", fontsize=14)
                plt.savefig(str(plotdir)+'T_r_'+str(r_limit)+'AU_'+str(title_point[i])+'.'+str(plot_form), \
                   format='eps', dpi=150)
                plt.clf()
#
    # Now generate exponent fits to the inbetween accretion (IA) snapshots for respective EF run
    # Values need to be read in similar to manner of non IA snapshots, however plotting routine not invoked
#
    if ( snaparr.ndim == 2 ):
#
        if (snaparr_IA.ndim == 1):
            file_n = len(snaparr_IA)
        elif (snaparr_IA.ndim == 2):
            file_n = len(snaparr_IA)*len(snaparr_IA[0])
#
        time = [[] for i in range(file_n)] ; r = [[] for i in range(file_n)]
        T = [[] for i in range(file_n)] ; sig = [[] for i in range(file_n)]
#
        file_list = []
        if (snaparr_IA.ndim == 1):
            snaparr_tmp_IA = np.array([0]*len(snaparr_IA))
            fcount = 0
            for a in range(0, len(snaparr_IA)):
                if (snaparr_IA[a] < (1000-70)):
                    file_list.append(dat_dir+'pdisc/DE05.du.00'+ \
                       str(snaparr_IA[a]+69)+'.pdisc.1')
                elif (snaparr_IA[a] > (1000-70)):
                    file_list.append(dat_dir+'pdisc/DE05.du.0'+ \
                       str(snaparr_IA[a]+69)+'.pdisc.1')
                snaparr_tmp_IA[a] = fcount
                fcount = fcount + 1
#
        elif (snaparr_IA.ndim == 2):
            snaparr_tmp_IA = np.array([[0]*len(snaparr_IA[0])]*len(snaparr_IA))
            fcount = 0
            for a in range(0, len(snaparr_IA)):
                for b in range(0, len(snaparr_IA[0])):
                    if (snaparr_IA[a][b] < (1000-70)):
                        file_list.append(dat_dir+'pdisc/DE05.du.00'+ \
                           str(snaparr_IA[a][b]+69)+'.pdisc.1')
                    elif (snaparr_IA[a][b] > (1000-70)):
                        file_list.append(dat_dir+'pdisc/DE05.du.0'+ \
                           str(snaparr_IA[a][b]+69)+'.pdisc.1')
                    snaparr_tmp_IA[a][b] = fcount
                    fcount = fcount + 1
#
        for i in range(0, file_n):
            f = open(file_list[i], 'r')
            header = f.readline()
#
            for line in f:
                line = line.strip() ; columns = line.split()
                time[i].append(float(columns[0])/1000.) ; r[i].append(float(columns[1]))
                T[i].append(float(columns[3])) ; sig[i].append(float(columns[4]))
            f.close()
#
            r_fit_S_IA = d['r_fitS_'+str(ea_run)+'_IA'] ; r_fit_S_IA = np.array(r_fit_S_IA)
            r_fit_E_IA = d['r_fitE_'+str(ea_run)+'_IA'] ; r_fit_E_IA = np.array(r_fit_E_IA)
            exp_err_IA = d['exp_err_'+str(ea_run)+'_IA']
#
        time = np.array(time) ; r = np.array(r) ; T = np.array(T) ; sig = np.array(sig)
#
    # Now compute exponents
#
        if (snaparr_IA.ndim == 1):
#
            coeffs = [0]*len(snaparr_IA) ; matcov = [0]*len(snaparr_IA)
            p_err = [0]*len(snaparr_IA) ; q_err = [0]*len(snaparr_IA)
            y_fitted = [[] for i in range(len(snaparr_IA))]
#
            for i in range(0, len(snaparr_IA)):
#
    # For SD
#
                coeffs[i], matcov[i] = curve_fit(func, np.log10(r[snaparr_tmp_IA[i],r_fit_S_IA[i]:r_fit_E_IA[i]]), \
                   np.log10(sig[snaparr_tmp_IA[i],r_fit_S_IA[i]:r_fit_E_IA[i]]), [1, 1])
                p_err[i] = np.sqrt(np.diag(matcov[i]))
                y_fitted[i] = func(np.log10(r[snaparr_tmp_IA[i],r_fit_S_IA[i]:r_fit_E_IA[i]]), \
                   coeffs[i][0], coeffs[i][1])
#
                f = open(arch_dir+'SD_exps.dat','a')
                f.write(str(ea_run)+' '+str(snaparr_IA[i])+' '+str(time[snaparr_tmp_IA[i],0])+ \
                   ' '+str( round(coeffs[i][0], 2) )+' '+str(p_err[i][0])+' '+str(exp_err_IA[i]) \
                   +' '+str(np.round(np.mean(H[i][2:25]), 2))+'\n' )
                f.close()
#
    # For T
#
                coeffs[i], matcov[i] = curve_fit(func, np.log10(r[snaparr_tmp_IA[i],r_fit_S_IA[i]:r_fit_E_IA[i]]), \
                   np.log10(T[snaparr_tmp_IA[i],r_fit_S_IA[i]:r_fit_E_IA[i]]), [1, 1])
                q_err[i] = np.sqrt(np.diag(matcov[i]))
                y_fitted[i] = func(np.log10(r[snaparr_tmp_IA[i],r_fit_S_IA[i]:r_fit_E_IA[i]]), \
                   coeffs[i][0], coeffs[i][1])
#
                f = open(arch_dir+'T_exps.dat','a')
                f.write(str(ea_run)+' '+str(snaparr_IA[i])+' '+str(time[snaparr_tmp_IA[i],0])+ \
                   ' '+str( round(coeffs[i][0], 2) )+' '+str(q_err[i][0])+'\n' )
                f.close()
#
    # Now for case where > 2 outbursts are analysed
#
        if (snaparr_IA.ndim == 2):
            file_n = len(snaparr_IA)*len(snaparr_IA[0])
#
            for i in range(0, len(snaparr_IA)):
#
                coeffs = [0]*len(snaparr_tmp_IA[0]) ; matcov = [0]*len(snaparr_tmp_IA[0])
                p_err = [0]*len(snaparr_tmp_IA[0]) ; q_err = [0]*len(snaparr_tmp_IA[0])
                y_fitted = [[] for j in range(len(snaparr_tmp_IA[0]))]
#
                for j in range(0, len(snaparr_IA[0])):
#
                    coeffs[j], matcov[j] = curve_fit(func, np.log10(r[snaparr_tmp_IA[i][j],r_fit_S_IA[i][j]:r_fit_E_IA[i][j]]), \
                       np.log10(sig[snaparr_tmp_IA[i][j],r_fit_S_IA[i][j]:r_fit_E_IA[i][j]]), [1, 1])
                    p_err[j] = np.sqrt(np.diag(matcov[j]))
                    y_fitted[j] = func(np.log10(r[snaparr_tmp_IA[i][j],r_fit_S_IA[i][j]:r_fit_E_IA[i][j]]), \
                       coeffs[j][0], coeffs[j][1])
#
                    f = open(arch_dir+'SD_exps.dat','a')
                    f.write(str(ea_run)+' '+str(snaparr_IA[i][j])+' '+str(time[snaparr_tmp_IA[i][j],0])+ \
                       ' '+str( round(coeffs[j][0], 2) )+' '+str(p_err[j][0])+' '+str(exp_err_IA[i][j]) \
                       +' '+str(np.round(np.mean(H[i][2:25]), 2))+'\n' )
                    f.close()
#
                    coeffs[j], matcov[j] = curve_fit(func, np.log10(r[snaparr_tmp_IA[i][j],r_fit_S_IA[i][j]:r_fit_E_IA[i][j]]), \
                       np.log10(T[snaparr_tmp_IA[i][j],r_fit_S_IA[i][j]:r_fit_E_IA[i][j]]), [1, 1])
                    q_err[j] = np.sqrt(np.diag(matcov[j]))
                    y_fitted[j] = func(np.log10(r[snaparr_tmp_IA[i][j],r_fit_S_IA[i][j]:r_fit_E_IA[i][j]]), \
                       coeffs[j][0], coeffs[j][1])
#
                    f = open(arch_dir+'T_exps.dat','a')
                    f.write(str(ea_run)+' '+str(snaparr_IA[i][j])+' '+str(time[snaparr_tmp_IA[i][j],0])+ \
                       ' '+str( round(coeffs[j][0], 2) )+' '+str(q_err[j][0])+'\n' )
                    f.close()
