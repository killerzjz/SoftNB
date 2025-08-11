import numpy as np
from numpy import fft
from numpy import *
import matplotlib.pyplot as plt
import scipy.io as scio
import scipy.signal as signal
import pandas as pd
from scipy import interpolate
from utils import *



xlabel_font = {
    'fontsize': 30,
    'fontweight': 'light',
    'color': 'black',
}

ylabel_font = {
    'fontsize': 30,
    'fontweight': 'light',
    'color': 'black',
}

if __name__=='__main__':
    # find npss start
    data = scio.loadmat('file/nb_signal_DL.mat')
    nb_sig = np.array(data['nb_signal_DL'])
    data = scio.loadmat('file/ref_grid.mat')
    ref_grid = np.array(data['ref_grid'])
    data = scio.loadmat('file/npss.mat')
    npss_standard = np.array(data['npss'])
    npss_standard = npss_standard.reshape(1920)
    cfo = 2000
    nb_sig = nb_sig.reshape(38400)
    at_cor = np.zeros((int(19200*1.8), 1), dtype=complex)
    at_cor5db = np.zeros((int(19200*1.8), 1), dtype=complex)
    at_cor_5db = np.zeros((int(19200*1.8), 1), dtype=complex)
    at_cor_15db = np.zeros((int(19200*1.8), 1), dtype=complex)
    at_cor_25db = np.zeros((int(19200*1.8), 1), dtype=complex)
    nb_sig5db = nb_sig + wgn(nb_sig,-1)
    nb_sig_5db = nb_sig + wgn(nb_sig,-5)
    nb_sig_15db = nb_sig + wgn(nb_sig,-15)
    nb_sig_25db = nb_sig + wgn(nb_sig,-30)
    nb_sig = nb_sig*[np.exp(1j*2*np.pi*idx*cfo/1920000) for idx in range(38400)]
    nb_sig5db = nb_sig5db*[np.exp(1j*2*np.pi*idx*cfo/1920000) for idx in range(38400)]
    nb_sig_5db = nb_sig_5db*[np.exp(1j*2*np.pi*idx*cfo/1920000) for idx in range(38400)]
    nb_sig_15db = nb_sig_15db*[np.exp(1j*2*np.pi*idx*cfo/1920000) for idx in range(38400)]
    for i in range(int(19200*1.8)):
        at_cor[i] = auto_correlation(nb_sig[i:], 3, npss_zc_seq, 137)
    for i in range(int(19200*1.8)):
        at_cor5db[i] = auto_correlation(nb_sig5db[i:], 3, npss_zc_seq, 137)
    for i in range(int(19200*1.8)):
        at_cor_5db[i] = auto_correlation(nb_sig_5db[i:], 3, npss_zc_seq, 137)
    for i in range(int(19200*1.8)):
        at_cor_15db[i] = auto_correlation(nb_sig_15db[i:], 3, npss_zc_seq, 137)
    for i in range(int(19200*1.8)):
        at_cor_25db[i] = auto_correlation(nb_sig_25db[i:], 3, npss_zc_seq, 137)

    plt.plot(abs(at_cor[0::8]))
    ax = plt.gca()
    ax.ticklabel_format(style='sci', scilimits=(-1,2), axis='x')
    plt.rc('font',size=16)
    plt.xlabel('Time offset(samples)',fontdict={'family':'Times New Roman', 'size':20,'weight':'bold'})
    plt.ylabel('correlation level',fontdict={'family':'Times New Roman', 'size':20,'weight':'bold'})
    plt.tick_params(labelsize=18)
    plt.show()

    plt.figure(figsize=(7,7))
    plt.plot(abs(at_cor5db[0::8]))
    ax = plt.gca()
    ax.ticklabel_format(style='sci', scilimits=(-1,2), axis='x')
    plt.xlabel('Time offset(samples)',fontdict={'family':'Times New Roman', 'size':24,'weight':'bold'})
    plt.ylabel('Correlation level',fontdict={'family':'Times New Roman', 'size':24,'weight':'bold'})
    plt.yticks(np.arange(0,0.4,0.1))
    # plt.title("SNR = 5dB",fontdict={'family':'Times New Roman', 'size':20,'weight':'bold'})
    plt.tick_params(labelsize=18)
    plt.show()

    plt.plot(abs(at_cor_5db[0::8]))
    ax = plt.gca()
    ax.ticklabel_format(style='sci', scilimits=(-1,2), axis='x')
    plt.xlabel('Time offset(samples)',fontdict={'family':'Times New Roman', 'size':20,'weight':'bold'})
    plt.ylabel('Correlation level',fontdict={'family':'Times New Roman', 'size':20,'weight':'bold'})
    plt.yticks(np.arange(0,0.4,0.1))
    plt.title("SNR = -5dB",fontdict={'family':'Times New Roman', 'size':20,'weight':'bold'})
    plt.tick_params(labelsize=18)
    plt.show()

    plt.plot(abs(at_cor_15db[0::8]))
    ax = plt.gca()
    ax.ticklabel_format(style='sci', scilimits=(-1,2), axis='x')
    plt.xlabel('Time offset(samples)',fontdict={'family':'Times New Roman', 'size':20,'weight':'bold'})
    plt.ylabel('Correlation level',fontdict={'family':'Times New Roman', 'size':20,'weight':'bold'})
    plt.yticks(np.arange(0,0.4,0.1))
    plt.title("SNR = -15dB",fontdict={'family':'Times New Roman', 'size':20,'weight':'bold'})
    plt.tick_params(labelsize=18)
    plt.show()

    plt.plot(abs(at_cor_25db[0::8]))
    ax = plt.gca()
    ax.ticklabel_format(style='sci', scilimits=(-1,2), axis='x')
    plt.yticks(np.arange(0,0.8,0.2))
    plt.xlabel('Time offset(samples)',fontdict={'family':'Times New Roman', 'size':20,'weight':'bold'})
    plt.ylabel('Correlation level',fontdict={'family':'Times New Roman', 'size':20,'weight':'bold'})
    plt.title("SNR = -25dB",fontdict={'family':'Times New Roman', 'size':20,'weight':'bold'})
    plt.tick_params(labelsize=18)
    plt.show()
    

    npss_start = 9600
    cfo_est = cfo_estimation(nb_sig[0:19200],npss_standard,100)
    nb_sig_t = nb_sig
    nb_sig = nb_sig*[np.exp(1j*idx*cfo_est) for idx in range(38400)]
    print(cfo_est)

    cfo_est5db = cfo_estimation(nb_sig5db[0:19200],npss_standard,5)
    nb_sig5db = nb_sig5db*[np.exp(1j*idx*cfo_est5db) for idx in range(38400)]
    print(cfo_est5db)

    cfo_est_5db = cfo_estimation(nb_sig_5db[0:19200],npss_standard,-5)
    nb_sig_5db = nb_sig_5db*[np.exp(1j*idx*cfo_est_5db) for idx in range(38400)]
    print(cfo_est_5db)

    cfo_est_15db = cfo_estimation(nb_sig_15db[0:19200],npss_standard,-15)
    nb_sig_15db = nb_sig_15db*[np.exp(1j*idx*cfo_est_15db) for idx in range(38400)]
    print(cfo_est_15db)

    # get 2 subframes
    standard_nsss = np.zeros((504,4,132,1),dtype=complex)
    for cell_id in range(504):
        for frame_number in range(4):
            standard_nsss[cell_id][frame_number] = generate_standerd_nsss(frame_number, cell_id)
    f1 = np.array(nb_sig[0:19200])
    f1 = f1.reshape(19200)
    f2 = np.array(nb_sig[19200:38400])
    f2 = f2.reshape(19200)
    # get nsss_seq
    nsss1 = f1[9*1920:10*1920]
    nsss15db = nsss1+wgn(nsss1,5)
    nsss1_5db = nsss1+wgn(nsss1,-5)
    nsss1_15db = nsss1+wgn(nsss1,-15)
    nsss1_25db = nsss1+wgn(nsss1,-25)
    nsss1 = subframe_remove_cp(nsss1)
    nsss15db = subframe_remove_cp(nsss15db)
    nsss1_5db = subframe_remove_cp(nsss1_5db)
    nsss1_15db = subframe_remove_cp(nsss1_15db)
    nsss1_25db = subframe_remove_cp(nsss1_25db)
    # nsss_seq_1 = ofdm_decode(nsss1)
    nsss_seq = ofdm_decode(nsss1,ref_grid,1)
    nsss_seq_15db = ofdm_decode(nsss15db,ref_grid,1)
    nsss_seq_1_5db = ofdm_decode(nsss1_5db,ref_grid,1)
    nsss_seq_1_15db = ofdm_decode(nsss1_15db,ref_grid,1)
    nsss_seq_1_25db = ofdm_decode(nsss1_25db,ref_grid,1)
    nsss2 = f2[9*1920:10*1920]
    nsss2 = subframe_remove_cp(nsss2)
    nsss_seq_2 = ofdm_decode(nsss2,ref_grid,1)
    

    cr_cor = np.zeros((504*4,1),dtype=complex)
    cr_cor5db = np.zeros((504*4,1),dtype=complex)
    cr_cor_5db = np.zeros((504*4,1),dtype=complex)
    cr_cor_15db = np.zeros((504*4,1),dtype=complex)
    cr_cor_25db = np.zeros((504*4,1),dtype=complex)
    max_cell_id = 0
    max_frame_num = 0
    max_cr_cor = 0
    for cell_id in range(504):
        for frame_number in range(4):
            cr_cor[4*(cell_id)+frame_number] = abs(cross_correlation(standard_nsss[cell_id][frame_number], nsss_seq[36:168], 132))
            if cr_cor[4*(cell_id)+frame_number] > max_cr_cor:
                max_cell_id = cell_id
                max_frame_num = frame_number
                max_cr_cor = cr_cor[4*(cell_id)+frame_number]
    max_cr_cor = 0
    for cell_id in range(504):
        for frame_number in range(4):
            cr_cor5db[4*(cell_id)+frame_number] = abs(cross_correlation(standard_nsss[cell_id][frame_number], nsss_seq_15db[36:168], 132))
            if cr_cor5db[4*(cell_id)+frame_number] > max_cr_cor:
                max_cell_id = cell_id
                max_frame_num = frame_number
                max_cr_cor = cr_cor5db[4*(cell_id)+frame_number]
    max_cr_cor = 0
    for cell_id in range(504):
        for frame_number in range(4):
            cr_cor_5db[4*(cell_id)+frame_number] = abs(cross_correlation(standard_nsss[cell_id][frame_number], nsss_seq_1_5db[36:168], 132))
            if cr_cor_5db[4*(cell_id)+frame_number] > max_cr_cor:
                max_cell_id = cell_id
                max_frame_num = frame_number
                max_cr_cor = cr_cor_5db[4*(cell_id)+frame_number]
    max_cr_cor = 0
    for cell_id in range(504):
        for frame_number in range(4):
            cr_cor_15db[4*(cell_id)+frame_number] = abs(cross_correlation(standard_nsss[cell_id][frame_number], nsss_seq_1_15db[36:168], 132))
            if cr_cor_15db[4*(cell_id)+frame_number] > max_cr_cor:
                max_cell_id = cell_id
                max_frame_num = frame_number
                max_cr_cor = cr_cor_15db[4*(cell_id)+frame_number]
    max_cr_cor = 0
    for cell_id in range(504):
        for frame_number in range(4):
            cr_cor_25db[4*(cell_id)+frame_number] = abs(cross_correlation(standard_nsss[cell_id][frame_number], nsss_seq_1_25db[36:168], 132))
            if cr_cor_25db[4*(cell_id)+frame_number] > max_cr_cor:
                max_cell_id = cell_id
                max_frame_num = frame_number
                max_cr_cor = cr_cor_25db[4*(cell_id)+frame_number]
    max_cr_cor = 0

    plt.figure(figsize=(8,6))
    plt.plot(cr_cor)
    plt.xlabel('Generated NSSS ID',fontdict={'family':'Times New Roman', 'size':20,'weight':'bold'})
    plt.ylabel('Correlation level',fontdict={'family':'Times New Roman', 'size':20,'weight':'bold'})
    #plt.title("SNR = 5dB",fontdict={'family':'Times New Roman', 'size':20,'weight':'bold'})
    plt.tick_params(labelsize=18)
    plt.show()
    plt.figure(figsize=(8,6))
    plt.plot(cr_cor5db)
    plt.xlabel('Generated NSSS ID',fontdict={'family':'Times New Roman', 'size':20,'weight':'bold'})
    plt.ylabel('Correlation level',fontdict={'family':'Times New Roman', 'size':20,'weight':'bold'})
    #plt.title("SNR = 5dB",fontdict={'family':'Times New Roman', 'size':20,'weight':'bold'})
    plt.tick_params(labelsize=18)
    plt.show()
    # plt.subplot(424)
    plt.figure(figsize=(8,6))
    plt.plot(cr_cor_5db)
    plt.xlabel('Generated NSSS ID',fontdict={'family':'Times New Roman', 'size':20,'weight':'bold'})
    plt.ylabel('Correlation level',fontdict={'family':'Times New Roman', 'size':20,'weight':'bold'})
    #plt.title("SNR = -5dB",fontdict={'family':'Times New Roman', 'size':20,'weight':'bold'})
    plt.tick_params(labelsize=18)
    plt.show()
    # plt.subplot(426)
    plt.figure(figsize=(8,6))
    plt.plot(cr_cor_15db)
    plt.xlabel('Generated NSSS ID',fontdict={'family':'Times New Roman', 'size':20,'weight':'bold'})
    plt.ylabel('Correlation level',fontdict={'family':'Times New Roman', 'size':20,'weight':'bold'})
    #plt.title("SNR = -15dB",fontdict={'family':'Times New Roman', 'size':20,'weight':'bold'})
    plt.tick_params(labelsize=18)
    plt.show()
    # plt.subplot(428)
    plt.figure(figsize=(8,6))
    plt.plot(cr_cor_25db)
    plt.xlabel('Generated NSSS ID',fontdict={'family':'Times New Roman', 'size':20,'weight':'bold'})
    plt.ylabel('Correlation level',fontdict={'family':'Times New Roman', 'size':20,'weight':'bold'})
    #plt.title("SNR = -25dB",fontdict={'family':'Times New Roman', 'size':20,'weight':'bold'})
    plt.tick_params(labelsize=18)
    # fig.tight_layout(h_pad=-2)
    plt.show()

    print("cell_id:{},max_cr_cor:{}".format(max_cell_id,max_cr_cor))
    plt.title("real NSSS and standard NSSS")
    plt.scatter(real(nsss_seq[36:168]), imag(nsss_seq[36:168]))
    plt.scatter(real(standard_nsss[max_cell_id][max_frame_num]), imag(standard_nsss[max_cell_id][max_frame_num]))
    plt.show()
    npbch1 = nb_sig_25db[0*1920:1*1920]
    # npbch2 = f2[5*1920:6*1920]
    npbch1 = subframe_remove_cp(npbch1)
    npbch1_seq = ofdm_decode(npbch1,ref_grid[0:12,0:14],0)
    npbch1_seq,nrs1_seq = npbch_seq_process(npbch1_seq)
    plt.scatter(real(npbch1_seq), imag(npbch1_seq))
    npbch1_seq = ofdm_decode(npbch1,ref_grid[0:12,0:14],1)
    npbch1_seq,nrs1_seq = npbch_seq_process(npbch1_seq)
    plt.scatter(real(npbch1_seq), imag(npbch1_seq))
    plt.show()

    npdsch1 = nb_sig5db[1*1920:2*1920]
    # npdsch1 = npdsch1+wgn(npdsch1,-15)
    npdsch1_seq = subframe_remove_cp(npdsch1)
    npdsch1_sym = ofdm_decode(npdsch1_seq,ref_grid[0:12,14:28],0)
    npdsch1_sym = npdsch_seq_process(npdsch1_sym)
    plt.scatter(real(npdsch1_sym), imag(npdsch1_sym))
    plt.show()




    
    






