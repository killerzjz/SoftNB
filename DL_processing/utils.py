import numpy as np
from numpy import fft
from scipy import interpolate
from numpy import *
import scipy.signal as signal

NPSS_SYMBOLS = 11
npss_zc_seq = [1,1,1,1,-1,-1,1,1,1,-1,1]
max_magnitude = 0
nrs_y = [5,6,5,6,12,13,12,13]
nrs_x = [0,3,6,9,0,3,6,9]
nrs = [nrs_x,nrs_y]

def npbch_seq_process(npbch_seq):
    npbch_seq_real = np.concatenate((npbch_seq[36:48],
                        npbch_seq[49:51],
                        npbch_seq[52:54],
                        npbch_seq[55:57],
                        npbch_seq[58:60],
                        npbch_seq[61:63],
                        npbch_seq[64:66],
                        npbch_seq[67:69],
                        npbch_seq[70:72],
                        npbch_seq[73:75],
                        npbch_seq[76:78],
                        npbch_seq[79:81],
                        npbch_seq[82:84],
                        npbch_seq[85:87],
                        npbch_seq[88:90],
                        npbch_seq[91:93],
                        npbch_seq[94:96],
                        npbch_seq[97:99],
                        npbch_seq[100:102],
                        npbch_seq[103:105],
                        npbch_seq[106:108],
                        npbch_seq[108:132],
                        npbch_seq[133:135],
                        npbch_seq[136:138],
                        npbch_seq[139:141],
                        npbch_seq[142:144],
                        npbch_seq[145:147],
                        npbch_seq[148:150],
                        npbch_seq[151:153],
                        npbch_seq[154:156],
                        npbch_seq[157:159],
                        npbch_seq[160:162],
                        npbch_seq[163:165],
                        npbch_seq[165:168],),axis=0)
    nrs_seq = np.concatenate((
        npbch_seq[60],
        npbch_seq[63],
        npbch_seq[66],
        npbch_seq[69],
        npbch_seq[72],
        npbch_seq[75],
        npbch_seq[78],
        npbch_seq[81],
        npbch_seq[144],
        npbch_seq[147],
        npbch_seq[150],
        npbch_seq[153],
        npbch_seq[156],
        npbch_seq[159],
        npbch_seq[162],
        npbch_seq[165],
    ),axis = 0)
    return npbch_seq_real,nrs_seq

def npdsch_seq_process(npdsch_seq):
    npdsch_seq_real = np.concatenate((npdsch_seq[0:60],
                                      npdsch_seq[61:63],
                                      npdsch_seq[64:66],
                                      npdsch_seq[67:69],
                                      npdsch_seq[70:72],
                                      npdsch_seq[73:75],
                                      npdsch_seq[76:78],
                                      npdsch_seq[79:81],
                                      npdsch_seq[82:144],
                                      npdsch_seq[145:147],
                                      npdsch_seq[148:150],
                                      npdsch_seq[151:153],
                                      npdsch_seq[154:156],
                                      npdsch_seq[157:159],
                                      npdsch_seq[160:162],
                                      npdsch_seq[163:165],
                                      npdsch_seq[166:168]
                                      ))
    return npdsch_seq_real

def cross_correlation(input1,input2,length):
    seq1 = input1[0:length]
    seq2 = input2[0:length]
    seq = seq1 * seq2.conjugate()
    cr_cor = sum(seq)
    return cr_cor

def auto_correlation(input_data, lag, zc_seq, length):
    at_cor = 0
    for i in range(NPSS_SYMBOLS-lag):
        seq1 = np.array(input_data[i*length:(i+1)*length])
        seq2 = np.array(input_data[(i+lag)*length:(i+lag+1)*length])
        seq = seq1 * seq2.conjugate()
        tmp_sum = sum(seq)
        if zc_seq[i] == zc_seq[i+lag]:
            at_cor = at_cor + tmp_sum
        else:
            at_cor = at_cor - tmp_sum
    return at_cor

def generate_standerd_nsss(frame_number, cell_id):
    nsss_seq = np.zeros((132,1),dtype=complex)
    b = [
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        [1,-1,-1, 1,-1, 1, 1,-1,-1, 1, 1,-1, 1,-1,-1, 1,-1, 1, 1,-1, 1,-1,-1, 1, 1,-1,-1, 1,-1, 1, 1,-1,
          1,-1,-1, 1,-1, 1, 1,-1,-1, 1, 1,-1, 1,-1,-1, 1,-1, 1, 1,-1, 1,-1,-1, 1, 1,-1,-1, 1,-1, 1, 1,-1,
          1,-1,-1, 1,-1, 1, 1,-1,-1, 1, 1,-1, 1,-1,-1, 1,-1, 1, 1,-1, 1,-1,-1, 1, 1,-1,-1, 1,-1, 1, 1,-1,
          1,-1,-1, 1,-1, 1, 1,-1,-1, 1, 1,-1, 1,-1,-1, 1,-1, 1, 1,-1, 1,-1,-1, 1, 1,-1,-1, 1,-1, 1, 1,-1],
        [1,-1,-1, 1,-1, 1, 1,-1,-1, 1, 1,-1, 1,-1,-1, 1,-1, 1, 1,-1, 1,-1,-1, 1, 1,-1,-1, 1,-1, 1, 1,-1,
         -1, 1, 1,-1, 1,-1,-1, 1, 1,-1,-1, 1,-1, 1, 1,-1, 1,-1,-1, 1,-1, 1, 1,-1,-1, 1, 1,-1, 1,-1,-1, 1,
          1,-1,-1, 1,-1, 1, 1,-1,-1, 1, 1,-1, 1,-1,-1, 1,-1, 1, 1,-1, 1,-1,-1, 1, 1,-1,-1, 1,-1, 1, 1,-1,
         -1, 1, 1,-1, 1,-1,-1, 1, 1,-1,-1, 1,-1, 1, 1,-1, 1,-1,-1, 1,-1, 1, 1,-1,-1, 1, 1,-1, 1,-1,-1, 1],
        [1,-1,-1, 1,-1, 1, 1,-1,-1, 1, 1,-1, 1,-1,-1, 1,-1, 1, 1,-1, 1,-1,-1, 1, 1,-1,-1, 1,-1, 1, 1,-1,
         -1, 1, 1,-1, 1,-1,-1, 1, 1,-1,-1, 1,-1, 1, 1,-1, 1,-1,-1, 1,-1, 1, 1,-1,-1, 1, 1,-1, 1,-1,-1, 1,
         -1, 1, 1,-1, 1,-1,-1, 1, 1,-1,-1, 1,-1, 1, 1,-1, 1,-1,-1, 1,-1, 1, 1,-1,-1, 1, 1,-1, 1,-1,-1, 1,
          1,-1,-1, 1,-1, 1, 1,-1,-1, 1, 1,-1, 1,-1,-1, 1,-1, 1, 1,-1, 1,-1,-1, 1, 1,-1,-1, 1,-1, 1, 1,-1]
    ]
    b = np.array(b)
    pi = np.pi
    q = (int)(cell_id/126)
    u = (cell_id % 126) + 3
    n_prime = 0
    m = 0
    theta = 33 / 132 *frame_number
    frame_shift = -1j * 2 *pi *theta
    root_index = -1j * pi * u / 131
    for n in range(132):
        n_prime = n % 131
        m = n % 128
        nsss_seq[n] = np.exp(frame_shift * n) * np.exp(root_index * n_prime *(n_prime + 1)) *b[q][m]
    return nsss_seq

def channel_estimation(grid,ref_grid_local):
    H_old = np.zeros((8,1),dtype=complex)
    for index in range(len(nrs_x)):
        nrs_local = grid[nrs_x[index],nrs_y[index]]
        ref = ref_grid_local[nrs_x[index],nrs_y[index]]
        H_ele = ref/nrs_local
        H_old[index] = H_ele
    gridx = range(12)
    gridy = range(14)
    H_est = interpolate.interp2d(nrs_x,nrs_y,angle(H_old),kind='linear')(gridy,gridx)
    H_est = np.exp(1j*H_est/(2*pi))
    return H_est

def ofdm_decode(input_seq,ref_grid_local,nsss_flag):
    res = np.zeros((12,14),dtype=complex)
    cp_fraction = 0.55
    for i in range(14):
        if i == 0 or i == 7:
            cp_len = 10
        else:
            cp_len = 9
        fft_start = np.int8(fix(cp_fraction*cp_len))
        fft_start = 0
        phasecorrection = np.array([exp(-1j*2*np.pi*(cp_len-fft_start)/128*idx) for idx in range(128)])
        phasecorrection = 1
        half_sc = np.array([exp(-1j*np.pi*(idx+fft_start-cp_len)/128) for idx in range(128+cp_len)])
        seq = input_seq[i]
        for_fft = seq[0:128+cp_len]*half_sc*phasecorrection
        seq = fft.fftshift(fft.fft(for_fft[cp_len:]))
        seq = seq[58:70]
        res[0:12,i] = seq
    if nsss_flag == 0:
        H_est = channel_estimation(res,ref_grid_local)
        res = res*H_est
    res = res.T
    return res.reshape(12*14,1)

def subframe_remove_cp(subframe_samples):
    symbol = np.zeros((14,138),dtype=complex)
    symbol[0] = subframe_samples[0:138]
    symbol[1][0:-1] = subframe_samples[138:275]
    symbol[2][0:-1] = subframe_samples[275:412]
    symbol[3][0:-1] = subframe_samples[412:549]
    symbol[4][0:-1] = subframe_samples[549:686]
    symbol[5][0:-1] = subframe_samples[686:823]
    symbol[6][0:-1] = subframe_samples[823:960]
    symbol[7] = subframe_samples[960:1098]
    symbol[8][0:-1] = subframe_samples[1098:1235]
    symbol[9][0:-1]= subframe_samples[1235:1372]
    symbol[10][0:-1] = subframe_samples[1372:1509]
    symbol[11][0:-1] = subframe_samples[1509:1646]
    symbol[12][0:-1] = subframe_samples[1646:1783]
    symbol[13][0:-1] = subframe_samples[1783:1920]
    return symbol

def wgn(x, snr):
    snr = 10**(snr/10.0)
    xpower = np.sum(x**2)/len(x)
    npower = xpower / snr
    return np.random.randn(len(x)) * np.sqrt(npower)

def cfo_estimation(NB_frame,ref,dB):
    npss = NB_frame[9600:9600+1920]
    npss = npss.conjugate()
    cfo_seq_raw = angle(npss*ref)
    cfo_seq_raw = cfo_seq_raw[549:1783]
    # plt.plot(np.arange(len(cfo_seq_raw)),cfo_seq_raw,'o',color = 'lightsteelblue',label = 'angle')
    cfo_seq_raw = signal.medfilt(cfo_seq_raw,9)
    for id in range(len(cfo_seq_raw)-8):
        if mean(cfo_seq_raw[id:id+4])-mean(cfo_seq_raw[id+4:id+8]) > 1.4*np.pi:
            cfo_seq_raw[id+4:] = cfo_seq_raw[id+4:] + 2*np.pi
            break
        if mean(cfo_seq_raw[id:id+4])-mean(cfo_seq_raw[id+4:id+8]) < -1.4*np.pi:
            cfo_seq_raw[id+4:] = cfo_seq_raw[id+4:] - 2*np.pi
            break
    cfo_seq = signal.medfilt(cfo_seq_raw,3)
    cfo_len = len(cfo_seq)
    cfo_est,bias = np.polyfit(np.array(range(cfo_len))+1,cfo_seq,1)
    cfo_seq_raw = signal.medfilt(cfo_seq_raw,5)
    # data1 = pd.DataFrame(cfo_seq_raw)
    # data1.to_csv(f'cfo_raw_seq_{dB}.csv')
    # data1 = pd.DataFrame(cfo_seq)
    # data1.to_csv(f'cfo_seq_{dB}.csv')
    # plt.plot(np.arange(len(cfo_seq_raw)),cfo_seq_raw,'o',color = 'royalblue',label = 'angle')
    # plt.plot(cfo_seq,color = 'midnightblue')
    # plt.plot(np.arange(len(cfo_seq)),cfo_est*np.arange(len(cfo_seq))+bias,color = 'midnightblue',label = 'liner regression',linewidth = 5)
    # plt.xlabel('Samples',fontdict={'family':'Times New Roman', 'size':24,'weight':'bold'})
    # plt.ylabel('Phase',fontdict={'family':'Times New Roman', 'size':24,'weight':'bold'})
    # plt.show()
    return cfo_est