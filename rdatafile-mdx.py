#!/usr/bin/env python

# This code is written for python >= 3.6 and require matplotlib and numpy
#
# Usage: python rdatafile-mdx.py inputfile prefix
#
# inputfile = file name of the CADI data input file (can be either md4 or md2)
# prefix = prefix of the resulting ionogram files
#
# Generate ionograms to file or screen
# The number of ionograms is specified from the variable number_of_ionograms_to_plot below.
# The filetype of the resulting files is spesified in the variable output_file_format below in addition to an integer
# starting with 0 for the first file.
#
# The number of ionograms to plot is used in case the data contain several ionograms. If it only contain one,or fewer
# than the number of ionograms to plot, the program will only plot the ones it can.
#
# The code is written based on IDL code obtained from the Canadian High Arctic Ionospheric Network (CHAIN) web pages;
# http://chain.physics.unb.ca/chain/.
# IDL code originally written by Ian Grant and modified by the same and JWM.
#

import math
from statistics import median
import sys
import struct
import datetime
from time import strptime
import numpy as np
import matplotlib
import copy
# matplotlib.use('Agg')       # use the 'Agg' backend when $DISPLAY environment is not defined
import matplotlib.pyplot as plt

plt.rcParams["figure.figsize"] = (12, 8)

max_ntimes = 256
max_ndopbins = 300000
dheight = 3.0  # not defined in data file

output_file_format = '.png'
number_of_ionograms_to_plot = 5

f = open(sys.argv[1], "rb")
#f = open('testdata.md2', "rb")

f.seek(-1,2)     # go to the file end.
eof = f.tell()   # get the end of file location
f.seek(0,0)      # go back to file beginning

try:
    # 1) read header information as described in the documentation p. 26-27
    site = f.read(3).decode("utf-8")
    ascii_datetime = f.read(22).decode("utf-8")
    filetype = f.read(1).decode("utf-8")
    nfreqs = struct.unpack("<H", f.read(2))[0]
    ndops = struct.unpack("<B", f.read(1))[0]
    minheight = struct.unpack("<H", f.read(2))[0]
    maxheight = struct.unpack("<H", f.read(2))[0]
    pps = struct.unpack("<B", f.read(1))[0]
    npulses_avgd = struct.unpack("<B", f.read(1))[0]
    base_thr100 = struct.unpack("<H", f.read(2))[0]
    noise_thr100 = struct.unpack("<H", f.read(2))[0]
    min_dop_forsave = struct.unpack("<B", f.read(1))[0]
    dtime = struct.unpack("<H", f.read(2))[0]
    gain_control = f.read(1).decode("utf-8")
    sig_process = f.read(1).decode("utf-8")
    noofreceivers = struct.unpack("<B", f.read(1))[0]
    spares = f.read(11).decode("utf-8")

    month = ascii_datetime[1:4]
    day = int(ascii_datetime[5:7])
    hour = int(ascii_datetime[8:10])
    minute = int(ascii_datetime[11:13])
    sec = int(ascii_datetime[14:16])
    year = int(ascii_datetime[17:21])

    month_number = strptime(month, '%b').tm_mon
    mydate = datetime.date(year, month_number, day)
    jd = mydate.toordinal() + 1721424.5
    jd0jd = datetime.date(1986, 1, 1)
    jd0 = jd0jd.toordinal() + 1721424.5

    time_header = (jd - jd0) * 86400 + hour * 3600 + minute * 60 + sec
    time_hour = 3600 * (time_header / 3600)

    # 2) read all frequencies used

    freqs = [struct.unpack("<f", f.read(4))[0] for i in range(nfreqs)]

    if filetype == 'I':
        max_nfrebins = nfreqs
    else:
        max_nfrebins = min(max_ntimes * nfreqs, max_ndopbins)

    # 3) read rawdata

    nheights = int(maxheight / dheight + 1)

    times = []
    frebins = []
    frebins_x = []
    frebins_gain_flag = []
    frebins_noise_flag = []
    frebins_noise_power10 = []
    time_min = 0
    time_sec = 0
    timex = -1
    freqx = nfreqs - 1
    dopbinx = -1
    frebinx = -1
    iq_bytes = np.zeros((noofreceivers, 2))
    dopbin_x_timex = []
    dopbin_x_freqx = []
    dopbin_x_hflag = []
    dopbin_x_dop_flag = []
    dopbin_iq = []
    hflag = 0

    time_min = struct.unpack("<B", f.read(1))[0]
    while time_min != 255:
        time_sec = struct.unpack("<B", f.read(1))[0]
        flag = struct.unpack("<B", f.read(1))[0]  # gainflag
        timex += 1
        times.append(time_hour + 60 * time_min + time_sec)
        for freqx in range(nfreqs):
            noise_flag = struct.unpack("<B", f.read(1))[0]  # noiseflag
            noise_power10 = struct.unpack("<H", f.read(2))[0]
            frebinx += 1
            frebins_gain_flag.append(flag)
            frebins_noise_flag.append(noise_flag)
            frebins_noise_power10.append(noise_power10)
            flag = struct.unpack("<B", f.read(1))[0]
            while flag < 224:
                ndops_oneh = struct.unpack("<B", f.read(1))[0]
                hflag = flag
                if ndops_oneh >= 128:
                    ndops_oneh = ndops_oneh - 128
                    hflag = hflag + 200
                for dopx in range(ndops_oneh):
                    dop_flag = struct.unpack("<B", f.read(1))[0]
                    for rec in range(noofreceivers):
                        iq_bytes[rec, 0] = struct.unpack("<B", f.read(1))[0]
                        iq_bytes[rec, 1] = struct.unpack("<B", f.read(1))[0]
                    dopbinx += 1
                    dopbin_iq.append(copy.deepcopy(iq_bytes))
                    dopbin_x_timex.append(timex)
                    dopbin_x_freqx.append(freqx)
                    dopbin_x_hflag.append(hflag)
                    dopbin_x_dop_flag.append(dop_flag)
                flag = struct.unpack("<B", f.read(1))[0]  # next hflag/gainflag/FF
        time_min = flag
        if ((f.tell() - 1) != eof):
            time_min = struct.unpack("<B", f.read(1))[0]  # next record
finally:
    f.close()

# frebins_gain_flag[] is the gainflag for each frequency
# frebins_noise_flag[] is the noiseflag for each frequency
# frebins_noise_power10[] is the 10 x averagenoisepower for each frequency

# dopbin_iq[] is an array of real and imaginary components for the receivers

num_of_ionograms = min(timex, number_of_ionograms_to_plot - 1) + 1

for j in range(num_of_ionograms):
    # calculate average power
    block = timex - j
    meanpower = []
    for idx in range(len(dopbin_iq)):
        if dopbin_x_timex[idx] == int(block):
            absvalue1 = math.sqrt((dopbin_iq[idx][0][0] - 128) ** 2 + (dopbin_iq[idx][0][1] - 128) ** 2)
            absvalue2 = math.sqrt((dopbin_iq[idx][1][0] - 128) ** 2 + (dopbin_iq[idx][1][1] - 128) ** 2)
            absvalue3 = math.sqrt((dopbin_iq[idx][2][0] - 128) ** 2 + (dopbin_iq[idx][2][1] - 128) ** 2)
            absvalue4 = math.sqrt((dopbin_iq[idx][3][0] - 128) ** 2 + (dopbin_iq[idx][3][1] - 128) ** 2)
            mvalue = median([absvalue1, absvalue2, absvalue3,absvalue4])  # median
            # mvalue = (absvalue1 + absvalue2 + absvalue3 + absvalue4) / 4. # mean
            if mvalue == 0:
                power = 0
            else:
                power = 20 * math.log10(mvalue)
            meanpower.append(power)

    # frequency in MHz
    frequency = []
    height = []
    for i in range(len(dopbin_x_freqx)):
        if dopbin_x_timex[i] == int(block):
            frequency.append(freqs[dopbin_x_freqx[i]] / 1000000.0)
            height.append(dopbin_x_hflag[i] * 3)

    # plotting

    title = 'CADI ionogram ' + site + ' - ' + ascii_datetime[1:11] + '{0}'.format(str(block).zfill(2)) + ':' + '{0}'.format(str(time_sec).zfill(2)) + ascii_datetime[16:21] + ' UTC'

    fig, ax = plt.subplots()

    cm = plt.cm.get_cmap('jet')
    sc = plt.scatter(frequency, height, s=1, c=meanpower, cmap=cm)
    ax.grid(True, which='both')

    ax.set_xlim(0, math.floor(max(frequency) + 1))  # set X limits (min and max frequency in MHz)
    ax.set_ylim(0, math.floor((maxheight + 100)/100) * 100)  # set Y limits (min and max height in km)
    ax.set_title(title)
    ax.set_xlabel('Frequency (Mhz)')
    ax.set_ylabel('Virtual height (km)')
    # plt.clim(43,45)
    cbar = plt.colorbar(sc)
    cbar.set_label('Power (dB)')

    # toggle output to file, comment/uncomment the next two lines
    # if you don't want output to file: usage: python rdatafile.py inputfile

    filename = sys.argv[2] + str(j) + output_file_format
    fig.savefig(filename)

# toggle output to screen, comment/uncomment the next line (note you can have both!)
plt.show()