from scipy.io import wavfile  # scipy library to read wav files
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import pdb
import hashlib

def stereoToMono(audiodata):
    newaudiodata = []
    d = audiodata.sum(axis=1) / 2
    newaudiodata.append(d)
    return np.array(newaudiodata, dtype='int16')

def spectrogram(file):
    sample_rate, X = wavfile.read(file)
    B = stereoToMono(X)
    A = B[0].tolist()
    facecolor = '#D0ECE6'
    picname = "C:/Users/User/Desktop/mp3/sp.png"
    a = X.shape[0] / sample_rate
    spectr, frecs, tims, img = plt.specgram(A, Fs=sample_rate, xextent=(0, a))
    band = spectr[15:250]
    sss = np.argmax(band.transpose(), 1)
    first_point = 0
    second_point = 0
    max_arr = []
    for element in range(len(sss)):
        if first_point == 0:
            if len(max_arr) > 0 and \
                    element - max_arr[-1] > 100 or \
                    len(max_arr) == 0:
                first_point = element
        elif second_point == 0:
            second_point = element
        if second_point != 0:
            if second_point - first_point < 100:
                if band[sss[first_point]][first_point] >= \
                        band[sss[second_point]][second_point]:
                    max_arr.append(first_point)
                else:
                    max_arr.append(second_point)
            else:
                max_arr.append(first_point)
            first_point = 0
            second_point = 0
    hashes = []
    for element in max_arr:
        element2 = len(max_arr) - 1
        while max_arr[element2] != element:
            if tims[max_arr[element2]] < tims[element] + 400 and \
                    tims[max_arr[element2]] > tims[element] and \
                    frecs[sss[max_arr[element2]]] < frecs[sss[element]] + 500 and \
                    frecs[sss[max_arr[element2]]] > frecs[sss[element]] - 500:
                one_hash = (frecs[sss[element]],
                            frecs[sss[max_arr[element2]]],
                            round(tims[max_arr[element2]] - tims[element], 3))
                hashes.append([hash(one_hash), round(tims[element], 3)])
            element2 -= 1
    return hashes

def main():
    c, i = 0, 1
    a = spectrogram("C:/Users/User/Desktop/mp3/5.wav")

    while i <= 5:
        max_gist = []
        j = 0
        x = []
        min_s = 0
        f = open('C:/Users/User/Desktop/mp3/' + str(i) + '.txt', 'r')
        for line in f:
            l = ''
            l1 = ''
            flag = 0
            for string in line:
                if string != ' ':
                    l += string
                elif string == ' ' and flag == 0:
                    flag = 1
                    x.append([])
                    x[j].append(int(l))
                    l = ''
                if flag == 1:
                    l1 += string
            x[j].append(float(l1))
            l1 = ''
            j += 1
        print('1')
        j = 0
        first_time_a = 0
        time_a = 0
        for element in a:
            if time_a == 0:
                time_a = element[1]
            else:
                first_time_a = time_a
                time_a = element[1]
            first_time_x = 0
            time_x = 0
            s = min_s
            flag_s = False
            while s < len(x):
                if time_x == 0:
                    time_x = x[s][1]
                else:
                    first_time_x = time_x
                    time_x = x[s][1]
                    flag_c = True
                if element[0] == x[s][0]:
                    if len(max_gist) == 0:
                        min_s = s
                        max_gist.append([])
                        max_gist[j].append(x[s][1])
                        max_gist[j].append(1)
                        j += 1
                    else:
                        i_m = 0
                        flag = False
                        while i_m < len(max_gist):
                            if max_gist[i_m][0] == x[s][1]:
                                max_gist[i_m][1] += 1
                                flag = True
                            i_m += 1
                        if (not flag) and ((time_a - first_time_a) == (time_x - first_time_x)):
                            max_gist[j - 1][1] += 1
                        elif not flag:
                            max_gist.append([])
                            max_gist[j].append(x[s][1])
                            max_gist[j].append(1)
                            flag_c = False
                            j += 1
                    if flag_s == False:
                        min_s = s
                        flag_s = True
                s += 1
        print(max_gist)
        i += 1

main()
