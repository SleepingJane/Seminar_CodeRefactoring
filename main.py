from scipy.io import wavfile  # scipy library to read wav files
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

def stereoToMono(audiodata):
    newaudiodata = []
    d = audiodata.sum(axis=1) / 2
    newaudiodata.append(d)
    return np.array(newaudiodata, dtype='int16')

def findLocalMaxSpectr(triple_spectrogram, bandwidth):
    first_point = 0
    second_point = 0
    spectr_local_max = []

    for element in range(len(triple_spectrogram)):
        if first_point == 0:
            if len(spectr_local_max) > 0 and \
                    element - spectr_local_max[-1] > 100 or \
                    len(spectr_local_max) == 0:
                first_point = element
        elif second_point == 0:
            second_point = element
        if second_point != 0:
            if second_point - first_point < 100:
                if bandwidth[triple_spectrogram[first_point]][first_point] >= \
                        bandwidth[triple_spectrogram[second_point]][second_point]:
                    spectr_local_max.append(first_point)
                else:
                    spectr_local_max.append(second_point)
            else:
                spectr_local_max.append(first_point)
            first_point = 0
            second_point = 0
    return spectr_local_max

def createHash(triple_spectrogram, spectr_local_max, frecs, tims):
    hashes = []
    time_inaccuracy = 400
    frec_inaccuracy = 500

    for element in spectr_local_max:
        element2 = len(spectr_local_max) - 1
        while spectr_local_max[element2] != element:
            # Условие проверки единственности локального максимума в заданных окрестностях
            if tims[spectr_local_max[element2]] < tims[element] + time_inaccuracy and \
                    tims[spectr_local_max[element2]] > tims[element] and \
                    frecs[triple_spectrogram[spectr_local_max[element2]]] < frecs[
                triple_spectrogram[element]] + frec_inaccuracy and \
                    frecs[triple_spectrogram[spectr_local_max[element2]]] > frecs[
                triple_spectrogram[element]] - frec_inaccuracy:
                one_hash = (frecs[triple_spectrogram[element]],
                            frecs[triple_spectrogram[spectr_local_max[element2]]],
                            round(tims[spectr_local_max[element2]] - tims[element], 3))
                hashes.append([hash(one_hash), round(tims[element], 3)])
            element2 -= 1

    return hashes


def spectrogram(file):
    sample_rate, wav_file_data = wavfile.read(file)
    mono_audio = stereoToMono(wav_file_data)
    audio_data = mono_audio[0].tolist()

    xextent_max = wav_file_data.shape[0] / sample_rate
    spectr, frecs, tims, img = plt.specgram(audio_data, Fs=sample_rate, xextent=(0, xextent_max))
    bandwidth = spectr[15:250]
    triple_spectrogram = np.argmax(bandwidth.transpose(), 1)

    spectr_local_max = findLocalMaxSpectr(triple_spectrogram, bandwidth)
    hashes = createHash(triple_spectrogram, spectr_local_max, frecs, tims)

    return hashes


def readHash(hashes_file):
    hash_array = []
    line_counter = 0

    for line in hashes_file:

        hash_line = ''
        full_hash_line = ''
        flag_lineend = False

        for string in line:
            if string != ' ':
                hash_line += string
            elif string == ' ' and not flag_lineend:
                flag_lineend = True
                hash_array.append([])
                hash_array[line_counter].append(int(hash_line))
                hash_line = ''
            if flag_lineend:
                full_hash_line += string
        hash_array[line_counter].append(float(full_hash_line))
        full_hash_line = ''
        line_counter += 1
    return hash_array


def findMinimumConstellation(hash_array, element, spectr_time):
    max_gist = []
    min_constellation, constellation, gist_element = 0, 0, 0
    first_time_hash, time_hash = 0, 0
    flag_constellation = False

    while constellation < len(hash_array):

        if time_hash != 0:
            first_time_hash = time_hash

        time_hash = hash_array[constellation][1]

        if element[0] == hash_array[constellation][0]:
            if len(max_gist) == 0:
                flag_constellation = False
                min_constellation = constellation
                max_gist.append([])

            else:
                if (not flag_constellation) and (spectr_time == time_hash):
                    max_gist[gist_element - 1][1] += 1
                elif not flag_constellation:
                    max_gist.append([])

            max_gist[gist_element].append(hash_array[constellation][1])
            max_gist[gist_element].append(1)
            gist_element += 1

            min_constellation = constellation
            flag_constellation = True

        constellation += 1
    return max_gist

def main():
    c, i = 0, 1
    file_name = "C:/Users/User/Desktop/mp3/5.wav"
    folder_name = "C:/Users/User/Desktop/mp3/"
    audio_spectr = spectrogram(file_name)

    audio_list_size = 5

    while i <= audio_list_size:
        max_gist = []
        min_s, j = 0, 0
        hashes_file = open(folder_name + str(i) + '.txt', 'r')
        hash_array = readHash(hashes_file)

        first_spectr_time = 0
        spectr_time = 0

        for element in audio_spectr:
            if spectr_time != 0:
                first_spectr_time = spectr_time
            spectr_time = element[1]
            max_gist = findMinimumConstellation(hash_array, element, spectr_time)
        i += 1

main()
