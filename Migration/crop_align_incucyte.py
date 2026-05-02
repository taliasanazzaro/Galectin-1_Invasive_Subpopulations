import pandas as pd
import matplotlib.pyplot as plt
import cv2 as cv
import os
import numpy as np
from multiprocessing import Manager, Process
import pickle

cell_line = 'GS54'

confluence_path = 'GS54_confluence_2'
phase_path = 'GS54_phase'

save_path = 'pickles'

confluence_save_file = f'{save_path}/{cell_line}_confluence.pkl'
phase_save_file = f'{save_path}/{cell_line}_phase.pkl'
adjust_save_file = f'{save_path}/{cell_line}_adjust.pkl'
times_save_file = f'{save_path}/{cell_line}_times.pkl'

median_cutoff = 0.1
correction_max_shift = 150

def import_sequence(directory, prefix=''):
    imgs = {}
    count = 0

    print(f'loading {directory}')

    for file in os.listdir(directory):
        f = os.path.join(directory, file)

        if os.path.isdir(f) or not f.endswith('.tif') or not file.startswith(prefix):
            continue

        well = file.split('_')[1]+'_'+file.split('_')[2]
        time = file.split('_')[3][:-7]
        
        if well not in imgs.keys():
            imgs[well] = {} 

        imgs[well][time] = cv.imread(f, 0)
        count += 1
    
    print(f'\t{count} images loaded')

    return imgs

def sequence_to_3d_matrix(sequence, times, warn=False):
    not_found_times = []
    missed_times = []

    mat = np.zeros((len(times), list(sequence.values())[0].shape[0], list(sequence.values())[0].shape[1]))

    for i, time in enumerate(times):
        if time in sequence.keys():
            mat[i] = sequence[time]
        else:
            not_found_times.append(time)

    for seq_time in sequence.keys():
        if seq_time not in times:
            missed_times.append(times)

    if warn:
        if len(not_found_times) > 0:
            print(f'{len(not_found_times)} defined times not found in sequence: {not_found_times}')
        if len(missed_times) > 0:
            print(f'{len(missed_times)} times in sequence not defined: {missed_times}')
    
    return mat

confluence = import_sequence(confluence_path)
phase = import_sequence(phase_path)

times = [np.sort(np.array([time for time in confluence[well].keys()])) for well in confluence.keys()]
assert np.all(times[i] == times[i-1] for i in range(len(times)-1)), 'Times do not match across wells'
times = times[0]

focus = pd.read_table(f'Focusposition_{cell_line}.txt', header=1, index_col=1).drop(columns='Date Time')
focus.index = times
focus.columns = [s.split(',')[0]+'_'+s[-1] for s in focus.columns]

fig, axs = plt.subplots(4, 12, figsize=(12*2, 4*2))

median = {}
skip_times = {}

for i, well in enumerate(focus.columns):
    ax = axs[i%4, i//4]
    
    median[well] = np.median(focus[well])
    skip_times[well] = focus.index[((median[well] + median[well]*median_cutoff) < focus[well]) | ((median[well] - median[well]*median_cutoff) > focus[well])]

    ax.plot(focus.index, focus[well], color='black')
    ax.scatter(skip_times[well], focus.loc[skip_times[well], well], color='red')
    ax.set_title(well)
    ax.set_ylim([3000, 5000])
    ax.set_xticks([])
    if i//5 > 0:
        ax.set_yticks([])
well_times = {well:np.array([t for t in times if t not in skip_times[well]], dtype=str) for well in confluence.keys()}
confluence = {well:sequence_to_3d_matrix(confluence[well], well_times[well]) for well in confluence.keys()}
phase = {well:sequence_to_3d_matrix(phase[well], well_times[well]) for well in phase.keys()}

for well in confluence.keys():
    confluence[well] = confluence[well] != 0

def calc_overlap(well, max_shift, overlap):
    m = np.zeros((confluence[well].shape[0], max_shift*2+1, max_shift*2+1))

    compare_to = confluence[well][:-1, correction_max_shift:confluence[well].shape[1]-correction_max_shift, correction_max_shift:confluence[well].shape[2]-correction_max_shift]

    for i, x in enumerate(range(-max_shift, max_shift+1)):
        print(f'\t{well}: {x}')
        for j, y in enumerate(range(-max_shift, max_shift+1)):
            m[1:, i, j] = np.sum(compare_to == confluence[well][1:, i:confluence[well].shape[1]-2*correction_max_shift+i, j:confluence[well].shape[2]-2*correction_max_shift+j], axis=(1,2))

    overlap[well] = m

    return overlap

overlap = {}

print('Calculating Overlaps')

with Manager() as manager:
    overlap_m = manager.dict()
    processes = []

    # creating processes
    for well in confluence.keys():
        p = Process(target=calc_overlap, args=(well, correction_max_shift, overlap_m)) # was using ~outer_area
        processes.append(p)
        p.start()

    # completing process
    for p in processes:
        p.join()

    for well in overlap_m.keys():
        overlap[well] = overlap_m[well]
offset_x = {}
offset_y = {}

print('Calculating Offset Coordinates')
for well in overlap.keys():
    print('\t'+well)
    offset_x[well] = np.zeros(confluence[well].shape[0], dtype=int)
    offset_y[well] = np.zeros(confluence[well].shape[0], dtype=int)

    for t in range(1, len(well_times[well])):
        x = np.where(overlap[well][t] == np.max(overlap[well][t]))[0]-correction_max_shift
        y = np.where(overlap[well][t] == np.max(overlap[well][t]))[1]-correction_max_shift

        xy = x + y
        offset_x[well][t] = x[np.argmin(xy)]+offset_x[well][t-1]
        offset_y[well][t] = y[np.argmin(xy)]+offset_y[well][t-1]

    # adjust the values so they are centered around the average movement
    offset_x[well] = offset_x[well]-int(np.mean(offset_x[well], dtype=int))
    offset_y[well] = offset_y[well]-int(np.mean(offset_y[well], dtype=int))

phase_crop = {}
confluence_crop = {}

print('Cropping')
for well in phase.keys():
    print(f'\t{well}')
    dim_x = confluence[well].shape[1]-offset_x[well].max()+offset_x[well].min()
    dim_y = confluence[well].shape[2]-offset_y[well].max()+offset_y[well].min()

    if dim_x > 0 and dim_y > 0:
        confluence_crop[well] = np.zeros((len(well_times[well]), dim_x, dim_y), dtype=bool)
        phase_crop[well] = np.zeros((len(well_times[well]), dim_x, dim_y), dtype=np.uint8)

        for t in range(len(well_times[well])):
            confluence_crop[well][t] = confluence[well][t, offset_x[well][t]-offset_x[well].min():dim_x+offset_x[well][t]-offset_x[well].min(), offset_y[well][t]-offset_y[well].min():dim_y+offset_y[well][t]-offset_y[well].min()]
            phase_crop[well][t] = phase[well][t, offset_x[well][t]-offset_x[well].min():dim_x+offset_x[well][t]-offset_x[well].min(), offset_y[well][t]-offset_y[well].min():dim_y+offset_y[well][t]-offset_y[well].min()]
    else:
        print(f'Negative dimensions for {well}: ({dim_x}, {dim_y})')


print('Writing Files')

with open(phase_save_file, 'wb') as handle:
    pickle.dump(phase_crop, handle)
with open(confluence_save_file, 'wb') as handle:
    pickle.dump(confluence_crop, handle)
with open(adjust_save_file, 'wb') as handle:
    pickle.dump((offset_x, offset_y), handle)
with open(times_save_file, 'wb') as handle:
    pickle.dump(well_times, handle)