# Import libraries
import os
import pandas as pd
import numpy as np
from audiolazy import str2midi 

# Section below takes a CSV file and creates a dictionary of data frames

# Directory path
path = '/mnt/c/Users/kchen/Documents/PD_diffusion_project/CSV_files6/'
#need to replace this with a user-prompt

# Initialize empty dictionary. Each key represents the time
# the spindle is alive, which maps to a df of
# 10 spindles with the highest radii.
spdl = {}
# update the for loop below to find the largest volume and highest y-value of all spindles.
# have a function that takes the highest radii or radii from first spindle and puts it in a list.

# empty variables to store highest values of y-coordinate or radii
max_rad = 0
max_y_c = 0

# For-loop to iterate over csv files
for time, files in enumerate(os.listdir(path)):
    if files.endswith('.csv'):
#         print(files)
        df = pd.read_csv(path + files) 
        df = df.sort_values(by=['x1'], ascending=False) #sort data from highest to lowest radius
        df = df.drop(['x2', 'x3', 'x4'], axis = 1)
        t = time-1 # set time so first frame starts at time = 0.
        spdl[t] = df.head(10)
        
        # extract max radii
        n_rad = spdl[t]['x1'].max()
        if n_rad > max_rad:
            max_rad = n_rad
        # extract max y-coordinate
        n_y_c = spdl[t]['x6'].max()
        if n_y_c > max_y_c:
            max_y_c = n_y_c
        
    else:
        continue
        
# General Mapping Function
def map_value(value, min_value, max_value, min_result, max_result):
    '''maps value (or array of values) from one range to another'''
    
    result = min_result + (value - min_value)/(max_value - min_value)*(max_result - min_result)
    return result

# Replacement Functions
def replace(dic, df_col, m_args): # map_args is a list
    results = []
    for v in (dic.values()):
        for old in v[df_col]:
            new = map_value(old, m_args[0], m_args[1], m_args[2], m_args[3])
            results.append(new)
    return results

def replace2(lists, m_args): # map_args is a list
    results = []
    for old in lists:
        new = map_value(old, m_args[0], m_args[1], m_args[2], m_args[3])
        results.append(new)
    return results

#Maximum and Minimum value functions
def extrema(row, col, which, sort): # which should take a string: 'min' or 'max'
    # ex_list = [] # ex stands for extrema
    ex = 0
    
    for v in (spdl.values()):
        # if sort == "x1" or "x2" or "x3":
        v = v.sort_values(by=[sort], ascending=False)
        
        # print(f"dataframe: {v}")
        
        num = v.iat[row, col]
        # print(f"selected value: {num}")
        
        if which == 'max':
            if num > ex:
                ex = num
        
        elif which == 'min':
            if num < ex:
                ex = num
    
    return ex

#Code below outputs duration fo music for our midi file.

frames_per_beat = 16
#number of frames for each beat of music 

time_per_frame = list(spdl.keys())

t_data = [tpf/frames_per_beat for tpf in time_per_frame]
t_data = list(np.repeat(t_data, 10))

# print(time_per_frame)
# print(t_data)
#print(max(t_data))

#Code below stores a list of all the notes of a grand piano
#with each pitch name being a string.

notes = ['C', 'C#', 'D', 'D#','E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B']

keyboard = []

for i in range(1, 8):
    # octave = []
    for n in notes:
        keyboard.append(str(n + str(i)))
    # keyboard.append(octave)

# insert lowest notes on the keyboard
keyboard = ['A0', 'A#0', 'B0'] + keyboard

# # insert highest note on the keyboard
# keyboard.insert(8, 'C8')
keyboard = keyboard + ['C8']


#Code block below maps data to MIDI note numbers

# first, let's find the minimum and maximum x and y values
x_max = extrema(0, 1, 'max', 'x5')
x_min = extrema(9, 1, 'min', 'x5')
y_max = extrema(0, 2, 'max', 'x6')
y_min = extrema(9, 2, 'min', 'x6')

# second, let's normalize the x and y values

x_norm = replace(spdl, 'x5', [x_min, x_max, 0, 1])

y_norm = replace(spdl, 'x6', [y_min, y_max, 0, 1])
# we also need the max and min of these normed values
x_n_max = max(x_norm)
x_n_min = min(x_norm)
y_n_max = max(y_norm)
y_n_min = min(y_norm)

# third, let's map just the y values to each pitch
y_pitch = replace2(y_norm, [y_n_min, y_n_max, 21, 108]) # refers to pitch class set; I can't use octave twice as a variable name...
# print(y_pitch)
choice = y_pitch
midi_data = [round(p) for p in choice]


#Code block below maps data to MIDI note velocities (dynamics)
# first, normalize radii
r_max = extrema(0, 0, 'max', 'x1')
# r_min = 0
r_min = extrema(9, 0, 'min', 'x1')
# Let's try to set r_min as 0 for now
print(f'max radii: {r_max}')
print(f'min radii: {r_min}')

r_norm = replace(spdl, 'x1', [r_min, r_max, 0, 1])
# print(r_norm)

# second, map normalized radii to velocities
vel_min,vel_max = 35,127   #minimum and maximum note velocity
vel_data = replace2(r_norm, [min(r_norm), max(r_norm), vel_min, vel_max]) #(vel_min, vel_max)
vel_data = [int(i) for i in vel_data]


#Code block below saves data as MIDI File
from midiutil import MIDIFile #import library to make midi file, https://midiutil.readthedocs.io/en/1.2.1/
    
#create midi file object, add tempo
my_midi_file = MIDIFile(1) #one track 
my_midi_file.addTempo(track=0, time=0, tempo=60)
my_midi_file.addProgramChange(0, 0 , 0, 48) # 48 tells midiutil to have a string orchestra play this.
# addProgramChange doesn't work if you use arguments like 'track=0'. So the arguments for this method are: track=0, channel=0, time=0, program=48.
# We want a string orchestra to be playing this, so the program number is 48. Note that it is generally 49, since the program numbers start on 1.
# But for MIDIutil, the program numbers start on 0.

#add midi notes
for i in range(len(midi_data)):
    my_midi_file.addNote(track=0, channel=0, pitch=midi_data[i], time=t_data[i], duration=2, volume=vel_data[i])
# what does duration mean here? The time of each note? And why set it to 2? What if I set it to 1? 
# apparently the music was shortened by a second?

#create and save the midi file itself
with open(path + 'Music/' + 'new_random_music.mid', "wb") as output:
    my_midi_file.writeFile(output) 
