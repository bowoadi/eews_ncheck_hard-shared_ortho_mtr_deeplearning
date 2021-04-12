#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import math
import os
import h5py
import datetime
from sklearn.model_selection  import train_test_split,StratifiedKFold,cross_val_score
pd.set_option('display.max_columns', None)
from obspy.core import UTCDateTime
from obspy import read,read_inventory
import warnings
warnings.filterwarnings("ignore")


from scipy.signal import butter, lfilter, filtfilt
def butter_bandpass(lowcut, highcut, fs, order):
  nyq = 0.5 * fs
  low = lowcut / nyq
  high = highcut / nyq
  b, a = butter(order, [low, high], btype='band')
  return b, a
 
#band pass filter with filtfilt
def butter_bandpass_filter(data, lowcut, highcut, fs, order):
  b, a = butter_bandpass(lowcut, highcut, fs, order=order)
  # y = lfilter(b, a, data)
  y =filtfilt(b, a, data)
  return y

file_dir = "/E/Gempa/EarthquakeCataloge/Jatim"
pick_dir = "/E/Gempa/EarthquakeCataloge/Jatim/Picking Jatim"

pick_list = []
for file in os.listdir(pick_dir):
    if file.endswith(".pick"):
        pick_list+=[file[:15]]
print(f'total events : {len(pick_list)}')


file_in_folder_path = "/E/Gempa/Dataset-Master/plot-folder"

p3_list = [file.split(".")[0] for file in os.listdir(file_in_folder_path)]
print(f'total events : {len(p3_list)}')

ev_list = []
for file in os.listdir(file_dir):
    if file.endswith(".mseed"):
        if (file.split(".")[0] in p3_list):
            ev_list+=[file.split(".")[0]]
print(f'total events : {len(ev_list)}')

def generate_dataset(filename):
    st = read(f"{file_dir}/{filename}.mseed")
    column = ["id","network","station","location","channel","starttime","endtime","sampling_rate",
              "delta","npts","calib","_format","dataquality","number_of_records","encoding","byteorder",
             "record_length","filesize","data"]
    df = pd.DataFrame(columns=column)
    for i in range(len(st)):
        tr = st[i]
        data = []
        for j in range(1,12):
            data+=[tr.stats[column[j]]]
        for j in range(12,18):
            data+=[tr.stats.mseed[column[j]]]
        data += [tr.data]
        #with id
        data = [filename + "_" +tr.stats.station] + data
        df.loc[i] = data
        
    #file pick to csv
    column = ["station", "c1" , "channel","c2","wave","c3","time","c4","c5","c6","c7","c8"]
    df_pick = pd.DataFrame(columns=column)
    f = open(f"{pick_dir}/{filename}.pick", "r")
    pick_list = f.read().split("\n")
    for i,pick in enumerate(pick_list):
        if(pick !=""):
            pick = pick.split(" ")
            data = pick[0:6]
            data += [UTCDateTime(pick[6][0:4]+"-"+pick[6][4:6]+'-'+pick[6][6:8]+"T"+pick[7][0:2]+":"+pick[7][2:4]+":"+pick[8] + "Z")]
            data +=pick[9:14]
            df_pick.loc[i] = data
    #get earthquake data
    list_station = df_pick.station.unique()
    list_station = sorted(list_station) # sort by name
    #check pick, if empty, copy pick from coupled event
    for stasiun_check in ["GMJI","JAGI","PWJI"]: 
        if stasiun_check not in list_station: # check station
            data_pick = df_pick[df_pick.station==[ station_pick for station_pick in list_station if station_pick in ["GMJI","JAGI","PWJI"]][0] ]
            for s_iter in range(len(data_pick)):
                data_pick.station = stasiun_check
            df_pick = pd.concat([df_pick,data_pick]).reset_index(drop=True)
    list_station = df_pick.station.unique()
    df_gempa = df.copy()
    df_gempa["gempa"] = df.station.apply(lambda p: p in list_station)
    df_gempa = df_gempa[df_gempa.gempa == True].sort_values(["station","channel"])
    #create dataset       
    column = ["starttime","endtime","npts",
             "E_waveform","N_waveform","Z_waveform","p_arrival_time","s_arrival_time","p_arrival_sample","s_arrival_sample"]
    dataset = pd.DataFrame(columns=column)
    df_gempa = df_gempa[(df_gempa.station == "GMJI") | (df_gempa.station == "PWJI") | (df_gempa.station == "JAGI")]
    for i in df_gempa.index[::3]:
        #get all
        starttime = max(df.starttime[i],df.starttime[i+1],df.starttime[i+2])
        endtime = min(df.endtime[i],df.endtime[i+1],df.endtime[i+2])
        if(df.channel[i] == "BHE"):
          E_waveform = df.data[i]
        else:
          continue
        if(df.channel[i+1] == "BHN"):
          N_waveform = df.data[i+1]
        else:
          continue
        if(df.channel[i+2] == "BHZ"):
          Z_waveform = df.data[i+2]
        else:
          continue
        #set starttime dan endtime
        E_waveform = E_waveform[int(round((starttime - df.starttime[i])*df.sampling_rate[i])):int(len(E_waveform) - (round((df.endtime[i] - endtime)*df.sampling_rate[i])))]
        N_waveform = N_waveform[int(round((starttime - df.starttime[i+1])*df.sampling_rate[i+1])):int(len(N_waveform) - (round((df.endtime[i+1] - endtime)*df.sampling_rate[i+1])))]
        Z_waveform = Z_waveform[int(round((starttime - df.starttime[i+2])*df.sampling_rate[i+2])):int(len(Z_waveform) - (round((df.endtime[i+2] - endtime)*df.sampling_rate[i+2])))]
        
        #bandpass filter
        fs = df.sampling_rate[i]
        lowcut = 1.0
        highcut = 5.0
        order = 5
        E_waveform = butter_bandpass_filter(E_waveform, lowcut, highcut, fs,order)                                  
        N_waveform = butter_bandpass_filter(N_waveform, lowcut, highcut, fs,order)
        Z_waveform = butter_bandpass_filter(Z_waveform, lowcut, highcut, fs,order)
        
        #normalization
        E_waveform = E_waveform/np.amax(np.abs(E_waveform))
        N_waveform = N_waveform/np.amax(np.abs(N_waveform))
        Z_waveform = Z_waveform/np.amax(np.abs(Z_waveform))
        
        #get p_arrival dan s_arrival
        s_arrival_time = UTCDateTime(str(df_pick[(df_pick.station == df_gempa.station[i]) & (df_pick.wave == "S")].time)[5:32])
        p_arrival_time = UTCDateTime(str(df_pick[(df_pick.station == df_gempa.station[i]) & (df_pick.wave == "P")].time)[5:32])
        s_arrival_sample = int(round((s_arrival_time - starttime)*df.sampling_rate[i]))
        p_arrival_sample = int(round((p_arrival_time - starttime)*df.sampling_rate[i]))
        if((len(E_waveform) == len(N_waveform)) & (len(E_waveform) == len(Z_waveform)) &(len(E_waveform)!=0)):
          dataset.loc[i] = [starttime,endtime,len(E_waveform),E_waveform,N_waveform,Z_waveform,p_arrival_time,s_arrival_time,
                          p_arrival_sample,s_arrival_sample]
        else:
          # print(df_gempa.station[i] , len(E_waveform),len(N_waveform),len(Z_waveform))
          length = min(len(E_waveform),len(N_waveform),len(Z_waveform))
          E_waveform = E_waveform[:length]
          N_waveform = N_waveform[:length]
          Z_waveform = Z_waveform[:length]
          dataset.loc[i] = [starttime,endtime,len(E_waveform),E_waveform,N_waveform,Z_waveform,p_arrival_time,s_arrival_time,
                          p_arrival_sample,s_arrival_sample]
    dataset[['id', 'network', 'station', 'location','sampling_rate', 'delta',
                             'calib']] = df_gempa[['id', 'network', 'station', 'location','sampling_rate', 'delta',
                             'calib']].loc[::3]
    dataset = dataset.set_index("id")
    return dataset


# CSV Signal + Plotting

plot_dir = "./plot-waveform"

window = 40 #detik
start_plot= 20 #seconds before p arrival
end_plot = 20  #seconds after p arrival 
geser =  30# w
number_plot= ((start_plot + end_plot)/ geser ) + 1
start_file = 0
length = 1

for j,filename in enumerate(ev_list[start_file:length+start_file]):
    try:
      dataset_ = generate_dataset('20170428_014929')
      if(dataset_.shape[0] == 3):
        for time in range(int(number_plot)):
            dataset = dataset_.copy()
            starttime = dataset.p_arrival_time.min() - datetime.timedelta(seconds =start_plot - (geser *time))  
            endtime = starttime + datetime.timedelta(seconds = window)
            for i in range(3):
                dataset.E_waveform[i] = dataset.E_waveform[i][int(round((starttime - dataset.starttime[i])*dataset.sampling_rate[i])):int(len(dataset.E_waveform[i]) - (round((dataset.endtime[i] - endtime)*dataset.sampling_rate[i])))]
                dataset.N_waveform[i] = dataset.N_waveform[i][int(round((starttime - dataset.starttime[i])*dataset.sampling_rate[i])):int(len(dataset.N_waveform[i]) - (round((dataset.endtime[i] - endtime)*dataset.sampling_rate[i])))]
                dataset.Z_waveform[i] = dataset.Z_waveform[i][int(round((starttime - dataset.starttime[i])*dataset.sampling_rate[i])):int(len(dataset.Z_waveform[i]) - (round((dataset.endtime[i] - endtime)*dataset.sampling_rate[i])))]
                dataset.npts[i] = len(dataset.E_waveform[i])
            #         dataset.s_arrival_sample[i] = int(round((dataset.s_arrival_time[i] - starttime)*dataset.sampling_rate[i]))
            #         dataset.p_arrival_sample[i] = int(round((dataset.p_arrival_time[i] - starttime)*dataset.sampling_rate[i]))

            fig, axs = plt.subplots(3, 1, figsize=(6, 8));
            fig.subplots_adjust(hspace=0);
            axs[0].plot(dataset.E_waveform[0],color = "red",alpha = 0.3);
            axs[0].plot(dataset.N_waveform[0],color = "green",alpha = 0.3);
            axs[0].plot(dataset.Z_waveform[0],color = "blue",alpha = 0.3);
            axs[0].set_xlabel('');
            axs[0].set_facecolor((0, 0, 0));
            axs[0].set_ylim(ymin = -1, ymax =1);
            axs[0].set_ylabel('')

            axs[1].plot(dataset.E_waveform[1],color = "red",alpha = 0.3);
            axs[1].plot(dataset.N_waveform[1],color = "green",alpha = 0.3);
            axs[1].plot(dataset.Z_waveform[1],color = "blue",alpha = 0.3);
            axs[1].set_xlabel('');
            axs[1].set_facecolor((0, 0, 0));
            axs[1].set_ylim(ymin = -1, ymax =1);
            axs[1].set_ylabel('')

            axs[2].plot(dataset.E_waveform[2],color = "red",alpha = 0.3);
            axs[2].plot(dataset.N_waveform[2],color = "green",alpha = 0.3);
            axs[2].plot(dataset.Z_waveform[2],color = "blue",alpha = 0.3);
            axs[2].set_xlabel('');
            axs[2].set_facecolor((0, 0, 0));
            axs[2].set_ylim(ymin = -1, ymax =1);
            axs[2].set_ylabel('')
            if(not(filename in os.listdir(plot_dir))):
                os.mkdir(f"{plot_dir}\\{filename}")
            plt.savefig(f"{plot_dir}\\{filename}\\{str(time).zfill(6)}_{filename}.png");
            plt.close()
        print(j,filename)
      else:
        print("Skip",j,filename)
    except:
      print("Error",j,filename)



