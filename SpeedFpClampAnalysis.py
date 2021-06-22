"""
SpeedFpClamp Data Analysis 
Ricky Pimentel
UNC Chapel Hill Applied Biomechanics Laboratory
2021

Run the script to perform data analysis and generate all article figures
Data avaible at https://drive.google.com/file/d/1PrpgwxUbaDNYojghtbIORW3qLK66NI31/view?usp=sharing
"""

import pandas as pd
import numpy as np
import os
import fnmatch
import matplotlib as mpl
import matplotlib.pyplot as plt
import datetime as dt
import scipy.io as sio
import scipy.stats as stats
import pingouin as pg

# identify folder with data files - RENAME to your folder path!
folder = r'E:\UNC_ABL\FpMetabolics_2020\MetData'

files = os.listdir(folder)
os.chdir(folder)
sty = 'seaborn'
mpl.style.use(sty)

#%% RMR Analysis
# initialize dict that will hold metabolic data
Subjects = {
    's001': {'RMR Avg':'NaN'},
    's002': {'RMR Avg':'NaN'},
    's003': {'RMR Avg':'NaN'},
    's004': {'RMR Avg':'NaN'},
    's005': {'RMR Avg':'NaN'},
    's006': {'RMR Avg':'NaN'},
    's007': {'RMR Avg':'NaN'},
    's008': {'RMR Avg':'NaN'},
    's009': {'RMR Avg':'NaN'},
    's010': {'RMR Avg':'NaN'},
    's011': {'RMR Avg':'NaN'},
    's012': {'RMR Avg':'NaN'},
    's013': {'RMR Avg':'NaN'},
    's014': {'RMR Avg':'NaN'},
    's015': {'RMR Avg':'NaN'},
    's016': {'RMR Avg':'NaN'},
    's017': {'RMR Avg':'NaN'},
    's018': {'RMR Avg':'NaN'},
    's019': {'RMR Avg':'NaN'},
    's020': {'RMR Avg':'NaN'},
    }

SubjNames = list(Subjects.keys())
print(SubjNames)

TrialAvg = {}
v = plt.get_cmap('jet')
cNorm  = mpl.colors.Normalize(vmin=0, vmax=len(Subjects))
scalarMap = mpl.cm.ScalarMappable(norm=cNorm, cmap='jet')

# get RMR file 
pattern = '*REE*'
matching = fnmatch.filter(files, pattern)
d = 0
for i in matching:
    for s in SubjNames:
        if s in i:
            Subj = s
            break 
        
    F = folder + '\\' + i
    RMR = pd.read_excel(F)
    
    # pull VO2 and time data
    VO2_kg = RMR.loc[2:len(RMR),'VO2/Kg']
    VO2 = RMR.loc[2:len(RMR),'VO2']
    VCO2 = RMR.loc[2:len(RMR),'VCO2']
    t = RMR.loc[2:len(RMR),'t'].values
    W = (VO2 /1000 * 16.5 + VCO2 /1000 * 4.51) * 1000 / 60 
    
    # find rows after 3 min
    T = []
    c = 0
    for i in t:
        c = c + 1
        if i.minute >=3:
             T.append(c)  
             
    # calculate average RMR and make array
    AvgVO2_kg = np.mean(VO2_kg[T])
    AvgVO2 = np.mean(VO2[T])
    AvgVCO2 = np.mean(VCO2[T])
    AvgW = np.mean(W[T])
    Seq = list(range(0, len(T)))
    A = np.ones((len(T)))
    AvgRMR_array = A * AvgW
    colorVal = scalarMap.to_rgba(d)
    plt.plot(W[T], color=colorVal, lw=2)
    plt.plot(T, AvgRMR_array, color=colorVal, lw=2)
    plt.show()
            
    Subjects[Subj]['RMR VO2_kg Avg'] = AvgVO2_kg
    Subjects[Subj]['RMR VO2 Avg'] = AvgVO2
    Subjects[Subj]['RMR VCO2 Avg'] = AvgVCO2
    Subjects[Subj]['RMR W Avg'] = AvgW

    print(s + ' Average RMR = ')
    print(AvgVO2_kg)
    print('mL/kg/min')
    d = d+1
    
plt.title('Resting Metabolic Rate')
plt.ylabel('W')
plt.xlabel('Record #')
plt.savefig('RMR.jpg', dpi=300)

#%% Load and extract Active VO2 data
AvgS = []
AvgF = []
plt.close('all')
# get active trial file 
pattern = '*CPET*'
matching = fnmatch.filter(files, pattern)
fig = plt.figure(figsize=[12,10])
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8]) # main axes
MkrSz = 10

d = 0
for i in matching:
    colorVal = scalarMap.to_rgba(d)
    
    for s in SubjNames:
        if s in i:
            Subj = s
            break 
        
    if '$' in i:
        break
        
    print(s)
    Folder = folder + '\\' + i
    Active = pd.read_excel(Folder)

    # pull VO2 and time data
    VO2_kg = Active.loc[2:len(Active),'VO2/Kg'].values.tolist()
    VO2 = Active.loc[2:len(Active),'VO2'].values.tolist()
    VCO2 = Active.loc[2:len(Active),'VO2'].values.tolist()
    t = Active.loc[2:len(Active),'t'].values.tolist()
    Mkr = Active.loc[2:len(Active),'Marker'].values.tolist()
    
    
    # convert start and end times to datetimes
    today = dt.datetime.today()
    ProtocolStartTime = dt.datetime.combine(today, t[0])
    ProtocolEndTime = dt.datetime.combine(today, t[-1])
    
    # Active Trial Analysis
    def TrlMetAnalysis(TrlStr, t, VO2_kg, VO2, VCO2, Mkr):
        today = dt.datetime.today()
        
        # add start and end times
        StartStr = 'start' + TrlStr
        EndStr = 'end' + TrlStr
            
        # pull start and end times of trial
        StartInd = Mkr.index(StartStr)
        EndInd = Mkr.index(EndStr)
        
        # define VO2 during trial
        TrialVO2_kg = VO2_kg[StartInd:EndInd]
        TrialVO2 = VO2[StartInd:EndInd]
        TrialVCO2 = VCO2[StartInd:EndInd]
        TrialW = (np.multiply(np.divide(TrialVO2, 1000), 16.58) + 
                  np.multiply(np.divide(TrialVCO2, 1000), 4.51)) * 1000 / 60 
        
        # create start and end times as datetimes
        StartTime = dt.datetime.combine(today, t[StartInd])
        # EndTime = dt.datetime.combine(today, t[EndInd])
        #  TrialTime = EndTime-StartTime
        
        # convert to seconds
        TrlSec = []
        for i in t[StartInd:EndInd]:
            TrlTime = dt.datetime.combine(today, i)
            ts = TrlTime - StartTime
            TrlSec.append(ts.total_seconds())
        
        # find final 2 min of trial
        Final2Min = [x for x in TrlSec if x >= 180]
        Final2MinInd = TrlSec.index(Final2Min[0])
        
        # average VO2 over final 2 min 
        TrialVO2_kgAvg = np.mean(TrialVO2_kg[Final2MinInd:EndInd])
        TrialVO2Avg = np.mean(TrialVO2[Final2MinInd:EndInd])
        TrialVCO2Avg = np.mean(TrialVCO2[Final2MinInd:EndInd])
        TrialWAvg = np.mean(TrialW[Final2MinInd:EndInd])
        
        VO2Data = {
            'Trial Name' : TrlStr,
            'VO2_kg Avg' : TrialVO2_kgAvg,
            'All VO2_kg Data' : TrialVO2,
            'VO2 Avg' : TrialVO2,
            'All VO2 Data' : TrialVO2Avg,
            'VCO2 Avg' : TrialVCO2Avg,
            'All VCO2 Data' : TrialVCO2,
            'Time Values' : TrlSec,
            'W Avg' : TrialWAvg,
            'All TrialW Data' : TrialW,
            }
        
        return VO2Data
    
    # extract & analyze rest times   
    StartInds = []
    i = 0
    for v in Mkr:
        if 'start' in str(v):
            StartInds.append(i)
        i = i + 1
        
    EndInds = []
    i = 0
    for v in Mkr:
        if 'end' in str(v):
            EndInds.append(i)
        i = i + 1
        
    del(StartInds[0])
    del(EndInds[-1])
    RestTime = []
    for i in range(len(StartInds)):
        starts = dt.datetime.combine(today, t[StartInds[i]])
        ends = dt.datetime.combine(today, t[EndInds[i]])
        NumSec = starts - ends
        RestTime.append(NumSec.seconds)

    Subjects[s]['RestTime'] = RestTime
    
    
    # analyze each active trial individually
    S_M20 = TrlMetAnalysis('S_M20', t, VO2_kg, VO2, VCO2, Mkr)
    S_M10 = TrlMetAnalysis('S_M10', t, VO2_kg, VO2, VCO2, Mkr)
    S_Norm = TrlMetAnalysis('S_Norm', t, VO2_kg, VO2, VCO2, Mkr)
    S_P10 = TrlMetAnalysis('S_P10', t, VO2_kg, VO2, VCO2, Mkr)
    S_P20 = TrlMetAnalysis('S_P20', t, VO2_kg, VO2, VCO2, Mkr)
    
    F_M20 = TrlMetAnalysis('F_M20', t, VO2_kg, VO2, VCO2, Mkr)
    F_M10 = TrlMetAnalysis('F_M10', t, VO2_kg, VO2, VCO2, Mkr)
    F_Norm = TrlMetAnalysis('F_Norm', t, VO2_kg, VO2, VCO2, Mkr)
    F_P10 = TrlMetAnalysis('F_P10', t, VO2_kg, VO2, VCO2, Mkr)
    F_P20 = TrlMetAnalysis('F_P20', t, VO2_kg, VO2, VCO2, Mkr)
        
    # get subject mass
    Subjects[s]['Mass'] = Active.loc[5, 'Unnamed: 1']
    
    # create arrays with all 5 trials
    # net VO2/kg 
    Y_S = [S_M20['VO2_kg Avg'], 
           S_M10['VO2_kg Avg'], 
           S_Norm['VO2_kg Avg'],
           S_P10['VO2_kg Avg'], 
           S_P20['VO2_kg Avg']]
    Y_S_net = Y_S - Subjects[Subj]['RMR VO2_kg Avg']
    Y_F = [F_M20['VO2_kg Avg'], 
           F_M10['VO2_kg Avg'], 
           F_Norm['VO2_kg Avg'],
           F_P10['VO2_kg Avg'], 
           F_P20['VO2_kg Avg'] ]
    Y_F_net = Y_F - Subjects[Subj]['RMR VO2_kg Avg']
    
    # net watts
    TrialW_S = [S_M20['W Avg'],
                S_M10['W Avg'], 
                S_Norm['W Avg'],
                S_P10['W Avg'], 
                S_P20['W Avg']]
    TrialW_S_net = (TrialW_S - Subjects[Subj]['RMR W Avg']) / Subjects[s]['Mass']
    TrialW_F = [F_M20['W Avg'], 
                F_M10['W Avg'], 
                F_Norm['W Avg'],
                F_P10['W Avg'], 
                F_P20['W Avg']]
    TrialW_F_net = (TrialW_F - Subjects[Subj]['RMR W Avg']) / Subjects[s]['Mass']
    

    
    # save in dict
    Subjects[s]['Trial_S_VO2'] = Y_S
    Subjects[s]['Trial_S_VO2net'] = Y_S_net 
    Subjects[s]['Trial_F_VO2'] = Y_F
    Subjects[s]['Trial_F_VO2net'] = Y_F_net
    Subjects[s]['TrialW_S_net'] = TrialW_S_net
    Subjects[s]['TrialW_F_net'] = TrialW_F_net

    # plot VO2 data by condition
    X = [1, 2, 3, 4, 5]
    # plt.plot([x - 0.1 for x in X] ,Y_S_net, '.', color=colorVal, lw=5, ms=MkrSz, alpha=0.6,
    #          label = Subj + ' Speed')
    # plt.plot([x + 0.1 for x in X], Y_F_net, '^', color=colorVal, lw=5, ms=MkrSz, alpha=0.6,
    #          label = Subj + ' Force')
    
    # save subject data to aggregate later
    TrialAvg['Subj'] = Subj
    AvgS.append(TrialW_S_net)
    AvgF.append(TrialW_F_net)

    d = d+1
    
# change shape of output into NxCondition numpy array
WAvg_S = np.reshape(AvgS, [len(Subjects), 5])
WAvg_F = np.reshape(AvgF, [len(Subjects), 5])

# calculate averages
TrialAvg['Trial_S_Wnet'] = np.mean(AvgF,axis=0)
TrialAvg['Trial_S_Wnet_sd'] = np.std(AvgF,axis=0)
TrialAvg['Trial_F_Wnet'] = np.mean(AvgF,axis=0)
TrialAvg['Trial_F_Wnet_sd'] = np.std(AvgF,axis=0)

# boxplot
c = 'red'
box = plt.boxplot(WAvg_S, positions=[x - 0.1 for x in X], 
            widths=0.16, patch_artist=True,
            boxprops=dict(facecolor=c, color=c),
            capprops=dict(color=c),
            whiskerprops=dict(color=c),
            flierprops=dict(color=c, markeredgecolor=c),
            medianprops=dict(color=c),
            )

c = 'c'
plt.boxplot(WAvg_F, positions=[x + 0.1 for x in X], 
            widths=0.16, patch_artist=True,
            boxprops=dict(color=c, facecolor=c),
            capprops=dict(color=c),
            whiskerprops=dict(color=c),
            flierprops=dict(color=c, markeredgecolor=c),
            medianprops=dict(color=c),
            )

# add labels and such to plot
ax.set_title('Metabolic Cost Across All Trials', fontsize=20)
ax.set_xticks([1, 2, 3, 4, 5])
ax.set_xticklabels(['-20','-10','Norm','+10', '+20'], fontsize=15)
plt.text(1,25, '')
plt.xlabel('Trial', fontsize=15)
plt.ylabel('Net W/kg', fontsize=15)   
plt.show()
# plt.savefig('Normalized Trial VO2.jpg', dpi=300)

#%% Rest Time Analysis
RestTimes = []
for s in Subjects:
    for x in Subjects[s]['RestTime']:
        RestTimes.append(x) 
        
AvgRestTimes = np.mean(RestTimes)
SDRestTimes = np.std(RestTimes)
MinRestTimes = np.min(RestTimes)
MaxRestTimes = np.max(RestTimes)

print('Rest Times')
print('Avg = ' + str(AvgRestTimes))
print('SD = ' + str(SDRestTimes))
print('Min = ' + str(MinRestTimes))
print('Max = ' + str(MaxRestTimes))

#%% Load Matlab Data
plt.close('all')
SubjNamez = []
for v in SubjNames:
    SubjNamez.append(v.replace('s','Subj'))

Levels = ['M20', 'M10', 'Norm', 'P10', 'P20']
SubjData = {}

pattern = '*.mat'
matching = fnmatch.filter(files, pattern)
for i in matching:
    for s in SubjNamez:
        if s in i:
            Subj = s
            break 
        
    print('Loading' + s)
    F = folder + '\\' + i
    Dict = {}
    MAT = sio.loadmat(F, mdict = Dict, squeeze_me=1)
    # plt.figure(figsize=[12,10])

    NormSpd = MAT['normSpeed']
    SpdTargets = MAT['speedTargets']
    FpTargets = MAT['FpTargets'][0,:]
    
    def SpdAnalysis(Cond, MAT, Color):
        
        # get variables
        Data = pd.DataFrame()
        Spd = MAT['FpTarget'][Cond]['Data'][:]['Speed'].tolist()
        Time = MAT['FpTarget'][Cond]['Data'][:]['Time'].tolist()
        FpTarget = MAT['FpTarget'][Cond]['TargetFp']
        MeanPkFp = MAT['FpTarget'][Cond]['Data'][:]['MeanPeakFp'].tolist()
        
        # get indicies of final 2 min
        Final2Min = [x for x in Time if x >= 180]
        Final2MinInd = Time.index(Final2Min[Cond])
        EndInd = len(Time)
        
        # calc avg speed over final 2 min
        # Last2 = Time[Final2MinInd:EndInd]
        AvgSpd = np.mean(Spd[Final2MinInd:EndInd])
        A = np.ones(len(Time))
        AvgSpd_array = A * AvgSpd
        
        # cacl avg Fp over final 2 min 
        AvgFp = np.mean(MeanPkFp[Final2MinInd:EndInd])
        AvgFp_array = A * AvgFp
        
        # create array of target speed & Fp
        # f = np.ones(len(Time))
        # SpdTarget = f * SpdTargets[Cond]
        # FpTarget = f * FpTargets[Cond]
    
        
        # plot
        # plt.plot(Time, Spd, '-', c=C, lw=4)
        # plt.plot(Last2, AvgSpd_array, '--', c=C, lw=4)
        # plt.plot(Time, SpdTarget, '-', c=C, lw=2)
        
        # colName = Subj + Levels[Cond]
        Z = np.zeros(len(Time))
        Z[Final2MinInd:EndInd] = 1
        Data['Final2min'] = Z
        Data['F_Time'] = Time
        Data['F_Fp'] = MeanPkFp
        Data['F_AvgFp'] = AvgFp_array
        Data['F_Spd'] = Spd
        Data['F_AvgSpd'] = AvgSpd_array
        Data['F_Target'] = FpTarget
        
        
        # analyze speed targeting trial   
        SpdData = {}
        SpdData['S_Time'] = MAT['NewSpeedTarget'][Cond]['Data'][:]['Time'].tolist()       
        SpdData['S_Fp'] = MAT['NewSpeedTarget'][Cond]['Data'][:]['MeanPeakFp'].tolist()
        
        
        Data['S_AvgFp'] = MAT['SpeedTarget'][Cond]['FpData']['Mean']
        Data['S_AvgSpd'] = MAT['SpeedTarget'][Cond]['Speed']
        
        return Data, SpdData
     
        
    DataM20, SpdDataM20 = SpdAnalysis(0, MAT, 'blue')
    DataM10, SpdDataM10 = SpdAnalysis(1, MAT, 'cornflowerblue')
    DataNorm, SpdDataNorm = SpdAnalysis(2, MAT, 'black')
    DataP10, SpdDataP10 = SpdAnalysis(3, MAT, 'orange')
    DataP20, SpdDataP20 = SpdAnalysis(4, MAT, 'orangered')
    
    
    SubjData.update({s+'M20': DataM20, 
                s+'M10': DataM10, 
                s+'Norm': DataNorm, 
                s+'P10': DataP10, 
                s+'P20': DataP20, 
                s+'SpdM20': SpdDataM20, 
                s+'SpdM10': SpdDataM10, 
                s+'SpdNorm': SpdDataNorm, 
                s+'SpdP10': SpdDataP10, 
                s+'SpdP20': SpdDataP20})

del DataM20, DataM10, DataNorm, DataP10, DataP20, Dict, MAT, F
del SpdDataM20, SpdDataM10, SpdDataNorm, SpdDataP10, SpdDataP20


#%% Sampling Frequency Analysis
SampFreq = []
for s in SubjData:
    t = SubjData[s]['F_Time']
    SampFreq.append(np.mean(np.diff(t))) 
    
SamplingMean = np.mean(SampFreq)
SamplingStd = np.std(SampFreq)

print('Avg Sampling Intermission = ' + str(SamplingMean))
print('Std Sampling Intermission = ' + str(SamplingStd))

#%% Plot Fps during Speed Clamp
sty = 'default'
mpl.style.use(sty)
plt.close('all')
fnt = 16
plt.rcParams.update({'font.size': fnt})
fig = plt.figure(figsize=[20,10])
# ax = plt.subplot(131)
ax = plt.axes((0.06, 0.075, 0.27, 0.75))
from matplotlib.patches import Rectangle
Levels = ['M20', 'M10', 'Norm', 'P10', 'P20']
Lvl = ['-20%', '-10%', 'Norm', '+10%', '+20%']
Colors = [mpl._color_data.CSS4_COLORS['blue'], 
          mpl._color_data.CSS4_COLORS['cornflowerblue'], 
          mpl._color_data.CSS4_COLORS['black'],
          mpl._color_data.CSS4_COLORS['orange'], 
          mpl._color_data.CSS4_COLORS['orangered']]
A1 = 0.15
A2 = 1
LW1 = 2
LW2 = 2
counter = 0
SubjFps = np.zeros([len(SubjNamez), 5])
    
for s in SubjNamez:
    for t in [0, 1, 2, 3, 4]:
    # calculate targeting accuracy
        Time = SubjData[s+'Spd'+Levels[t]]['S_Time']
        NormTarget = SubjData[s+'Norm']['F_Target'][2]
        TrlFp = SubjData[s+'Spd'+Levels[t]]['S_Fp']
        Fp = list(filter(None, TrlFp))
        Ind = TrlFp.index(Fp[0])
        
        plt.plot(Time[Ind:], Fp/NormTarget, 
                  c=Colors[t], alpha=A1)
        
         # get average Fp over final 2 min
        Final2Min = [x for x in Time if x >= 180]
        a = Time.index(Final2Min[0])
        b = len(Time)
        AvgFp = np.ones(len(TrlFp[a:b])) * np.mean(TrlFp[a:b])
        
        SubjFps[counter, t] = AvgFp[0] / NormTarget
    counter = counter + 1
    
for i in range(5):
    ax.add_patch(Rectangle(((180,np.mean(SubjFps[:,i])-0.025)), 120, 0.05,
              edgecolor = Colors[i],
              facecolor = Colors[i],
              alpha = A2,
              fill=True, lw=0))
    val = int(np.mean(SubjFps[:,i])*100)
    font = mpl.font_manager.FontProperties()
    font.set_weight('bold')
    font.set_size(fnt)
    plt.text(240, np.mean(SubjFps[:,i]), 
              Lvl[i] + ' Avg: ' + str(val) + '%', 
              va='center', ha='center', c='w', 
              fontproperties=font)
    
plt.axvline(x = 180, color='k', lw=2)
plt.text(240, 0.62, 'Final 2 Minutes', c='k', 
         fontsize=fnt, ha='center')
plt.xlim([0, 300])
plt.ylim([0.6, 1.4])
ax.set_xticks([0, 100, 200, 300])
ax.set_yticks([0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4])
ax.set_yticklabels(['60%', '70%', '80%', '90%', '100%',
                  '110%', '120%', '130%', '140%'])

plt.xlabel('Time(s)')
plt.ylabel('Relative Fp')
plt.title('Fp during Speed Clamp')
plt.text(-25, 1.42, 'A', fontsize=30)

#%% Plot Fps during Fp Clamp
# ax = plt.subplot(132)
ax = plt.axes((0.39, 0.075, 0.27, 0.75))
from matplotlib.patches import Rectangle
Levels = ['M20', 'M10', 'Norm', 'P10', 'P20']
Lvl = ['-20%', '-10%', 'Norm', '+10%', '+20%']
Colors = [mpl._color_data.CSS4_COLORS['blue'], 
          mpl._color_data.CSS4_COLORS['cornflowerblue'], 
          mpl._color_data.CSS4_COLORS['black'],
          mpl._color_data.CSS4_COLORS['orange'], 
          mpl._color_data.CSS4_COLORS['orangered']]
A1 = 0.15
A2 = 1
LW1 = 2
LW2 = 2
counter = 0
SubjFps = np.zeros([len(SubjNamez), 5])
def MovingAvg(Vals, Window):
    Filt = Vals
    Win = int((Window-1) / 2)
    for x in np.arange(Win,len(Vals)-Win):
        Filt[x] = np.nanmean([Vals[x-Win:x+Win]])
        
    return Filt
    
for s in SubjNamez:
    for t in [0, 1, 2, 3, 4]:
    # calculate targeting accuracy
        Time = SubjData[s+Levels[t]]['F_Time'].to_list()
        NormTarget = SubjData[s+'Norm']['F_Target'][2]
        TrlFp = SubjData[s+Levels[t]]['F_Fp'].values.tolist() 
        Fp = list(filter(None, TrlFp))
        Ind = TrlFp.index(Fp[0])
       
        plt.plot(Time[Ind:], Fp/NormTarget, 
                  c=Colors[t], alpha=A1)
        
         # get average Fp over final 2 min
        Final2Min = [x for x in Time if x >= 180]
        a = Time.index(Final2Min[0])
        b = len(Time)
        AvgFp = np.ones(len(TrlFp[a:b])) * np.mean(TrlFp[a:b])
        
        SubjFps[counter, t] = AvgFp[0] / NormTarget
    counter = counter + 1
    
for i in range(5):
    ax.add_patch(Rectangle(((180,np.mean(SubjFps[:,i])-0.025)), 120, 0.05,
              edgecolor = Colors[i],
              facecolor = Colors[i],
              alpha = A2,
              fill=True, lw=0))
    val = int(np.mean(SubjFps[:,i])*100)
    font = mpl.font_manager.FontProperties()
    font.set_weight('bold')
    font.set_size(fnt)
    plt.text(240, np.mean(SubjFps[:,i]), 
              Lvl[i] + ' Avg: ' + str(val) + '%', 
              va='center', ha='center', c='w', 
              fontproperties=font)
    
plt.axvline(x = 180, color='k', lw=2)
plt.text(240, 0.62, 'Final 2 Minutes', c='k', 
         fontsize=fnt, ha='center')
plt.xlim([0, 300])
plt.ylim([0.6, 1.4])
ax.set_xticks([0, 100, 200, 300])
ax.set_yticks([0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4])
ax.set_yticklabels(['60%', '70%', '80%', '90%', '100%',
                  '110%', '120%', '130%', '140%'])

plt.xlabel('Time(s)')
plt.ylabel('Relative Fp')
plt.title('Fp during Fp Clamp')
plt.text(-25, 1.42, 'B', fontsize=30)

#%% Plot General Speeds
# from matplotlib.patches import Rectangle
Levels = ['M20', 'M10', 'Norm', 'P10', 'P20']
Lvl = ['-20%', '-10%', 'Norm', '+10%', '+20%']
Colors = [mpl._color_data.CSS4_COLORS['blue'], 
          mpl._color_data.CSS4_COLORS['cornflowerblue'], 
          mpl._color_data.CSS4_COLORS['black'],
          mpl._color_data.CSS4_COLORS['orange'], 
          mpl._color_data.CSS4_COLORS['orangered']]

# ax = plt.subplot(133)
ax = plt.axes((0.72, 0.075, 0.27, 0.75))
A1 = 0.3
A2 = 1
LW1 = 2
LW2 = 2
SubjSpds = np.zeros([len(SubjNamez), 5])
counter = 0
for s in SubjNamez:
    
    NormSpeed = SubjData[s+'Norm']['S_AvgSpd'][0]
    # time = SubjData[s+'Norm']['F_Time'].to_list()
    NormSpd = SubjData[s+'Norm']['F_Spd'].values / NormSpeed
    
    for t in [0, 1, 2, 3, 4]:
    
        # plot individual speed lines for final 4 min
        time = SubjData[s+Levels[t]]['F_Time'].to_list()
        TrlSpd = SubjData[s+Levels[t]]['F_Spd'].values / NormSpeed

        plt.plot(time, TrlSpd, lw=LW1,
                 c=Colors[t], alpha=A1)
        
        # plot final 2 min average speeds
        Final2Min = [x for x in time if x >= 180]
        a = time.index(Final2Min[0])
        b = len(time)
        AvgSpd = np.ones(len(TrlSpd[a:b])) * np.mean(TrlSpd[a:b])
        SubjSpds[counter, t] = AvgSpd[0]
    counter = counter + 1

for i in range(5):
    ax.add_patch(Rectangle(((180,np.mean(SubjSpds[:,i])-0.025)), 120, 0.05,
             edgecolor = Colors[i],
             facecolor = Colors[i],
             alpha = A2,
             fill=True, lw=0))
    val = int(np.mean(SubjSpds[:,i])*100)
    font = mpl.font_manager.FontProperties()
    font.set_weight('bold')
    font.set_size(fnt)
    plt.text(240, np.mean(SubjSpds[:,i]), 
             Lvl[i] + ' Avg: ' + str(val) + '%', 
             va='center', ha='center', c='w', 
             fontproperties=font)
    
plt.axvline(x = 180, color='k', lw=2)
plt.text(240, 0.62, 'Final 2 Minutes', c='k', 
         fontsize=fnt, ha='center')
plt.xlim([0, 300])
plt.ylim([0.6, 1.4])
ax.set_xticks([0, 100, 200, 300])
ax.set_yticks([0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4])
ax.set_yticklabels(['60%', '70%', '80%', '90%', '100%',
                  '110%', '120%', '130%', '140%'])

plt.xlabel('Time(s)')
plt.ylabel('Relative Walking Speed')
plt.text(-25, 1.42, 'C', fontsize=30)
plt.title('Walking Speed during Fp Clamp')

plt.savefig('BiofeedbackPerformance.jpg', dpi=300)
plt.savefig('BiofeedbackPerformance.pdf', dpi=300)
    

#%% plot Speed and Fp from each trial 
plt.close('all')
fig = plt.figure(figsize=[12,12])
Conditions = ['M20', 'M10', 'P10', 'P20']
Colors = [mpl._color_data.CSS4_COLORS['blue'], 
          mpl._color_data.CSS4_COLORS['cornflowerblue'], 
          mpl._color_data.CSS4_COLORS['orange'], 
          mpl._color_data.CSS4_COLORS['orangered']]
TrialInd = [0, 1, 3, 4]
AllSpd_F = []
AllFp_F = []
AllSpd_S = []
AllFp_S = []
N = len(Subjects)
Ones = np.ones(N)
Mass = [0]*N
MassTxt = [0]*N
Alls = [0]*N
C = 0
for s in Subjects:
    subjSpd_S = [0, 0, 0, 0, 0]
    subjFp_S = [0, 0, 0, 0, 0]
    subjSpd_F = [0, 0, 0, 0, 0]
    subjFp_F = [0, 0, 0, 0, 0]
    
    Mass[C] = Subjects[s]['Mass']
    MassTxt[C] = s

    Subj = s.replace('s', 'Subj')
    Trial = 'Norm'
    Key1 = Subj + Trial
    # Key2 = Subj + Trial
    # generate norms
    NormSpd_F = SubjData[Key1]['F_AvgSpd'][0]
    NormFp_F = SubjData[Key1]['F_AvgFp'][0]
    NormSpd_S = SubjData[Key1]['S_AvgSpd'][0]
    NormFp_S = SubjData[Key1]['S_AvgFp'][0]
    
    subjSpd_S[2] = NormSpd_S
    subjFp_S[2] = NormFp_S
    subjSpd_F[2] = NormSpd_F
    subjFp_F[2] = NormFp_F
    
    ax1 = fig.add_subplot(221)
    plt.scatter(1, Subjects[s]['TrialW_S_net'][2],
                    c='k', marker='.')
    ax2 = fig.add_subplot(222)
    plt.scatter(1, Subjects[s]['TrialW_S_net'][2],
                    c='k', marker='.')
    
    ax3 = fig.add_subplot(223)
    plt.scatter(1, Subjects[s]['TrialW_F_net'][2],
                    c='k', marker='.')
    ax4 = fig.add_subplot(224)
    plt.scatter(1, Subjects[s]['TrialW_F_net'][2],
                    c='k', marker='.')
    
    # loop through non-norm conditions
    for cond in [0, 1, 2, 3]:
        Trial = Conditions[cond]
        Key1 = Subj + Trial
        # Key2 = Subj + Trial
        
        # calculate & normalize variables
        Spd_S = SubjData[Key1]['S_AvgSpd'][0]
        Fp_S = SubjData[Key1]['S_AvgFp'][0]
        Spd_F = SubjData[Key1]['F_AvgSpd'][0]
        Fp_F = SubjData[Key1]['F_AvgFp'][0]

        subjSpd_S[TrialInd[cond]] = Spd_S
        subjFp_S[TrialInd[cond]] = Fp_S
        subjSpd_F[TrialInd[cond]] = Spd_F
        subjFp_F[TrialInd[cond]] = Fp_F
        
          # plot values being sure to normalize to Norm trial
        ax1 = fig.add_subplot(221)
        plt.scatter(Spd_S/NormSpd_S,
                    Subjects[s]['TrialW_S_net'][TrialInd[cond]], 
                    c=Colors[cond], marker='.', alpha=0.5)
        ax2 = fig.add_subplot(222)
        plt.scatter(Fp_S/NormFp_S,
                    Subjects[s]['TrialW_S_net'][TrialInd[cond]], 
                    c=Colors[cond], marker='.', alpha=0.5)
        

        ax3 = fig.add_subplot(223)
        plt.scatter(Spd_F/NormSpd_F, 
                        Subjects[s]['TrialW_F_net'][TrialInd[cond]], 
                        c=Colors[cond], marker='.', alpha=0.5)
        ax4 = fig.add_subplot(224)
        plt.scatter(Fp_F/NormFp_F, 
                        Subjects[s]['TrialW_F_net'][TrialInd[cond]], 
                        c=Colors[cond], marker='.', alpha=0.5)
        
        
    AllSpd_S.append(subjSpd_S)
    AllFp_S.append(subjFp_S)
    AllSpd_F.append(subjSpd_F)
    AllFp_F.append(subjFp_F)
    Alls[C] = Subj
    C = C+1
    
    
AllSpd_S = np.reshape(AllSpd_S, [len(Subjects), 5])
AllFp_S = np.reshape(AllFp_S, [len(Subjects), 5])
AllSpd_F = np.reshape(AllSpd_F, [len(Subjects), 5])
AllFp_F = np.reshape(AllFp_F, [len(Subjects), 5])

CoT_Fp_S = WAvg_S / AllSpd_S
CoT_Fp_F = WAvg_F / AllSpd_F
Fp_S = np.array(AllFp_S)
Fp_F = np.array(AllFp_S)

for i in range(len(Subjects)):
    Fp_S[i,:] = AllFp_S[i,:] / Mass[i]
    Fp_F[i,:] = AllFp_F[i,:] / Mass[i]

Fp_S = np.reshape(Fp_S, [len(Subjects), 5])
Fp_F = np.reshape(Fp_F, [len(Subjects), 5])

    
ax1.set_xlabel('Normalized Speed', fontsize=10)
ax1.set_ylabel('Net W/kg', fontsize=10)
ax1.set_title('Speed during Fixed Speed', fontsize=12)

ax2.set_xlabel('Normalized Fp', fontsize=10)
# ax2.set_ylabel(ax, 'Net W/kg', fontsize=10)
ax2.set_title('Fp during Fixed Speed', fontsize=12)

ax3.set_xlabel('Normalized Speed', fontsize=10)
ax3.set_ylabel('Net W/kg', fontsize=10)
ax3.set_title('Speed during Fp Targeting', fontsize=12)

ax4.set_xlabel('Normalized Fp', fontsize=10)
# ax4.set_ylabel(ax, 'Net W/kg', fontsize=10)
ax4.set_title('Fp during Fp Targeting', fontsize=12)

# plt.savefig('SpeedsFps.jpg', dpi=300)


#%% Speed Between Fixed and Targeting
plt.close('all')
Conditions = ['-20', '-10', 'Norm', '+10', '+20']
Colors = [mpl._color_data.CSS4_COLORS['blue'], 
          mpl._color_data.CSS4_COLORS['cornflowerblue'], 
          'k',
          mpl._color_data.CSS4_COLORS['orange'], 
          mpl._color_data.CSS4_COLORS['orangered']]
N = len(Subjects)
Ones = np.ones(N)
Mass = [0]*N
Ofst = 0.1
BarOfst = 0.2
Trans = 0.4
Trans2 = 1
MkrSz = 16
Fnt = 12
TFnt = 16

fig = plt.figure(figsize=[12,12])
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)
AllFp_S_kg = np.zeros_like(AllSpd_S)
AllFp_F_kg = np.zeros_like(AllSpd_S)

for x in range(5):
    for i in range(N):
        ax1.plot([X[x]-Ofst, X[x]+Ofst], 
                  [AllSpd_S[i, x], AllSpd_F[i,x]],
                  '-', c=Colors[x], alpha=Trans)
        
        Mass[i] = Subjects[SubjNames[i]]['Mass']
        AllFp_S_kg[i, x] = AllFp_S[i, x] / Mass[i]
        AllFp_F_kg[i, x] = AllFp_F[i, x] / Mass[i]
        
        ax2.plot([X[x]-Ofst, X[x]+Ofst], 
                  [AllFp_S[i, x] / Mass[i], AllFp_F[i,x] / Mass[i]],
                  '-', c=Colors[x], alpha=Trans)
        
        ax3.plot(AllSpd_S[i, x], AllSpd_F[i,x], 
                  '.', c=Colors[x], alpha=Trans2)
        
        ax4.plot(AllFp_S[i, x] / Mass[i], AllFp_F[i,x] / Mass[i], 
                  '.', c=Colors[x], alpha=Trans2)
        
    # plot group averages
    ax1.errorbar(X[x]-BarOfst, np.mean(AllSpd_S[:, x], axis=0), 
                  yerr=np.std(AllSpd_S[:, x], axis=0),      
                  marker='.', c=Colors[x], ecolor=Colors[x], markersize=MkrSz)
    ax1.errorbar(X[x]+BarOfst, np.mean(AllSpd_F[:, x], axis=0), 
                  yerr=np.std(AllSpd_F[:, x], axis=0),      
                  marker='^', c=Colors[x], ecolor=Colors[x], markersize=MkrSz)
    ax2.errorbar(X[x]-BarOfst, np.mean(AllFp_S_kg[:, x], axis=0), 
                  yerr=np.std(AllFp_S_kg[:, x], axis=0),      
                  marker='.', c=Colors[x], ecolor=Colors[x], markersize=MkrSz)
    ax2.errorbar(X[x]+BarOfst, np.mean(AllFp_F_kg[:, x], axis=0), 
                  yerr=np.std(AllFp_F_kg[:, x], axis=0),      
                  marker='^', c=Colors[x], ecolor=Colors[x], markersize=MkrSz)
    
ax1.set_xlabel('Condition', fontsize=Fnt)
ax1.set_ylabel('m/s', fontsize=Fnt)
ax1.set_title('Speed Across Conditions', fontsize=TFnt)
ax1.set_xticks(X)
ax1.set_xticklabels(Conditions)

ax2.set_xlabel('Condition', fontsize=Fnt)
ax2.set_ylabel('N / kg', fontsize=Fnt)
ax2.set_title('Fp Across Conditions', fontsize=TFnt)
ax2.set_xticks(X)
ax2.set_xticklabels(Conditions)

ax3.set_xlabel('Fixed Speed (m/s)', fontsize=Fnt)
ax3.set_ylabel('Fp Targeting (m/s)', fontsize=Fnt)
ax3.set_title('Speed by Condition', fontsize=TFnt)

ax4.set_xlabel('Fixed Speed (N/kg)', fontsize=Fnt)
ax4.set_ylabel('Fp Targeting (N/kg)', fontsize=Fnt)
ax4.set_title('Fp by Condition', fontsize=TFnt)

plt.savefig('SpeedsFpComp.jpg', dpi=300)


#%% Abstract Plot (Speed, Fp, and CoT)
plt.close('all')
Conditions = ['-20%', '-10%', 'Norm', '+10%', '+20%']
Colors = [mpl._color_data.CSS4_COLORS['blue'], 
          mpl._color_data.CSS4_COLORS['cornflowerblue'], 
          'k',
          mpl._color_data.CSS4_COLORS['orange'], 
          mpl._color_data.CSS4_COLORS['orangered']]
N = len(Subjects)
Ones = np.ones(N)
# Mass = [0]*N
Ofst = 0.12
Trans = 0.25
Full = 1
MkrSz = 10
MkrSz2 = 14
fnt = 15

fig = plt.figure(figsize=[12,12])
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)

for x in range(5):
        
    # 1st plot - speed
    ax1.plot([X[x]-Ofst, X[x]+Ofst], 
             [AllSpd_S[:, x], AllSpd_F[:,x]],
             '-', c=Colors[x], alpha=Trans)
    ax1.errorbar(X[x]-BarOfst, np.mean(AllSpd_S[:, x], axis=0), 
                 yerr=np.std(AllSpd_S[:, x], axis=0),      
                 marker='o', c=Colors[x], ecolor=Colors[x], markersize=MkrSz2)
    ax1.errorbar(X[x]+BarOfst, np.mean(AllSpd_F[:, x], axis=0), 
                 yerr=np.std(AllSpd_F[:, x], axis=0),      
                 marker='s', c=Colors[x], ecolor=Colors[x], markersize=MkrSz2)
    
    # 2nd plot - Fp
    ax2.plot([X[x]-Ofst, X[x]+Ofst], 
             [AllFp_S_kg[:, x]*9.81, AllFp_F_kg[:,x]*9.81],
             '-', c=Colors[x], alpha=Trans)
    ax2.errorbar(X[x]-BarOfst, np.mean(AllFp_S_kg[:, x]*9.81, axis=0), 
                 yerr=np.std(AllFp_S_kg[:, x]*9.81, axis=0),      
                 marker='o', c=Colors[x], ecolor=Colors[x], markersize=MkrSz2)
    ax2.errorbar(X[x]+BarOfst, np.mean(AllFp_F_kg[:, x]*9.81, axis=0), 
                 yerr=np.std(AllFp_F_kg[:, x]*9.81, axis=0),      
                 marker='s', c=Colors[x], ecolor=Colors[x], markersize=MkrSz2)
    
    # 3rd plot - Met Cost
    ax3.plot([X[x]-Ofst, X[x]+Ofst], 
             [WAvg_S[:,x], WAvg_F[:,x]],
             '-', c=Colors[x], alpha=Trans)
    ax3.errorbar(X[x]-BarOfst, np.mean(WAvg_S[:,x], axis=0), 
         yerr=np.std(WAvg_S[:,x], axis=0), 
         marker='o', ecolor=Colors[x], color=Colors[x], markersize=MkrSz2)
    ax3.errorbar(X[x]+BarOfst, np.mean(WAvg_F[:,x], axis=0), 
         yerr=np.std(WAvg_F[:,x], axis=0), 
         marker='s', ecolor=Colors[x], color=Colors[x], markersize=MkrSz2)
    
    # 4th plot - CoT
    ax4.plot([X[x]-Ofst, X[x]+Ofst], 
             [CoT_Fp_S[:,x], CoT_Fp_F[:,x]],
             '-', c=Colors[x], alpha=Trans)
    ax4.errorbar(X[x]-BarOfst, np.mean(CoT_Fp_S[:,x], axis=0), 
         yerr=np.std(CoT_Fp_S[:,x], axis=0), 
         marker='o', ecolor=Colors[x], color=Colors[x], markersize=MkrSz2)
    ax4.errorbar(X[x]+BarOfst, np.mean(CoT_Fp_F[:,x], axis=0), 
         yerr=np.std(CoT_Fp_F[:,x], axis=0), 
         marker='s', ecolor=Colors[x], color=Colors[x], markersize=MkrSz2)
    
    
# create legend
ax1.plot(1, 1.95, 'o', color='k', markersize=MkrSz2)
ax1.plot(1, 1.85, 's', color='k', markersize=MkrSz2)
ax1.text(1.12, 1.95, ' = Speed Clamp', color='k', fontsize=18, va='center')
ax1.text(1.12, 1.85, ' = Fp Clamp', color='k', fontsize=18, va='center')
ax1.text(0, 2.08, 'A', fontsize=fnt*2)
ax2.text(0, 31, 'B', fontsize=fnt*2)
ax3.text(0, 10.5, 'C', fontsize=fnt*2)
ax4.text(0, 6.2, 'D', fontsize=fnt*2)

# edit axes
# ax1.set_xlabel('Condition', fontsize=15)
ax1.set_ylabel('Walking Speed (m/s)', fontsize=15)
# ax1.set_title('A', fontsize=20, horizontalalignment='left')
ax1.set_xticks(X)
ax1.set_xticklabels(Conditions, fontsize=15)
ax1.tick_params(axis='y', labelsize=15) 
# plt.title(label='A', fontsize=20, Loc='left')
# ax1.set_ylim(1, 2)

# ax2.set_xlabel('Condition', fontsize=15)
ax2.set_ylabel('Fp (% body weight)', fontsize=15)
# ax2.set_title('B', fontsize=20, horizontalalignment='left')
ax2.set_xticks(X)
ax2.set_xticklabels(Conditions, fontsize=15)
ax2.tick_params(axis='y', labelsize=15) 
# ax2.set_ylim(1.3, 3)


ax3.plot(1, 9.2, 'o', color='k', markersize=MkrSz2)
ax3.plot(1, 8.2, 's', color='k', markersize=MkrSz2)
ax3.text(1.12, 9.2, ' = Speed Clamp', color='k', fontsize=18, va='center')
ax3.text(1.12, 8.2, ' = Fp Clamp', color='k', fontsize=18, va='center')

ax3.set_xticks(X)
ax3.set_xticklabels(Conditions, fontsize=15)
ax3.set_ylabel('Net Metabolic Power (W/kg)', fontsize=15)
# ax3.set_title('C', fontsize=20, Loc='left')
ax3.tick_params(axis='y', labelsize=15) 
ax3.tick_params(axis='x', labelsize=15)
# ax3.set_xlim(1.3, 3)
ax3.set_ylim(0, 10)

ax4.set_xticks(X)
ax4.set_xticklabels(Conditions, fontsize=15)
ax4.set_ylabel('CoT (J/kg/m)', fontsize=15)
# ax4.set_title('D', fontsize=20, Loc='left')
ax4.tick_params(axis='y', labelsize=15) 
ax4.set_ylim(1, 6)
ax4.tick_params(axis='x', labelsize=15)


#%% Run stats and add to plot
S = range(1,len(Subjects)+1)
Ones = np.ones(5)
fnt = 12
RMA = pd.DataFrame({'subjects': np.tile(np.repeat(S, len(X)), 2),
                   'condition': np.tile(X, len(Subjects)*2),
                   'clamp': np.repeat(np.hstack((Ones, Ones*2)), len(Subjects)),
                   'speed': np.reshape([AllSpd_S, AllSpd_F],
                                       [len(Subjects)*2*5, 1][0]),
                   'Fp': np.reshape([AllFp_S_kg*9.81, AllFp_F_kg*9.81],
                                    [len(Subjects)*2*5, 1][0]), 
                   'MetCost': np.reshape([WAvg_S, WAvg_F],
                                         [len(Subjects)*2*5, 1][0]), 
                   'CoT': np.reshape([CoT_Fp_S, CoT_Fp_F],
                                     [len(Subjects)*2*5, 1][0])}
                    )

AnovaNames = ['speed', 'Fp', 'MetCost', 'CoT']
aov = {}
for A in AnovaNames:
    aov[A] = pg.rm_anova(data=RMA, dv=A, within=['condition', 'clamp'], subject='subjects', detailed=True)
    print('\n\n' + A + '\n')
    print('P values: ')
    print(aov[A]['p-unc'])
    print('Partial Eta Sq: ')
    print(aov[A]['np2'])
    

#perform the repeated measures ANOVA
# print('Speeds')
# print(AnovaRM(data=RMA, depvar='speed', subject='subjects', 
#               within=['condition', 'clamp']).fit())
# AnovaRM(data=RMA, depvar='speed', subject='subjects', 
#               within=['condition', 'clamp'])


# #perform the repeated measures ANOVA
# print('Fps')
# print(AnovaRM(data=RMA, depvar='Fp', subject='subjects', 
#               within=['condition', 'clamp']).fit())

# # Met Cost repeated measures ANOVA
# print('MetCost')
# print(AnovaRM(data=RMA, depvar='MetCost', subject='subjects', 
#               within=['condition', 'clamp']).fit())

# # CoT repeated measures ANOVA
# print('CoT')
# print(AnovaRM(data=RMA, depvar='CoT', subject='subjects', 
#               within=['condition', 'clamp']).fit())

# place ANOVA values in fig
ax1.text(3.8, 1.15,'condition'+'\n'+'clamp'+'\n'+'condition x clamp', 
         va='top', fontsize = fnt, ha='right')
ax1.text(4.0, 1.19,'     p'+'\n'+'<0.001'+'\n'+'  0.004'+'\n'+'<0.001', 
         va='top', fontsize = fnt, ha='left')
np2 = np.round(aov['speed']['np2'].to_list(), 3)
ax1.text(4.75, 1.225,'    $\eta^2_p$\n'+str(np2[0])+'\n'+str(np2[1])+'\n'+str(np2[2]), 
         va='top', fontsize = fnt, ha='left')


ax2.text(3.8, 17.4,'condition'+'\n'+'clamp'+'\n'+'condition x clamp', 
         va='top', fontsize = fnt, ha='right')
ax2.text(4, 18,'     p'+'\n'+'<0.001'+'\n'+'  0.001'+'\n'+'<0.001', 
         va='top', fontsize = fnt, ha='left')
np2 = np.round(aov['Fp']['np2'].to_list(), 3)
ax2.text(4.75, 18.5,'    $\eta^2_p$\n'+str(np2[0])+'\n'+str(np2[1])+'\n'+str(np2[2]), 
         va='top', fontsize = fnt, ha='left')

ax3.text(3.8, 1.675,'condition'+'\n'+'clamp'+'\n'+'condition x clamp', 
          va='top', fontsize = fnt, ha='right')
ax3.text(4, 2,'     p'+'\n'+'<0.001'+'\n'+'  0.002'+'\n'+'  0.126', 
          va='top', fontsize = fnt, ha='left')
np2 = np.round(aov['MetCost']['np2'].to_list(), 3)
ax3.text(4.75, 2.325,'    $\eta^2_p$\n'+str(np2[0])+'\n'+str(np2[1])+'\n'+str(np2[2]), 
         va='top', fontsize = fnt, ha='left')

ax4.text(3.8, 1.8,'condition'+'\n'+'clamp'+'\n'+'condition x clamp', 
         va='top', fontsize = fnt, ha='right')
ax4.text(4, 2,'     p'+'\n'+'<0.001'+'\n'+'  0.010'+'\n'+'  0.313', 
         va='top', fontsize = fnt, ha='left')
np2 = np.round(aov['CoT']['np2'].to_list(), 3)
ax4.text(4.75, 2.175,'    $\eta^2_p$\n'+str(np2[0])+'\n'+str(np2[1])+'\n'+str(np2[2]), 
         va='top', fontsize = fnt, ha='left')

#%% Post hoc T-tests
Ast = 26
Hash = 18

# speed sub-analysis
# test across conditions
T_SCondSpeed = np.ones(5)
T_FCondSpeed = np.ones(5)
ES_SCondSpeed = np.ones(5)
ES_FCondSpeed = np.ones(5)
G = np.array([np.ones(20), 2*np.ones(20)]).reshape([40, 1])
for x in [0, 1, 3, 4]:
    d = np.reshape([AllSpd_S[:,x], AllSpd_S[:,2]], [40, 1])
    D = pd.DataFrame(np.hstack([d,G]), columns=['X','G'])
    Stats = pg.pairwise_tukey(D, dv='X', between='G', effsize='eta-square')
    T_SCondSpeed[x] = float(Stats['p-tukey'])
    ES_SCondSpeed[x] = float(Stats['eta-square'])
    if T_SCondSpeed[x] < 0.05 :
        ax1.text(x+1-BarOfst, np.mean(AllSpd_S[:,x])+0.18, '*', 
                 c = Colors[x], fontsize=Ast, ha='center')
        
    d = np.reshape([AllSpd_F[:,x], AllSpd_F[:,2]], [40, 1])
    D = pd.DataFrame(np.hstack([d,G]), columns=['X','G'])
    Stats = pg.pairwise_tukey(D, dv='X', between='G', effsize='eta-square')
    T_FCondSpeed[x] = float(Stats['p-tukey'])
    ES_FCondSpeed[x] = float(Stats['eta-square'])
    if T_FCondSpeed[x] < 0.05 :
        ax1.text(x+1+BarOfst, np.mean(AllSpd_F[:,x])+0.18, '*', 
                 c = Colors[x], fontsize=Ast, ha='center')
        
    
print('\nSpeed Speed Conditions Post Hoc')
print('p-values')
print(T_SCondSpeed)
print('effect sizes')
print(np.round(ES_SCondSpeed, decimals=5))
print('Speed Fp Conditions Post Hoc')
print('p-values')
print(np.round(T_FCondSpeed, decimals=5))
print('effect sizes')
print(np.round(ES_FCondSpeed, decimals=5))

# between clamps
T_Speed = np.ones(5)
ES_Speed = np.ones(5)
for x in range(5):
    d = np.reshape([AllSpd_S[:,x], AllSpd_F[:,x]], [40, 1])
    D = pd.DataFrame(np.hstack([d,G]), columns=['X','G'])
    Stats = pg.pairwise_tukey(D, dv='X', between='G', effsize='eta-square')
    T_Speed[x] = float(Stats['p-tukey'])
    ES_Speed[x] = float(Stats['eta-square'])
    if T_Speed[x] < 0.05 :
        y = np.mean([np.mean(AllSpd_S[:,x], axis=0), np.mean(AllSpd_F[:,x], axis=0)])
        ax1.text(x+1, y+0.2, '#', 
                 c = Colors[x], fontsize=Hash, ha='center')

print('Speed Between Clamp Post Hoc')
print('p-values')
print(np.round(T_Speed, decimals=5))
print('effect sizes')
print(np.round(ES_Speed, decimals=5))
print(' ')
print(' ')
        
# Fp sub-analysis
# between conditions
T_SCondFp = np.ones(5)
T_FCondFp = np.ones(5)
ES_SCondFp = np.ones(5)
ES_FCondFp = np.ones(5)
for x in [0, 1, 3, 4]:
    d = np.reshape([AllFp_S_kg[:,x]*9.81, AllFp_S_kg[:,2]*9.81], [40, 1])
    D = pd.DataFrame(np.hstack([d,G]), columns=['X','G'])
    Stats = pg.pairwise_tukey(D, dv='X', between='G', effsize='eta-square')
    T_SCondFp[x] = float(Stats['p-tukey'])
    ES_SCondFp[x] = float(Stats['eta-square'])
    if T_SCondFp[x] < 0.05 :
        ax2.text(x+1-BarOfst, np.mean(AllFp_S_kg[:,x]*9.81)+2.5, '*', 
                 c = Colors[x], fontsize=Ast, ha='center')

    d = np.reshape([AllFp_F_kg[:,x]*9.81, AllFp_F_kg[:,2]*9.81], [40, 1]) 
    D = pd.DataFrame(np.hstack([d,G]), columns=['X','G'])
    Stats = pg.pairwise_tukey(D, dv='X', between='G', effsize='eta-square')
    T_FCondFp[x] = float(Stats['p-tukey'])
    ES_FCondFp[x] = float(Stats['eta-square'])
    if T_FCondFp[x] < 0.05 :
        ax2.text(x+1+BarOfst, np.mean(AllFp_F_kg[:,x]*9.81, axis=0)+2.5, '*', 
                 c = Colors[x], fontsize=Ast, ha='center')
        
print('Fp Speed Conditions Post Hoc')
print('p-values')
print(T_SCondFp)
print('effect sizes')
print(np.round(ES_SCondFp, decimals=5))
print('Fp Fp Conditions Post Hoc')
print('p-values')
print(T_FCondFp)
print('effect sizes')
print(np.round(ES_FCondFp, decimals=5))

# between clamps
T_Fp = np.ones(5)
ES_Fp = np.ones(5)
for x in range(5):
    d = np.reshape([AllFp_S_kg[:,x]*9.81, AllFp_F_kg[:,x]*9.81], [40, 1])
    D = pd.DataFrame(np.hstack([d,G]), columns=['X','G'])
    Stats = pg.pairwise_tukey(D, dv='X', between='G', effsize='eta-square')
    T_Fp[x] = float(Stats['p-tukey'])
    ES_Fp[x] = float(Stats['eta-square'])
    if T_Fp[x] < 0.05 :
        y = np.mean([np.mean(AllFp_S_kg[:,x]*9.81, axis=0), 
                     np.mean(AllFp_F_kg[:,x]*9.81, axis=0)])
        ax2.text(x, y+3, '#', c = Colors[x], fontsize=Hash, ha='center')
print('Fp Between Clamp Post Hoc')
print('p-values')
print(np.round(T_Speed, decimals=5))
print('effect sizes')
print(np.round(ES_Speed, decimals=5))
print(' ')
print(' ')

# MetCost sub-analysis
# between conditions
T_SCondMetCost = np.ones(5)
T_FCondMetCost = np.ones(5)
ES_SCondMetCost = np.ones(5)
ES_FCondMetCost = np.ones(5)
for x in range(5): 
    d = np.reshape([WAvg_S[:,x], WAvg_S[:,2]], [40, 1])
    D = pd.DataFrame(np.hstack([d,G]), columns=['X','G'])
    Stats = pg.pairwise_tukey(D, dv='X', between='G', effsize='eta-square')
    T_SCondMetCost[x] = float(Stats['p-tukey'])
    ES_SCondMetCost[x] = float(Stats['eta-square'])
    if T_SCondMetCost[x] < 0.05 :
        ax3.text(x+1-BarOfst, np.mean(WAvg_S[:,x], axis=0)+1.5,
                 '*', c = Colors[x], fontsize=Ast, ha='center')
    
    d = np.reshape([WAvg_F[:,x], WAvg_F[:,2]], [40, 1]) 
    D = pd.DataFrame(np.hstack([d,G]), columns=['X','G'])
    Stats = pg.pairwise_tukey(D, dv='X', between='G', effsize='eta-square')
    T_FCondMetCost[x] = float(Stats['p-tukey'])
    ES_FCondMetCost[x] = float(Stats['eta-square'])
    if T_FCondMetCost[x] < 0.05 :
        ax3.text(x+1+BarOfst, np.mean(WAvg_F[:,x], axis=0)+1.5, 
                 '*', c = Colors[x], fontsize=Ast, ha='center')
        
print('MetCost Speed Conditions Post Hoc')
print('p-values')
print(T_SCondMetCost)
print('effect sizes')
print(np.round(ES_SCondMetCost, decimals=5))
print('MetCost Conditions Post Hoc')
print('p-values')
print(T_FCondMetCost)
print('effect sizes')
print(np.round(ES_FCondMetCost, decimals=5))

# post hoc difference in net metabolic cost for lowest condition intensity
# np.mean([WAvg_F[:,0]]) - np.mean([WAvg_S[:,0]])

# between clamps
T_MetCost = np.ones(5)
ES_MetCost = np.ones(5)
for x in range(5):
    d = np.reshape([WAvg_S[:,x], WAvg_F[:,x]], [40, 1])    
    D = pd.DataFrame(np.hstack([d,G]), columns=['X','G'])
    Stats = pg.pairwise_tukey(D, dv='X', between='G', effsize='eta-square')
    T_MetCost[x] = float(Stats['p-tukey'])
    ES_MetCost[x] = float(Stats['eta-square'])
    if T_MetCost[x] < 0.05 :
        y = np.mean([np.mean(WAvg_S[:,x], axis=0), 
                     np.mean(WAvg_F[:,x], axis=0)])
        ax3.text(x+1, y+1.7, '#', 
                 c = Colors[x], fontsize=Hash, ha='center')

print('MetCost Between Clamp Post Hoc')
print('p-values')
print(np.round(T_MetCost, decimals=5))
print('effect sizes')
print(np.round(ES_MetCost, decimals=5))
print(' ')
print(' ')

# CoT sub-analysis
# between conditions
T_SCondCoT = np.ones(5)
T_FCondCoT = np.ones(5)
ES_SCondCoT = np.ones(5)
ES_FCondCoT = np.ones(5)
for x in range(5):
    d = np.reshape([CoT_Fp_S[:,x], CoT_Fp_S[:,2]], [40, 1]) 
    D = pd.DataFrame(np.hstack([d,G]), columns=['X','G'])
    Stats = pg.pairwise_tukey(D, dv='X', between='G', effsize='eta-square')
    T_SCondCoT[x] = float(Stats['p-tukey'])
    ES_SCondCoT[x] = float(Stats['eta-square'])
    if T_SCondCoT[x] < 0.05 :
        ax4.text(x+1-BarOfst, np.mean(CoT_Fp_S[:,x], axis=0)+1,
                 '*', c = Colors[x], fontsize=Ast, ha='center')
    
    d = np.reshape([CoT_Fp_F[:,x], CoT_Fp_F[:,2]], [40, 1]) 
    D = pd.DataFrame(np.hstack([d,G]), columns=['X','G'])
    Stats = pg.pairwise_tukey(D, dv='X', between='G', effsize='eta-square')
    T_FCondCoT[x] = float(Stats['p-tukey'])
    ES_FCondCoT[x] = float(Stats['eta-square'])
    if T_FCondCoT[x] < 0.05 :
        ax4.text(x+1+BarOfst, np.mean(CoT_Fp_F[:,x], axis=0)+1, 
                 '*', c = Colors[x], fontsize=Ast, ha='center')
        
print('CoT Speed Conditions Post Hoc')
print('p-values')
print(T_SCondCoT)
print('effect sizes')
print(np.round(ES_SCondCoT, decimals=5))
print('CoT Fp Conditions Post Hoc')
print('p-values')
print(T_FCondCoT)
print('effect sizes')
print(np.round(ES_FCondCoT, decimals=5))

# between clamps
T_CoT = np.ones(5)
ES_CoT = np.ones(5)
for x in range(5):
    d = np.reshape([CoT_Fp_S[:,x], CoT_Fp_F[:,x]], [40, 1])
    D = pd.DataFrame(np.hstack([d,G]), columns=['X','G'])
    Stats = pg.pairwise_tukey(D, dv='X', between='G', effsize='eta-square')
    T_CoT[x] = float(Stats['p-tukey'])
    ES_CoT[x] = float(Stats['eta-square'])
    if T_CoT[x] < 0.05 :
        y = np.mean([np.mean(CoT_Fp_S[:,x], axis=0), 
                     np.mean(CoT_Fp_F[:,x], axis=0)])
        ax4.text(x+1, y+1.1, '#', 
                 c = Colors[x], fontsize=Hash, ha='center')

print('CoT Between Clamp Post Hoc')
print('p-values')
print(np.round(T_CoT, decimals=5))
print('effect sizes')
print(np.round(ES_CoT, decimals=5))
print(' ')
print(' ')   

plt.savefig('Clamps.png', dpi=300)
plt.savefig('Clamps.pdf', dpi=300)

#%% Correlation Plot
plt.close('all')
fig = plt.figure(figsize=[18,12])
sz = 50
sz2 = 100
A = 0.4
fnt = 15
txt = 13

# speed by Fp
plt1 = plt.subplot(231)
for i in range(5):
    plt1.scatter(AllSpd_S[:,i], AllFp_S_kg[:,i]*9.81, 
                c=Colors[i], marker='o', s=sz)
    plt1.scatter(AllSpd_F[:,i], AllFp_F_kg[:,i]*9.81, 
                c=Colors[i], marker='s', s=sz)
    
# calculate trendlines
z_S = np.polyfit(np.hstack(AllSpd_S), np.hstack(AllFp_S_kg), 1)
p_S = np.poly1d(z_S)
x = np.linspace(np.min(AllSpd_S), np.max(AllSpd_S), 25)
plt1.scatter(x,p_S(x)*9.81,c='k',marker='o', s=sz2, alpha = A)
plt1.plot(x,p_S(x)*9.81,'-k')
# the line equation and R
# C_S = np.corrcoef(np.hstack(AllSpd_S), np.hstack(AllFp_S_kg))[0,1]
[c_S, P_S] = stats.pearsonr(np.hstack(AllSpd_S), np.hstack(AllFp_S_kg))
LineEq_S = 'y = ' + str(round(z_S[0],2)) + 'x + ' + str(round(z_S[1],2))
plt1.text(2.0, 19, 'Speed Clamp', fontsize=txt, ha='right')
plt1.text(2.0, 18, LineEq_S, fontsize=txt, ha='right')
plt1.text(2.0, 17, 'R$^2$ = ' + str(round(c_S*c_S,3)), fontsize=txt, ha='right')
plt1.text(2.0, 16, 'p < 0.001', fontsize=txt, ha='right')
plt1.text(0.85, 31, 'A', fontsize=fnt*2)
    
z_F = np.polyfit(np.hstack(AllSpd_F), np.hstack(AllFp_F_kg), 1)
p_F = np.poly1d(z_F)
x = np.linspace(np.min(AllSpd_F), np.max(AllSpd_F), 25)
plt1.scatter(x,p_F(x)*9.81,c='k',marker='s', s=sz2, alpha = A)
plt1.plot(x,p_F(x)*9.81,'-k')
# the line equation and R
# C_F = np.corrcoef(np.hstack(AllSpd_F), np.hstack(AllFp_F_kg))[0,1]
[c_F, P_F] = stats.pearsonr(np.hstack(AllSpd_F), np.hstack(AllFp_F_kg))
LineEq_F = 'y = ' + str(round(z_F[0],2)) + 'x + ' + str(round(z_F[1],2))
plt1.text(1.02, 29.5, 'Fp Clamp', fontsize=txt)
plt1.text(1.02, 28.5, LineEq_F, fontsize=txt)
plt1.text(1.02, 27.5, 'R$^2$ = ' + str(round(c_F*c_F,3)), fontsize=txt)
plt1.text(1.02, 26.5, 'p < 0.001', fontsize=txt)

plt1.set_xlabel('Walking Speed (m/s)', fontsize=fnt)
plt1.set_xticks([1, 1.5, 2])
plt1.set_ylabel('Fp (% body weight)', fontsize=fnt)
plt1.set_yticks([15, 20, 25, 30])
plt1.tick_params(axis='y', labelsize=fnt) 


# speed by net metabolic cost
plt2 = plt.subplot(232)
for i in range(5):
    plt2.scatter(AllSpd_S[:,i], WAvg_S[:,i], 
                c=Colors[i], marker='o', s=sz)
    plt2.scatter(AllSpd_F[:,i], WAvg_F[:,i], 
                c=Colors[i], marker='s', s=sz)
    
# calculate trendlines
z_S = np.polyfit(np.hstack(AllSpd_S), np.hstack(WAvg_S), 1)
p_S = np.poly1d(z_S)
x = np.linspace(np.min(AllSpd_S), np.max(AllSpd_S), 25)
plt2.scatter(x,p_S(x),c='k',marker='o', s=sz2, alpha = A)
plt2.plot(x,p_S(x),'-k')
# the line equation and R
[c_S, P_S] = stats.pearsonr(np.hstack(AllSpd_S), np.hstack(WAvg_S))
LineEq_S = 'y = ' + str(round(z_S[0],2)) + 'x + ' + str(round(z_S[1],2))
plt2.text(2.0, 3, 'Speed Clamp', fontsize=txt, ha='right')
plt2.text(2.0, 2.5, LineEq_S, fontsize=txt, ha='right')
plt2.text(2.0, 2, 'R$^2$ = ' + str(round(c_S*c_S,3)), fontsize=txt, ha='right')
plt2.text(2.0, 1.5, 'p < 0.001', fontsize=txt, ha='right')
plt2.text(0.85, 10, 'B', fontsize=fnt*2)
z_F = np.polyfit(np.hstack(AllSpd_F), np.hstack(WAvg_F), 1)
p_F = np.poly1d(z_F)
x = np.linspace(np.min(AllSpd_F), np.max(AllSpd_F), 25)
plt2.scatter(x,p_F(x),c='k',marker='s', s=sz2, alpha = A)
plt2.plot(x,p_F(x),'-k')
# the line equation and R
[c_F, P_F] = stats.pearsonr(np.hstack(AllSpd_F), np.hstack(WAvg_F))
LineEq_F = 'y = ' + str(round(z_F[0],2)) + 'x + ' + str(round(z_F[1],2))
plt2.text(1.02, 9, 'Fp Clamp', fontsize=txt)
plt2.text(1.02, 8.5, LineEq_F, fontsize=txt)
plt2.text(1.02, 8, 'R$^2$ = ' + str(round(c_F*c_F,3)), fontsize=txt)
plt2.text(1.02, 7.5, 'p < 0.001', fontsize=txt)
plt2.set_xlabel('Walking Speed (m/s)', fontsize=fnt)
plt2.set_xticks([1, 1.5, 2])
plt2.set_ylabel('Net Metabolic Power (W/kg)', fontsize=fnt)
plt2.set_yticks([2, 4, 6, 8])
plt2.tick_params(axis='y', labelsize=fnt) 


# speed by CoT
plt3 = plt.subplot(233)
for i in range(5):
    plt3.scatter(AllSpd_S[:,i], CoT_Fp_S[:,i], 
                c=Colors[i], marker='o', s=sz)
    plt3.scatter(AllSpd_F[:,i], CoT_Fp_F[:,i], 
                c=Colors[i], marker='s', s=sz)


z_S = np.polyfit(np.hstack(AllSpd_S), np.hstack(CoT_Fp_S), 2)
p_S = np.poly1d(z_S)
x = np.linspace(np.min(AllSpd_S), np.max(AllSpd_S), 25)
plt3.scatter(x,p_S(x),c='k',marker='o', s=sz2, alpha = A)
plt3.plot(x,p_S(x),'-k')
# the line equation and R
def polyfit(x, y, degree):
    results = {}
    coeffs = np.polyfit(x, y, degree)
    p = np.poly1d(coeffs)
    #calculate r-squared
    yhat = p(x)
    ybar = np.sum(y)/len(y)
    ssreg = np.sum((yhat-ybar)**2)
    sstot = np.sum((y - ybar)**2)
    results['r_squared'] = ssreg / sstot

    return results

R = polyfit(np.hstack(AllSpd_S), np.hstack(CoT_Fp_S), 2)
[c_S, P_S] = stats.pearsonr(np.hstack(AllSpd_S), np.hstack(CoT_Fp_S))
LineEq_S = 'y = ' + str(round(z_S[0],2)) + 'x$^2$ + ' + str(round(z_S[1],2)) + 'x + ' + str(round(z_S[2],2))
plt3.text(2.0, 2.25, 'Speed Clamp', fontsize=txt, ha='right')
plt3.text(2.0, 2, LineEq_S, fontsize=txt, ha='right')
plt3.text(2.0, 1.75, 'R$^2$ = ' + str(round(R['r_squared'], 3)), fontsize=txt, ha='right')
plt3.text(2.0, 1.5, 'p < 0.001', fontsize=txt, ha='right')
plt3.text(0.85, 6.2, 'C', fontsize=fnt*2)

z_F = np.polyfit(np.hstack(AllSpd_F), np.hstack(CoT_Fp_F), 2)
p_F = np.poly1d(z_F)
x = np.linspace(np.min(AllSpd_F), np.max(AllSpd_F), 25)
plt3.scatter(x,p_F(x),c='k',marker='s', s=sz2, alpha = A)
plt3.plot(x,p_F(x),'-k')
# the line equation and R
R = polyfit(np.hstack(AllSpd_F), np.hstack(CoT_Fp_F), 2)
[c_F, P_F] = stats.pearsonr(np.hstack(AllSpd_F), np.hstack(CoT_Fp_F))
LineEq_F = 'y = ' + str(round(z_F[0],2)) + 'x$^2$ + ' + str(round(z_F[1],2)) + 'x + ' + str(round(z_F[2],2))
plt3.text(1.02, 5.8, 'Fp Clamp', fontsize=txt)
plt3.text(1.02, 5.5, LineEq_F, fontsize=txt)
plt3.text(1.02, 5.2, 'R$^2$ = ' + str(round(R['r_squared'],3)), fontsize=txt)
plt3.text(1.02, 4.9, 'p < 0.001', fontsize=txt)

plt3.set_xlabel('Walking Speed (m/s)', fontsize=fnt)
plt3.set_xticks([1, 1.5, 2])
plt3.set_ylabel('Cost of Transport (J/kg/m)', fontsize=fnt)
plt3.set_yticks([2, 3, 4, 5, 6])
plt3.tick_params(axis='y', labelsize=fnt) 


# Fp by net metabolic cost
plt5 = plt.subplot(235)
for i in range(5):
    plt5.scatter(AllFp_S_kg[:,i]*9.81, WAvg_S[:,i], 
                c=Colors[i], marker='o', s=sz)
    plt5.scatter(AllFp_F_kg[:,i]*9.81, WAvg_F[:,i], 
                c=Colors[i], marker='s', s=sz)
    
# calculate trendlines
z_S = np.polyfit(np.hstack(AllFp_S_kg)*9.81, np.hstack(WAvg_S), 1)
p_S = np.poly1d(z_S)
x = np.linspace(np.min(AllFp_S_kg)*9.81, np.max(AllFp_S_kg)*9.81, 25)
plt5.scatter(x,p_S(x),c='k',marker='o', s=sz2, alpha = A)
plt5.plot(x,p_S(x),'-k')
# the line equation and R
# C_S = np.corrcoef(np.hstack(AllFp_S_kg)*9.81, np.hstack(WAvg_S))[0,1]
[c_S, P_S] = stats.pearsonr(np.hstack(AllFp_S_kg*9.81), np.hstack(WAvg_S))
LineEq_S = 'y = ' + str(round(z_S[0],2)) + 'x + ' + str(round(z_S[1],2))
plt5.text(30, 3, 'Speed Clamp', fontsize=txt, ha='right')
plt5.text(30, 2.5, LineEq_S, fontsize=txt, ha='right')
plt5.text(30, 2, 'R$^2$ = ' + str(round(c_S*c_S,3)), fontsize=txt, ha='right')
plt5.text(30, 1.5, 'p < 0.001', fontsize=txt, ha='right')
plt5.text(13.5, 9.7, 'D', fontsize=fnt*2)

z_F = np.polyfit(np.hstack(AllFp_F_kg)*9.81, np.hstack(WAvg_F), 1)
p_F = np.poly1d(z_F)
x = np.linspace(np.min(AllFp_F_kg)*9.81, np.max(AllFp_F_kg)*9.81, 25)
plt5.scatter(x,p_F(x),c='k',marker='s', s=sz2, alpha = A)
plt5.plot(x,p_F(x),'-k')
# the line equation and R
# C_F = np.corrcoef(np.hstack(AllFp_F_kg)*9.81, np.hstack(WAvg_F))[0,1]
[c_F, P_F] = stats.pearsonr(np.hstack(AllFp_F_kg*9.81), np.hstack(WAvg_F))
LineEq_F = 'y = ' + str(round(z_F[0],2)) + 'x + ' + str(round(z_F[1],2))
plt5.text(15.5, 9, 'Fp Clamp', fontsize=txt)
plt5.text(15.5, 8.5, LineEq_F, fontsize=txt)
plt5.text(15.5, 8, 'R$^2$ = ' + str(round(c_F*c_F,3)), fontsize=txt)
plt5.text(15.5, 7.5, 'p < 0.001', fontsize=txt)

plt5.set_xlabel('Fp (% body weight)', fontsize=fnt)
# plt5.set_xticks([1, 1.5, 2])
plt5.set_ylabel('Net Metabolic Power (W/kg)', fontsize=fnt)
plt5.set_yticks([2, 4, 6, 8])
plt5.tick_params(axis='y', labelsize=fnt) 


# Fp by CoT
plt6 = plt.subplot(236)
for i in range(5):
    plt6.scatter(AllFp_S_kg[:,i]*9.81, CoT_Fp_S[:,i], 
                c=Colors[i], marker='o', s=sz)
    plt6.scatter(AllFp_F_kg[:,i]*9.81, CoT_Fp_F[:,i], 
                c=Colors[i], marker='s', s=sz)
    
# calculate trendlines
z_S = np.polyfit(np.hstack(AllFp_S_kg)*9.81, np.hstack(CoT_Fp_S), 2)
p_S = np.poly1d(z_S)
x = np.linspace(np.min(AllFp_S_kg)*9.81, np.max(AllFp_S_kg)*9.81, 25)
plt6.scatter(x,p_S(x),c='k',marker='o', s=sz2, alpha = A)
plt6.plot(x,p_S(x),'-k')
# the line equation and R
R = polyfit(np.hstack(AllFp_S_kg)*9.81, np.hstack(CoT_Fp_S), 2)
[c_S, P_S] = stats.pearsonr(np.hstack(AllFp_S*9.81), np.hstack(CoT_Fp_S))
LineEq_S = 'y = ' + str(round(z_S[0],2)) + 'x$^2$ + ' + str(round(z_S[1],2)) + 'x + ' + str(round(z_S[2],2))
plt6.text(30, 2.3, 'Speed Clamp', fontsize=txt, ha='right')
plt6.text(30, 2, LineEq_S, fontsize=txt, ha='right')
plt6.text(30, 1.7, 'R$^2$ = ' + str(round(R['r_squared'],3)), fontsize=txt, ha='right')
plt6.text(30, 1.4, 'p < 0.001', fontsize=txt, ha='right')
plt6.text(13.5, 6.2, 'E', fontsize=fnt*2)

z_F = np.polyfit(np.hstack(AllFp_F_kg)*9.81, np.hstack(CoT_Fp_F), 2)
p_F = np.poly1d(z_F)
x = np.linspace(np.min(AllFp_F_kg)*9.81, np.max(AllFp_F_kg)*9.81, 25)
plt6.scatter(x,p_F(x),c='k',marker='s', s=sz2, alpha = A)
plt6.plot(x,p_F(x),'-k')
# the line equation and R
R = polyfit(np.hstack(AllFp_F_kg)*9.81, np.hstack(CoT_Fp_F), 2)
[c_F, P_F] = stats.pearsonr(np.hstack(AllFp_F*9.81), np.hstack(CoT_Fp_F))
# C_F = np.corrcoef(np.hstack(AllFp_F_kg)*9.81, np.hstack(CoT_Fp_F))[0,1]
LineEq_F = 'y = ' + str(round(z_F[0],2)) + 'x$^2$ + ' + str(round(z_F[1],2)) + 'x + ' + str(round(z_F[2],2))
plt6.text(15.5, 5.8, 'Fp Clamp', fontsize=txt)
plt6.text(15.5, 5.5, LineEq_F, fontsize=txt)
plt6.text(15.5, 5.2, 'R$^2$ = ' + str(round(R['r_squared'], 3)), fontsize=txt)
plt6.text(15.5, 4.9, 'p < 0.001', fontsize=txt)

plt6.set_xlabel('Fp (% body weight)', fontsize=fnt)
# plt6.set_xticks([1, 1.5, 2])
plt6.set_ylabel('Cost of Transport (J/kg/m)', fontsize=fnt)
plt6.set_yticks([2, 3, 4, 5, 6])
plt6.tick_params(axis='y', labelsize=fnt) 

# legend
plt4 = plt.subplot(234)
plt4.text(0.25, 0.6, 'Conditions', ha='center', weight='bold', style='italic', fontsize=18)
plt4.text(0.25, 0.55, '+20%', ha='center', c=Colors[4], weight='bold')
plt4.text(0.25, 0.5, '+10%', ha='center', c=Colors[3], weight='bold')
plt4.text(0.25, 0.45, 'Norm', ha='center', c=Colors[2], weight='bold')
plt4.text(0.25, 0.4, '-10%', ha='center', c=Colors[1], weight='bold')
plt4.text(0.25, 0.35, '-20%', ha='center', c=Colors[0], weight='bold')

plt4.plot(0.22, 0.85, 'ok', markersize=16)
plt4.plot(0.22, 0.75, 'sk', markersize=16)
plt4.text(0.23, 0.85, ' = Speed Clamp', va='center') 
plt4.text(0.23, 0.75, ' = Fp Clamp', va='center') 
plt4.set_xlim([0.15, 0.35])
plt4.set_ylim(0.3, 0.95)
plt4.axis('off')

plt.savefig('Corr.png', dpi=300)
plt.savefig('Corr.pdf', dpi=300)