%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate ECoG signal power for selected frequency bands
% Used for data collected during TDP surgeries (PT6 and PT7)

% Author: Jasmine (Zixuan) Zhou
% July, 2021

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize Variables and Load Signals for Monkey G

Fs = 1000; % downsampled frequency
channels = 1:32; % there are 32 channels

% X and Y locations of electrodes in array
x = [4 4 3 2 2 2 1 1 1 1 2 2 3 3 3 4 4 4 5 6 6 7 8 7 7 7 6 6 6 5 5 5]*0.75;
y = [3 1 2 1 5 3 4 2 8 6 9 7 6 4 8 5 7 9 8 9 7 8 5 6 2 4 3 5 1 6 4 2]*0.75;

load('C:\Users\jzhou33\Documents\MATLAB\PT7_NeurophysiologyData\DownsampledSignalsPT7.mat');
sigG_baseline = baselineDS;
sigG_PT = duringPTDS;
sigG_postPT = postPTDS;
sigG_stim_1 = duringStimDS_1;
sigG_stim_2 = duringStimDS_2;

disp('Loaded monkey G downsampled signals')

%% Prepare Baseline DownSampledSignals to Plot Power Spectral Density

% select freuqency range for the analysis
filt_rg = [4 7]; % theta
% filt_rg = [12 29]; % beta
% filt_rg = [30 59]; % low gamma
% filt_rg = [60 150]; % high gamma

% Design notch filters
d60 = designfilt('bandstopiir', 'FilterOrder', 2, 'HalfPowerFrequency1',...
    59.75, 'HalfPowerFrequency2', 60.25, 'DesignMethod', 'butter',...
    'SampleRate', Fs);
d120 = designfilt('bandstopiir', 'FilterOrder', 2, 'HalfPowerFrequency1',...
    119.75, 'HalfPowerFrequency2', 120.25, 'DesignMethod', 'butter',...
    'SampleRate', Fs);
d180 = designfilt('bandstopiir', 'FilterOrder', 2, 'HalfPowerFrequency1',...
    179.75, 'HalfPowerFrequency2', 180.25, 'DesignMethod', 'butter',...
    'SampleRate', Fs);
d240 = designfilt('bandstopiir', 'FilterOrder', 2, 'HalfPowerFrequency1',...
    239.75, 'HalfPowerFrequency2', 240.25, 'DesignMethod', 'butter',...
    'SampleRate', Fs);

% get signal for pre-PT baseline
for iChan = channels
    % 60Hz
    sigG_baseline.ipsi{iChan} = filtfilt(d60, sigG_baseline.ipsi{iChan});
    sigG_baseline.cont{iChan} = filtfilt(d60, sigG_baseline.cont{iChan});
    
    % 120 Hz
    sigG_baseline.ipsi{iChan} = filtfilt(d120, sigG_baseline.ipsi{iChan});
    sigG_baseline.cont{iChan} = filtfilt(d120, sigG_baseline.cont{iChan});
    
    % 180 Hz
    sigG_baseline.ipsi{iChan} = filtfilt(d180, sigG_baseline.ipsi{iChan});
    sigG_baseline.cont{iChan} = filtfilt(d180, sigG_baseline.cont{iChan});
    
    % 240 Hz
    sigG_baseline.ipsi{iChan} = filtfilt(d240, sigG_baseline.ipsi{iChan});
    sigG_baseline.cont{iChan} = filtfilt(d240, sigG_baseline.cont{iChan});
    
    sigG_baseline.ipsi{iChan} = bandpass(sigG_baseline.ipsi{iChan},filt_rg, Fs,...
        'ImpulseResponse', 'iir', 'Steepness', 0.95);
    sigG_baseline.cont{iChan} = bandpass(sigG_baseline.cont{iChan},filt_rg, Fs,...
        'ImpulseResponse', 'iir', 'Steepness', 0.95);
end
disp('Filtered Baseline Signal')

% get signal for during PT baseline
for iChan = channels
    % 60Hz
    sigG_PT.ipsi{iChan} = filtfilt(d60, sigG_PT.ipsi{iChan});
    sigG_PT.cont{iChan} = filtfilt(d60, sigG_PT.cont{iChan});
    
    % 120 Hz
    sigG_PT.ipsi{iChan} = filtfilt(d120, sigG_PT.ipsi{iChan});
    sigG_PT.cont{iChan} = filtfilt(d120, sigG_PT.cont{iChan});
    
    % 180 Hz
    sigG_PT.ipsi{iChan} = filtfilt(d180, sigG_PT.ipsi{iChan});
    sigG_PT.cont{iChan} = filtfilt(d180, sigG_PT.cont{iChan});
    
    % 240 Hz
    sigG_PT.ipsi{iChan} = filtfilt(d240, sigG_PT.ipsi{iChan});
    sigG_PT.cont{iChan} = filtfilt(d240, sigG_PT.cont{iChan});
    
    sigG_PT.ipsi{iChan} = bandpass(sigG_PT.ipsi{iChan},filt_rg, Fs,...
        'ImpulseResponse', 'iir', 'Steepness', 0.95);
    sigG_PT.cont{iChan} = bandpass(sigG_PT.cont{iChan},filt_rg, Fs,...
        'ImpulseResponse', 'iir', 'Steepness', 0.95);
end

% get signal for post PT
for iChan = channels
    % 60Hz
    sigG_postPT.ipsi{iChan} = filtfilt(d60, sigG_postPT.ipsi{iChan});
    sigG_postPT.cont{iChan} = filtfilt(d60, sigG_postPT.cont{iChan});
    
    % 120 Hz
    sigG_postPT.ipsi{iChan} = filtfilt(d120, sigG_postPT.ipsi{iChan});
    sigG_postPT.cont{iChan} = filtfilt(d120, sigG_postPT.cont{iChan});
    
    % 180 Hz
    sigG_postPT.ipsi{iChan} = filtfilt(d180, sigG_postPT.ipsi{iChan});
    sigG_postPT.cont{iChan} = filtfilt(d180, sigG_postPT.cont{iChan});
    
    % 240 Hz
    sigG_postPT.ipsi{iChan} = filtfilt(d240, sigG_postPT.ipsi{iChan});
    sigG_postPT.cont{iChan} = filtfilt(d240, sigG_postPT.cont{iChan});
    
    sigG_postPT.ipsi{iChan} = bandpass(sigG_postPT.ipsi{iChan},filt_rg, Fs,...
        'ImpulseResponse', 'iir', 'Steepness', 0.95);
    sigG_postPT.cont{iChan} = bandpass(sigG_postPT.cont{iChan},filt_rg, Fs,...
        'ImpulseResponse', 'iir', 'Steepness', 0.95);
end
disp('Filtered Post PT Signal')

% get signal during stimulation
for iChan = channels
    % 60Hz
    sigG_stim_1.ipsi{iChan} = filtfilt(d60, sigG_stim_1.ipsi{iChan});
    sigG_stim_1.cont{iChan} = filtfilt(d60, sigG_stim_1.cont{iChan});
    
    % 120 Hz
    sigG_stim_1.ipsi{iChan} = filtfilt(d120, sigG_stim_1.ipsi{iChan});
    sigG_stim_1.cont{iChan} = filtfilt(d120, sigG_stim_1.cont{iChan});
    
    % 180 Hz
    sigG_stim_1.ipsi{iChan} = filtfilt(d180, sigG_stim_1.ipsi{iChan});
    sigG_stim_1.cont{iChan} = filtfilt(d180, sigG_stim_1.cont{iChan});
    
    % 240 Hz
    sigG_stim_1.ipsi{iChan} = filtfilt(d240, sigG_stim_1.ipsi{iChan});
    sigG_stim_1.cont{iChan} = filtfilt(d240, sigG_stim_1.cont{iChan});
    
    sigG_stim_1.ipsi{iChan} = bandpass(sigG_stim_1.ipsi{iChan},filt_rg, Fs,...
        'ImpulseResponse', 'iir', 'Steepness', 0.95);
    sigG_stim_1.cont{iChan} = bandpass(sigG_stim_1.cont{iChan},filt_rg, Fs,...
        'ImpulseResponse', 'iir', 'Steepness', 0.95);
end
disp('Filtered Stim Signal')

% get signal for post stimulation baseline
for iChan = channels
    % 60Hz
    sigG_stim_2.ipsi{iChan} = filtfilt(d60, sigG_stim_2.ipsi{iChan});
    sigG_stim_2.cont{iChan} = filtfilt(d60, sigG_stim_2.cont{iChan});
    
    % 120 Hz
    sigG_stim_2.ipsi{iChan} = filtfilt(d120, sigG_stim_2.ipsi{iChan});
    sigG_stim_2.cont{iChan} = filtfilt(d120, sigG_stim_2.cont{iChan});
    
    % 180 Hz
    sigG_stim_2.ipsi{iChan} = filtfilt(d180, sigG_stim_2.ipsi{iChan});
    sigG_stim_2.cont{iChan} = filtfilt(d180, sigG_stim_2.cont{iChan});
    
    % 240 Hz
    sigG_stim_2.ipsi{iChan} = filtfilt(d240, sigG_stim_2.ipsi{iChan});
    sigG_stim_2.cont{iChan} = filtfilt(d240, sigG_stim_2.cont{iChan});
    
    sigG_stim_2.ipsi{iChan} = bandpass(sigG_stim_2.ipsi{iChan},filt_rg, Fs,...
        'ImpulseResponse', 'iir', 'Steepness', 0.95);
    sigG_stim_2.cont{iChan} = bandpass(sigG_stim_2.cont{iChan},filt_rg, Fs,...
        'ImpulseResponse', 'iir', 'Steepness', 0.95);
end
disp('Filtered Post-Stim Signal')


disp('Filtered All Signal')

%% Clean Signals by Removing Artifacts
% Samples with noise artifacts were found using FindArtifacts.m 
for iChan = channels
    
    cleanG_base.ipsi{iChan} = sigG_baseline.ipsi{iChan};
    cleanG_base.ipsi{iChan}(1:Fs) = nan; 
    
    cleanG_PT.ipsi{iChan} = sigG_PT.ipsi{iChan};
    cleanG_PT.ipsi{iChan}(1:Fs) = nan;
    cleanG_PT.ipsi{iChan}(131000:131000+Fs) = nan; 
    
    cleanG_postPT.ipsi{iChan} = sigG_postPT.ipsi{iChan};
    cleanG_postPT.ipsi{iChan}(1:Fs) = nan; 
    cleanG_postPT.ipsi{iChan}(1457000:1457000+Fs) = nan;
  
    cleanG_stim_1.ipsi{iChan} = sigG_stim_1.ipsi{iChan};
    
    cleanG_stim_2.ipsi{iChan} = sigG_stim_2.ipsi{iChan};
    cleanG_stim_2.ipsi{iChan}(2193000:2197000) = nan;
    cleanG_stim_2.ipsi{iChan}(2376000:2380000) = nan;
    cleanG_stim_2.ipsi{iChan}(2443000:2446000) = nan;

end
%% Check Cleaned Signals in Plots
figure;

% select channel to check
iChan = 12;

subplot(5,1,1)
t =  linspace(0,length(cleanG_base.ipsi{iChan})/60/Fs,length(cleanG_base.ipsi{iChan}));
plot(t,cleanG_base.ipsi{iChan})
title('baseline')
xlabel('time (min)')
ylim([-500 500])

subplot(5,1,2)
t =  linspace(0,length(cleanG_PT.ipsi{iChan})/60/Fs,length(cleanG_PT.ipsi{iChan}));
plot(t,cleanG_PT.ipsi{iChan})
title('during PT')
xlabel('time (min)')
ylim([-500 500])

subplot(5,1,3)
t =  linspace(0,length(cleanG_postPT.ipsi{iChan})/60/Fs,length(cleanG_PT.ipsi{iChan}));
plot(t,cleanG_postPT.ipsi{iChan})
title('post PT')
xlabel('time (min)')
ylim([-500 500])

subplot(5,1,4)
t =  linspace(0,length(cleanG_stim_1.ipsi{iChan})/60/Fs,length(cleanG_stim_1.ipsi{iChan}));
plot(t,cleanG_stim_1.ipsi{iChan})
title('during stim')
xlabel('time (min)')
ylim([-500 500])

subplot(5,1,5)
t =  linspace(0,length(cleanG_stim_2.ipsi{iChan})/60/Fs,length(cleanG_stim_2.ipsi{iChan}));
plot(t,cleanG_stim_2.ipsi{iChan})
title('after stim')
xlabel('time (min)')
ylim([-500 500])

%% Take Signal Power for All Recording Periods

% define the signal range for power analysis, in each recording period
prePT_baseline = 1 : 29*60*Fs;
duringPT_baseline = 1 : 25*60*Fs;
postPT_baseline_1 = 1 : 29*60*Fs;
postPT_baseline_2 = 30*60*Fs : 58*60*Fs;

stim_baseline_1 = 663*Fs : (663+115)*Fs;
stim_baseline_2 = 1383*Fs : (1383+115)*Fs;
stim_baseline_3 = 2103*Fs : (2103+115)*Fs;
stim_baseline_4 = 3000*Fs : (3000+115)*Fs;
stim_baseline_5 = 215*Fs : (215+115)*Fs;
postStim_baseline = 1020*Fs : 2730*Fs;

% initialize the variables

% For power averaged over each recording period
totPwrG_baseline = nan(1,length(channels));
totPwrG_duringPT = nan(1,length(channels));
totPwrG_postPT_1 = nan(1,length(channels));
totPwrG_postPT_2 = nan(1,length(channels));
totPwrG_stim_1 = nan(1,length(channels));
totPwrG_stim_2 = nan(1,length(channels));
totPwrG_stim_3 = nan(1,length(channels));
totPwrG_stim_4 = nan(1,length(channels));
totPwrG_stim_5 = nan(1,length(channels));
totPwrG_postStim = nan(1,length(channels));

% For power averaged every 10 seconds
G_base = nan(length(prePT_baseline),length(channels));
G_duringPT = nan(length(duringPT_baseline),length(channels));
G_postPT = nan(length(postPT_baseline_1(1):postPT_baseline_2(end)),length(channels));
G_stim = nan(3600*Fs,length(channels));
G_postStim = nan(length(postStim_baseline),length(channels));

% calculate signal power for each channel
for iChan = channels
    
    G_base (:,iChan) = cleanG_base.ipsi{iChan}(prePT_baseline);
    totPwrG_baseline(iChan) = SignalPower(cleanG_base.ipsi{iChan}(prePT_baseline),Fs);
    
    G_duringPT (:,iChan) = cleanG_PT.ipsi{iChan}(duringPT_baseline);
    totPwrG_duringPT(iChan) = SignalPower(cleanG_PT.ipsi{iChan}(duringPT_baseline),Fs);
    
    G_postPT (:,iChan) = cleanG_postPT.ipsi{iChan}(postPT_baseline_1(1):postPT_baseline_2(end));
    totPwrG_postPT_1(iChan) = SignalPower(cleanG_postPT.ipsi{iChan}(postPT_baseline_1),Fs);
    totPwrG_postPT_2(iChan) = SignalPower(cleanG_postPT.ipsi{iChan}(postPT_baseline_2),Fs);    
    
    G_stim (stim_baseline_1(1) : stim_baseline_1(1)+115*Fs,iChan) = cleanG_stim_1.ipsi{iChan}(stim_baseline_1);
    G_stim (stim_baseline_1(1)+12*60*Fs : stim_baseline_1(1)+115*Fs+12*60*Fs,iChan) = cleanG_stim_1.ipsi{iChan}(stim_baseline_2);
    G_stim (stim_baseline_1(1)+24*60*Fs : stim_baseline_1(1)+115*Fs+24*60*Fs,iChan) = cleanG_stim_1.ipsi{iChan}(stim_baseline_3);   
    G_stim (stim_baseline_1(1)+48*60*Fs : stim_baseline_1(1)+115*Fs+48*60*Fs,iChan) = cleanG_stim_2.ipsi{iChan}(stim_baseline_5);
    totPwrG_stim_1(iChan) = SignalPower(sigG_stim_1.ipsi{iChan}(stim_baseline_1),Fs);
    totPwrG_stim_2(iChan) = SignalPower(sigG_stim_1.ipsi{iChan}(stim_baseline_2),Fs);
    totPwrG_stim_3(iChan) = SignalPower(sigG_stim_1.ipsi{iChan}(stim_baseline_3),Fs);
    totPwrG_stim_5(iChan) = SignalPower(sigG_stim_2.ipsi{iChan}(stim_baseline_5),Fs);
    
    G_postStim (:,iChan) = cleanG_stim_2.ipsi{iChan}(postStim_baseline);
    totPwrG_postStim(iChan) = SignalPower(cleanG_stim_2.ipsi{iChan}(postStim_baseline),Fs);
end



%% Optional: Change Variable Names to Save Data

totPwrG_baseline_theta = totPwrG_baseline;
totPwrG_postPT_1_theta = totPwrG_postPT_1;
totPwrG_postPT_2_theta = totPwrG_postPT_2;
totPwrG_stim_1_theta = totPwrG_stim_1;
totPwrG_stim_2_theta = totPwrG_stim_2;
totPwrG_stim_3_theta = totPwrG_stim_3;
totPwrG_stim_5_theta = totPwrG_stim_5;
totPwrG_postStim_theta = totPwrG_postStim;

% totPwrG_baseline_beta = totPwrG_baseline;
% totPwrG_postPT_1_beta = totPwrG_postPT_1;
% totPwrG_postPT_2_beta = totPwrG_postPT_2;
% totPwrG_stim_1_beta = totPwrG_stim_1;
% totPwrG_stim_2_beta = totPwrG_stim_2;
% totPwrG_stim_3_beta = totPwrG_stim_3;
% totPwrG_stim_5_beta = totPwrG_stim_5;
% totPwrG_postStim_beta = totPwrG_postStim;

% totPwrG_baseline_lg = totPwrG_baseline;
% totPwrG_postPT_1_lg = totPwrG_postPT_1;
% totPwrG_postPT_2_lg = totPwrG_postPT_2;
% totPwrG_stim_1_lg = totPwrG_stim_1;
% totPwrG_stim_2_lg = totPwrG_stim_2;
% totPwrG_stim_3_lg = totPwrG_stim_3;
% totPwrG_stim_5_lg = totPwrG_stim_5;
% totPwrG_postStim_lg = totPwrG_postStim;

% totPwrG_baseline_hg = totPwrG_baseline;
% totPwrG_postPT_1_hg = totPwrG_postPT_1;
% totPwrG_postPT_2_hg = totPwrG_postPT_2;
% totPwrG_stim_1_hg = totPwrG_stim_1;
% totPwrG_stim_2_hg = totPwrG_stim_2;
% totPwrG_stim_3_hg = totPwrG_stim_3;
% totPwrG_stim_5_hg = totPwrG_stim_5;
% totPwrG_postStim_hg = totPwrG_postStim;
%}

%% Optional: Plot Heat Maps of Power for Visualization (Monkey G)
%{
% define color range
minF_hg = min([totPwrG_stim_1, totPwrG_stim_2,totPwrG_stim_3,...
    totPwrG_stim_5, totPwrG_postPT_1, totPwrG_postPT_2, ...
    totPwrG_baseline, totPwrG_postStim]);
maxF_hg = max([totPwrG_stim_1, totPwrG_stim_2,totPwrG_stim_3,...
    totPwrG_stim_5, totPwrG_postPT_1, totPwrG_postPT_2, ...
    totPwrG_baseline, totPwrG_postStim]);

% define color map
load('pink map.mat')
map = map_new;

figure, set(gcf, 'color', 'w', 'units', 'normalize', 'outerposition', [.3 0 .5 0.5])
subplot(2,4,1), mapArray_new(totPwrG_baseline(channels),x(channels),y(channels),[],[],map), axis off
title('Before stimulation', 'fontweight', 'bold')
caxis([minF_hg maxF_hg]);
colorbar(gca,'WestOutside')
set(gca,'fontsize',10)

subplot(2,4,2), mapArray_new(totPwrG_postPT_1(channels),x(channels),y(channels),[],[],map), axis off
title('Post-stroke, pre-stim', 'fontweight', 'bold')
caxis([minF_hg maxF_hg]);
colorbar(gca,'WestOutside')
set(gca,'fontsize',10)

subplot(2,4,3), mapArray_new(totPwrG_stim_1(channels),x(channels),y(channels),[],[],map), axis off
title('After 10min of stim', 'fontweight', 'bold')
caxis([minF_hg maxF_hg]);
colorbar(gca,'WestOutside')
set(gca,'fontsize',10)

subplot(2,4,4), mapArray_new(totPwrG_stim_2(channels),x(channels),y(channels),[],[],map), axis off
title('After 20min of stim', 'fontweight', 'bold')
caxis([minF_hg maxF_hg]);
colorbar(gca,'WestOutside')
set(gca,'fontsize',10)

subplot(2,4,5), mapArray_new(totPwrG_stim_3(channels),x(channels),y(channels),[],[],map), axis off
title('After 30min of stim', 'fontweight', 'bold')
caxis([minF_hg maxF_hg]);
colorbar(gca,'WestOutside')
set(gca,'fontsize',10)

subplot(2,4,6), mapArray_new(totPwrG_stim_5(channels),x(channels),y(channels),[],[],map), axis off
title('After 50min of stim', 'fontweight', 'bold')
caxis([minF_hg maxF_hg]);
colorbar(gca,'WestOutside')
set(gca,'fontsize',10)

subplot(2,4,7), mapArray_new(totPwrG_postStim(channels),x(channels),y(channels),[],[],map), axis off
title('After 60min of stim', 'fontweight', 'bold')
caxis([minF_hg maxF_hg]);
colorbar(gca,'WestOutside')
set(gca,'fontsize',10)

sgtitle('Ipsilesional power')

%}

%% Calculate Time-Series Power (dt = 10s) for Monkey G

dt = 10*Fs; % number of samples every 10 seconds

% initialize variables
lngth_allbase = size(G_base,1);
lngth_allduringPT = size(G_duringPT,1);
lngth_allpostPT = size(G_postPT,1);
lngth_allstim = size(G_stim,1);
lngth_allpostStim = size(G_postStim,1);

% number of 10 second intervals
Ndt_allbase = floor(lngth_allbase/dt); 
Ndt_allduringPT = floor(lngth_allduringPT/dt);
Ndt_allpostPT = floor(lngth_allpostPT/dt);
Ndt_allstim = floor(lngth_allstim/dt);
Ndt_allpostStim = floor(lngth_allpostStim/dt);

pwrG_allbase = nan(Ndt_allbase,length(channels));
pwrG_allduringPT = nan(Ndt_allduringPT,length(channels));
pwrG_allpostPT = nan(Ndt_allpostPT,length(channels));
pwrG_allstim = nan(Ndt_allstim,length(channels));
pwrG_allpostStim = nan(Ndt_allpostStim,length(channels));

for idt = 1:Ndt_allbase   
    interval = ((idt-1)*dt+1):(idt*dt); % time interval in samples
    for iChan = channels
        pwrG_allbase(idt,iChan) = SignalPower(G_base(interval,iChan),Fs);
    end
end

for idt = 1:Ndt_allduringPT   
    interval = ((idt-1)*dt+1):(idt*dt); % time interval in samples
    for iChan = channels
        pwrG_allduringPT(idt,iChan) = SignalPower(G_duringPT(interval,iChan),Fs);
    end
end


for idt = 1:Ndt_allpostPT
    interval = ((idt-1)*dt+1):(idt*dt); % time interval in samples   
    for iChan = channels
        pwrG_allpostPT(idt,iChan) = SignalPower(G_postPT(interval,iChan),Fs);
    end
end

for idt = 1:Ndt_allstim   
    interval = ((idt-1)*dt+1):(idt*dt); % time interval in samples
    for iChan = channels
        pwrG_allstim(idt,iChan) = SignalPower(G_stim(interval,iChan),Fs);       
    end
end

for idt = 1:Ndt_allpostStim  
    interval = ((idt-1)*dt+1):(idt*dt); % time interval in samples
    for iChan = channels
        pwrG_allpostStim(idt,iChan) = SignalPower(G_postStim(interval,iChan),Fs);      
    end
end

%% Optional: Change Variable Names to Save Data (Monkey G)

pwrG_allbase_th = pwrG_allbase;
pwrG_allduringPT_th = pwrG_allduringPT;
pwrG_allpostPT_th = pwrG_allpostPT;
pwrG_allstim_th = pwrG_allstim;
pwrG_allpostStim_th = pwrG_allpostStim;

%% Optional: Plot Time-Series Power (dt = 10s) for Monkey G 

%{
figure;
subplot(5,1,1)
plot(pwrG_allbase)

subplot(5,1,2)
plot(pwrG_allduringPT)

subplot(5,1,3)
plot(pwrG_allpostPT)

subplot(5,1,4)
plot(pwrG_allstim)

subplot(5,1,5)
plot(pwrG_allpostStim)
%}


%% Initialize Variables and Load Signals for Monkey F

good_chans_ipsi = [3:4 6 11:16 18:20 22 25 26:27 29:31];

% X and Y locations of electrodes in array
x = [4 4 3 2 2 2 1 1 1 1 2 2 3 3 3 4 4 4 5 6 6 7 8 7 7 7 6 6 6 5 5 5]*0.75;
y = [3 1 2 1 5 3 4 2 8 6 9 7 6 4 8 5 7 9 8 9 7 8 5 6 2 4 3 5 1 6 4 2]*0.75;

% load downsampled signals
load('C:\Users\jzhou33\Documents\MATLAB\PT6_NeurophysiologyData\DownsampledSignalsPT6.mat');
sigF_baseline = baselineDS;
sigF_PT = duringPTDS;
sigF_postPT = postPTDS;
sigF_stim = duringStimDS;
sigF_poststim = postStimDS;

disp('Loaded monkey F downsampled signals')

%% Prepare Baseline DownSampledSignals to Plot Power Spectral Density

% select freuqency range for the analysis

% filt_rg = [4 7]; % theta
% filt_rg = [12 29]; % beta
% filt_rg = [30 59]; % low gamma
filt_rg = [60 150]; % high gamma

% Design notch filters
d60 = designfilt('bandstopiir', 'FilterOrder', 2, 'HalfPowerFrequency1',...
    59.75, 'HalfPowerFrequency2', 60.25, 'DesignMethod', 'butter',...
    'SampleRate', Fs);
d120 = designfilt('bandstopiir', 'FilterOrder', 2, 'HalfPowerFrequency1',...
    119.75, 'HalfPowerFrequency2', 120.25, 'DesignMethod', 'butter',...
    'SampleRate', Fs);
d180 = designfilt('bandstopiir', 'FilterOrder', 2, 'HalfPowerFrequency1',...
    179.75, 'HalfPowerFrequency2', 180.25, 'DesignMethod', 'butter',...
    'SampleRate', Fs);
d240 = designfilt('bandstopiir', 'FilterOrder', 2, 'HalfPowerFrequency1',...
    239.75, 'HalfPowerFrequency2', 240.25, 'DesignMethod', 'butter',...
    'SampleRate', Fs);

% Filter signals for pre PT baseline
for iChan = channels
    % 60Hz
    sigF_baseline.ipsi{iChan} = filtfilt(d60, sigF_baseline.ipsi{iChan});
    sigF_baseline.cont{iChan} = filtfilt(d60, sigF_baseline.cont{iChan});
    
    % 120 Hz
    sigF_baseline.ipsi{iChan} = filtfilt(d120, sigF_baseline.ipsi{iChan});
    sigF_baseline.cont{iChan} = filtfilt(d120, sigF_baseline.cont{iChan});
    
    % 180 Hz
    sigF_baseline.ipsi{iChan} = filtfilt(d180, sigF_baseline.ipsi{iChan});
    sigF_baseline.cont{iChan} = filtfilt(d180, sigF_baseline.cont{iChan});
    
    % 240 Hz
    sigF_baseline.ipsi{iChan} = filtfilt(d240, sigF_baseline.ipsi{iChan});
    sigF_baseline.cont{iChan} = filtfilt(d240, sigF_baseline.cont{iChan});
    
    sigF_baseline.ipsi{iChan} = bandpass(sigF_baseline.ipsi{iChan},filt_rg, Fs,...
        'ImpulseResponse', 'iir', 'Steepness', 0.95);
    sigF_baseline.cont{iChan} = bandpass(sigF_baseline.cont{iChan},filt_rg, Fs,...
        'ImpulseResponse', 'iir', 'Steepness', 0.95);
end

% Filter signals for during PT
for iChan = channels
    % 60Hz
    sigF_PT.ipsi{iChan} = filtfilt(d60, sigF_PT.ipsi{iChan});
    sigF_PT.cont{iChan} = filtfilt(d60, sigF_PT.cont{iChan});
    
    % 120 Hz
    sigF_PT.ipsi{iChan} = filtfilt(d120, sigF_PT.ipsi{iChan});
    sigF_PT.cont{iChan} = filtfilt(d120, sigF_PT.cont{iChan});
    
    % 180 Hz
    sigF_PT.ipsi{iChan} = filtfilt(d180, sigF_PT.ipsi{iChan});
    sigF_PT.cont{iChan} = filtfilt(d180, sigF_PT.cont{iChan});
    
    % 240 Hz
    sigF_PT.ipsi{iChan} = filtfilt(d240, sigF_PT.ipsi{iChan});
    sigF_PT.cont{iChan} = filtfilt(d240, sigF_PT.cont{iChan});
    
    sigF_PT.ipsi{iChan} = bandpass(sigF_PT.ipsi{iChan},filt_rg, Fs,...
        'ImpulseResponse', 'iir', 'Steepness', 0.95);
    sigF_PT.cont{iChan} = bandpass(sigF_PT.cont{iChan},filt_rg, Fs,...
        'ImpulseResponse', 'iir', 'Steepness', 0.95);
end

% Filter signals for post PT
for iChan = channels
    % 60Hz
    sigF_postPT.ipsi{iChan} = filtfilt(d60, sigF_postPT.ipsi{iChan});
    sigF_postPT.cont{iChan} = filtfilt(d60, sigF_postPT.cont{iChan});
    
    % 120 Hz
    sigF_postPT.ipsi{iChan} = filtfilt(d120, sigF_postPT.ipsi{iChan});
    sigF_postPT.cont{iChan} = filtfilt(d120, sigF_postPT.cont{iChan});
    
    % 180 Hz
    sigF_postPT.ipsi{iChan} = filtfilt(d180, sigF_postPT.ipsi{iChan});
    sigF_postPT.cont{iChan} = filtfilt(d180, sigF_postPT.cont{iChan});
    
    % 240 Hz
    sigF_postPT.ipsi{iChan} = filtfilt(d240, sigF_postPT.ipsi{iChan});
    sigF_postPT.cont{iChan} = filtfilt(d240, sigF_postPT.cont{iChan});
    
    sigF_postPT.ipsi{iChan} = bandpass(sigF_postPT.ipsi{iChan},filt_rg, Fs,...
        'ImpulseResponse', 'iir', 'Steepness', 0.95);
    sigF_postPT.cont{iChan} = bandpass(sigF_postPT.cont{iChan},filt_rg, Fs,...
        'ImpulseResponse', 'iir', 'Steepness', 0.95);
end

% Filter signals for during stimulation
for iChan = channels
    % 60Hz
    sigFStim.ipsi{iChan} = filtfilt(d60, sigF_stim.ipsi{iChan});
    sigFStim.cont{iChan} = filtfilt(d60, sigF_stim.cont{iChan});
    
    % 120 Hz
    sigFStim.ipsi{iChan} = filtfilt(d120, sigFStim.ipsi{iChan});
    sigFStim.cont{iChan} = filtfilt(d120, sigFStim.cont{iChan});
    
    % 180 Hz
    sigFStim.ipsi{iChan} = filtfilt(d180, sigFStim.ipsi{iChan});
    sigFStim.cont{iChan} = filtfilt(d180, sigFStim.cont{iChan});
    
    % 240 Hz
    sigFStim.ipsi{iChan} = filtfilt(d240, sigFStim.ipsi{iChan});
    sigFStim.cont{iChan} = filtfilt(d240, sigFStim.cont{iChan});
    
    sigFStim.ipsi{iChan} = bandpass(sigFStim.ipsi{iChan},filt_rg, Fs,...
        'ImpulseResponse', 'iir', 'Steepness', 0.95);
    sigFStim.cont{iChan} = bandpass(sigFStim.cont{iChan},filt_rg, Fs,...
        'ImpulseResponse', 'iir', 'Steepness', 0.95);
end

% Filter signals for post stimulation baseline
for iChan = channels
    % 60Hz
    sigF_poststim.ipsi{iChan} = filtfilt(d60, sigF_poststim.ipsi{iChan});
    sigF_poststim.cont{iChan} = filtfilt(d60, sigF_poststim.cont{iChan});
    
    % 120 Hz
    sigF_poststim.ipsi{iChan} = filtfilt(d120, sigF_poststim.ipsi{iChan});
    sigF_poststim.cont{iChan} = filtfilt(d120, sigF_poststim.cont{iChan});
    
    % 180 Hz
    sigF_poststim.ipsi{iChan} = filtfilt(d180, sigF_poststim.ipsi{iChan});
    sigF_poststim.cont{iChan} = filtfilt(d180, sigF_poststim.cont{iChan});
    
    % 240 Hz
    sigF_poststim.ipsi{iChan} = filtfilt(d240, sigF_poststim.ipsi{iChan});
    sigF_poststim.cont{iChan} = filtfilt(d240, sigF_poststim.cont{iChan});
    
    sigF_poststim.ipsi{iChan} = bandpass(sigF_poststim.ipsi{iChan},filt_rg, Fs,...
        'ImpulseResponse', 'iir', 'Steepness', 0.95);
    sigF_poststim.cont{iChan} = bandpass(sigF_poststim.cont{iChan},filt_rg, Fs,...
        'ImpulseResponse', 'iir', 'Steepness', 0.95);
end

disp('Filtered All Signal')

%% Clean Signals by Removing Artifacts
% Samples with noise artifacts were found using FindArtifacts.m 

for iChan = good_chans_ipsi
    
    % baseline period
    cleanF_base.ipsi{iChan} = sigF_baseline.ipsi{iChan};
    cleanF_base.ipsi{iChan}(1:Fs) = nan;    
    
    % PT illumination period
    cleanF_PT.ipsi{iChan} = sigF_PT.ipsi{iChan};
    cleanF_PT.ipsi{iChan}(1:2*Fs) = nan;
    cleanF_PT.ipsi{iChan}(415000:440000+2*Fs) = nan;
    if iChan == 18 || iChan == 26 || iChan == 27
       cleanF_PT.ipsi{iChan}(268000:272000) = nan;   
    elseif iChan == 25
       cleanF_PT.ipsi{iChan}(1:273000) = nan;
    end
    
    % post PT period
    cleanF_postPT.ipsi{iChan} = sigF_postPT.ipsi{iChan};
    cleanF_postPT.ipsi{iChan}(1:Fs) = nan; 
    cleanF_postPT.ipsi{iChan}(403000:403000+Fs) = nan;
    if iChan == 25
       cleanF_postPT.ipsi{iChan}(3099000:3099000+Fs) = nan;
    end
    
    % post stimulation period
    cleanF_poststim.ipsi{iChan} = sigF_poststim.ipsi{iChan};
    cleanF_poststim.ipsi{iChan}(1:30*Fs) = nan;
    cleanF_poststim.ipsi{iChan}(1060000:1060000+2*Fs) = nan;
    if iChan == 6
        cleanF_poststim.ipsi{iChan}(2.6*60*Fs:3.9*60*Fs) = nan;
    end
    
end
%% Check Cleaned Signals in Plots 
figure;

% select channel to check
iChan = 25;

subplot(4,1,1)
t =  linspace(0,length(cleanF_base.ipsi{iChan})/60/Fs,length(cleanF_base.ipsi{iChan}));
plot(t,cleanF_base.ipsi{iChan})
title('pre-PT baseline')
xlabel('time (min)')
ylim([-500 500])

subplot(4,1,2)
t =  linspace(0,length(cleanF_PT.ipsi{iChan})/60/Fs,length(cleanF_PT.ipsi{iChan}));
plot(t,cleanF_PT.ipsi{iChan})
title('during PT')
xlabel('time (min)')
ylim([-500 500])

subplot(4,1,3)
t =  linspace(0,length(cleanF_postPT.ipsi{iChan})/60/Fs,length(cleanF_postPT.ipsi{iChan}));
plot(t,cleanF_postPT.ipsi{iChan})
title('post PT')
xlabel('time (min)')
ylim([-500 500])

subplot(4,1,4)
t =  linspace(0,length(cleanF_poststim.ipsi{iChan})/60/Fs,length(cleanF_poststim.ipsi{iChan}));
plot(t,cleanF_poststim.ipsi{iChan})
title('post stimulation')
xlabel('time (min)')
ylim([-500 500])
%% Take Signal Power for All Recording Periods (Monkey F)

% define the signal range for power analysis, in each recording period
prePT_baseline = 1 : 29*60*Fs;
duringPT_baseline = 1 : 29*60*Fs;
postPT_baseline_1 = 1 : 29*60*Fs;
postPT_baseline_2 = 30*60*Fs : 58*60*Fs;

stim_baseline_1 = 615*Fs : (615+115)*Fs;
stim_baseline_2 = 1335*Fs : (1335+115)*Fs;
stim_baseline_3 = 2055*Fs : (2055+115)*Fs;
stim_baseline_4 = 2775*Fs : (2775+115)*Fs;
stim_baseline_5 = 3495*Fs : (3495+115)*Fs;
postStim_baseline = 1 : 29*60*Fs;

totPwrF_baseline = nan(1,length(channels));
totPwrF_duringPT = nan(1,length(channels));
totPwrF_postPT_1 = nan(1,length(channels));
totPwrF_postPT_2 = nan(1,length(channels));
totPwrF_stim_1 = nan(1,length(channels));
totPwrF_stim_2 = nan(1,length(channels));
totPwrF_stim_3 = nan(1,length(channels));
totPwrF_stim_4 = nan(1,length(channels));
totPwrF_stim_5 = nan(1,length(channels));
totPwrF_postStim = nan(1,length(channels));

F_base = nan(length(prePT_baseline),length(channels));
F_duringPT = nan(length(duringPT_baseline),length(channels));
F_postPT = nan(length(postPT_baseline_1(1):postPT_baseline_2(end)),length(channels));
F_stim = nan(length(1:stim_baseline_5(end)),length(channels));
F_postStim = nan(length(postStim_baseline),length(channels));

% calculate signal power for each channel
for iChan = good_chans_ipsi
    
    F_base (:,iChan) = cleanF_base.ipsi{iChan}(prePT_baseline);
    totPwrF_baseline(iChan) = SignalPower(cleanF_base.ipsi{iChan}(prePT_baseline),Fs);
    
    F_duringPT (:,iChan) = cleanF_PT.ipsi{iChan}(duringPT_baseline);
    totPwrF_duringPT(iChan) = SignalPower(cleanF_PT.ipsi{iChan}(duringPT_baseline),Fs);
    
    F_postPT (:,iChan) = cleanF_postPT.ipsi{iChan}(postPT_baseline_1(1):postPT_baseline_2(end));
    totPwrF_postPT_1(iChan) = SignalPower(cleanF_postPT.ipsi{iChan}(postPT_baseline_1),Fs);
    totPwrF_postPT_2(iChan) = SignalPower(cleanF_postPT.ipsi{iChan}(postPT_baseline_2),Fs);
    
    F_stim (stim_baseline_1(1) : stim_baseline_1(1)+115*Fs,iChan) = sigFStim.ipsi{iChan}(stim_baseline_1);
    F_stim (stim_baseline_1(1)+12*60*Fs : stim_baseline_1(1)+115*Fs+12*60*Fs,iChan) = sigFStim.ipsi{iChan}(stim_baseline_2);
    F_stim (stim_baseline_1(1)+24*60*Fs : stim_baseline_1(1)+115*Fs+24*60*Fs,iChan) = sigFStim.ipsi{iChan}(stim_baseline_3);
    F_stim (stim_baseline_1(1)+36*60*Fs : stim_baseline_1(1)+115*Fs+36*60*Fs,iChan) = sigFStim.ipsi{iChan}(stim_baseline_4);
    F_stim (stim_baseline_1(1)+48*60*Fs : stim_baseline_1(1)+115*Fs+48*60*Fs,iChan) = sigFStim.ipsi{iChan}(stim_baseline_5);
    totPwrF_stim_1(iChan) = SignalPower(sigFStim.ipsi{iChan}(stim_baseline_1),Fs);
    totPwrF_stim_2(iChan) = SignalPower(sigFStim.ipsi{iChan}(stim_baseline_2),Fs);
    totPwrF_stim_3(iChan) = SignalPower(sigFStim.ipsi{iChan}(stim_baseline_3),Fs);
    totPwrF_stim_4(iChan) = SignalPower(sigFStim.ipsi{iChan}(stim_baseline_4),Fs);
    totPwrF_stim_5(iChan) = SignalPower(sigFStim.ipsi{iChan}(stim_baseline_5),Fs);
    
    F_postStim (:,iChan) = cleanF_poststim.ipsi{iChan}(postStim_baseline);
    totPwrF_postStim(iChan) = SignalPower(cleanF_poststim.ipsi{iChan}(postStim_baseline),Fs);
end


%% Optional: Change Variable Names to Save Data (Monkey F)

%
totPwrF_baseline_theta = totPwrF_baseline;
totPwrF_postPT_1_theta = totPwrF_postPT_1;
totPwrF_postPT_2_theta = totPwrF_postPT_2;
totPwrF_stim_1_theta = totPwrF_stim_1;
totPwrF_stim_2_theta = totPwrF_stim_2;
totPwrF_stim_3_theta = totPwrF_stim_3;
totPwrF_stim_4_theta = totPwrF_stim_4;
totPwrF_stim_5_theta = totPwrF_stim_5;
totPwrF_postStim_theta = totPwrF_postStim;

% totPwrF_baseline_beta = totPwrF_baseline;
% totPwrF_postPT_1_beta = totPwrF_postPT_1;
% totPwrF_postPT_2_beta = totPwrF_postPT_2;
% totPwrF_stim_1_beta = totPwrF_stim_1;
% totPwrF_stim_2_beta = totPwrF_stim_2;
% totPwrF_stim_3_beta = totPwrF_stim_3;
% totPwrF_stim_4_beta = totPwrF_stim_4;
% totPwrF_stim_5_beta = totPwrF_stim_5;
% totPwrF_postStim_beta = totPwrF_postStim;

% totPwrF_baseline_lg = totPwrF_baseline;
% totPwrF_postPT_1_lg = totPwrF_postPT_1;
% totPwrF_postPT_2_lg = totPwrF_postPT_2;
% totPwrF_stim_1_lg = totPwrF_stim_1;
% totPwrF_stim_2_lg = totPwrF_stim_2;
% totPwrF_stim_3_lg = totPwrF_stim_3;
% totPwrF_stim_4_lg = totPwrF_stim_4;
% totPwrF_stim_5_lg = totPwrF_stim_5;
% totPwrF_postStim_lg = totPwrF_postStim;

% totPwrF_baseline_hg = totPwrF_baseline;
% totPwrF_postPT_1_hg = totPwrF_postPT_1;
% totPwrF_postPT_2_hg = totPwrF_postPT_2;
% totPwrF_stim_1_hg = totPwrF_stim_1;
% totPwrF_stim_2_hg = totPwrF_stim_2;
% totPwrF_stim_3_hg = totPwrF_stim_3;
% totPwrF_stim_4_hg = totPwrF_stim_4;
% totPwrF_stim_5_hg = totPwrF_stim_5;
% totPwrF_postStim_hg = totPwrF_postStim;
%}

%% Optional: Plot Heat Maps of Power for Visualization (Monkey F)
%{
% define color range
minF_hg = min([totPwrF_stim_1, totPwrF_stim_2,totPwrF_stim_3,...
    totPwrF_stim_4, totPwrF_stim_5, totPwrF_postPT_1, totPwrF_postPT_2,...
    totPwrF_baseline, totPwrF_postStim]);
maxF_hg = max([totPwrF_stim_1, totPwrF_stim_2,totPwrF_stim_3,...
    totPwrF_stim_4, totPwrF_stim_5, totPwrF_postPT_1, totPwrF_postPT_2,...
    totPwrF_baseline, totPwrF_postStim]);

% define color map
load('pink map.mat')
map = map_new;

figure, set(gcf, 'color', 'w', 'units', 'normalize', 'outerposition', [.3 0 .5 0.5])
subplot(2,4,1), mapArray_new(totPwrF_baseline(good_chans_ipsi),x(good_chans_ipsi),y(good_chans_ipsi),[],[],map), axis off
title('Before stimulation', 'fontweight', 'bold')
caxis([minF_hg maxF_hg]);
colorbar(gca,'WestOutside')
set(gca,'fontsize',10)

subplot(2,4,2), mapArray_new(totPwrF_postPT_1(good_chans_ipsi),x(good_chans_ipsi),y(good_chans_ipsi),[],[],map), axis off
title('Post-stroke, pre-stim', 'fontweight', 'bold')
caxis([minF_hg maxF_hg]);
colorbar(gca,'WestOutside')
set(gca,'fontsize',10)

subplot(2,4,3), mapArray_new(totPwrF_stim_1(good_chans_ipsi),x(good_chans_ipsi),y(good_chans_ipsi),[],[],map), axis off
title('After 10min of stim', 'fontweight', 'bold')
caxis([minF_hg maxF_hg]);
colorbar(gca,'WestOutside')
set(gca,'fontsize',10)

subplot(2,4,4), mapArray_new(totPwrF_stim_2(good_chans_ipsi),x(good_chans_ipsi),y(good_chans_ipsi),[],[],map), axis off
title('After 20min of stim', 'fontweight', 'bold')
caxis([minF_hg maxF_hg]);
colorbar(gca,'WestOutside')
set(gca,'fontsize',10)

subplot(2,4,5), mapArray_new(totPwrF_stim_3(good_chans_ipsi),x(good_chans_ipsi),y(good_chans_ipsi),[],[],map), axis off
title('After 30min of stim', 'fontweight', 'bold')
caxis([minF_hg maxF_hg]);
colorbar(gca,'WestOutside')
set(gca,'fontsize',10)

subplot(2,4,6), mapArray_new(totPwrF_stim_4(good_chans_ipsi),x(good_chans_ipsi),y(good_chans_ipsi),[],[],map), axis off
title('After 40min of stim', 'fontweight', 'bold')
caxis([minF_hg maxF_hg]);
colorbar(gca,'WestOutside')
set(gca,'fontsize',10)

subplot(2,4,7), mapArray_new(totPwrF_stim_5(good_chans_ipsi),x(good_chans_ipsi),y(good_chans_ipsi),[],[],map), axis off
title('After 50min of stim', 'fontweight', 'bold')
caxis([minF_hg maxF_hg]);
colorbar(gca,'WestOutside')
set(gca,'fontsize',10)

subplot(2,4,8), mapArray_new(totPwrF_postStim(good_chans_ipsi),x(good_chans_ipsi),y(good_chans_ipsi),[],[],map), axis off
title('After 60min of stim', 'fontweight', 'bold')
caxis([minF_hg maxF_hg]);
colorbar(gca,'WestOutside')
set(gca,'fontsize',10)

sgtitle('Ipsilesional power')

%}

%% Calculate Time-Series Power (dt = 10s) for Monkey F

dt = 10*Fs; % number of samples every 10 seconds

% initialize variables
lngth_allbase = size(F_base,1);
lngth_allduringPT = size(F_duringPT,1);
lngth_allpostPT = size(F_postPT,1);
lngth_allstim = size(F_stim,1);
lngth_allpostStim = size(F_postStim,1);

% get the number of 10 second intervals
Ndt_allbase = floor(lngth_allbase/dt); 
Ndt_allduringPT = floor(lngth_allduringPT/dt);
Ndt_allpostPT = floor(lngth_allpostPT/dt);
Ndt_allstim = floor(lngth_allstim/dt);
Ndt_allpostStim = floor(lngth_allpostStim/dt);

pwrF_allbase = nan(Ndt_allbase,length(channels));
pwrF_allduringPT = nan(Ndt_allduringPT,length(channels));
pwrF_allpostPT = nan(Ndt_allpostPT,length(channels));
pwrF_allstim = nan(Ndt_allstim,length(channels));
pwrF_allpostStim = nan(Ndt_allpostStim,length(channels));

for idt = 1:Ndt_allbase  
    interval = ((idt-1)*dt+1):(idt*dt); % time interval in samples
    for iChan = good_chans_ipsi
        pwrF_allbase(idt,iChan) = SignalPower(F_base(interval,iChan),Fs);
    end
end

for idt = 1:Ndt_allduringPT
    interval = ((idt-1)*dt+1):(idt*dt); % time interval in samples
    for iChan = good_chans_ipsi
        pwrF_allduringPT(idt,iChan) = SignalPower(F_duringPT(interval,iChan),Fs);
    end
end


for idt = 1:Ndt_allpostPT
    interval = ((idt-1)*dt+1):(idt*dt); % time interval in samples
    for iChan = good_chans_ipsi
        pwrF_allpostPT(idt,iChan) = SignalPower(F_postPT(interval,iChan),Fs);      
    end
end

for idt = 1:Ndt_allstim
    interval = ((idt-1)*dt+1):(idt*dt); % time interval in samples
    for iChan = good_chans_ipsi
        pwrF_allstim(idt,iChan) = SignalPower(F_stim(interval,iChan),Fs);   
    end
end

for idt = 1:Ndt_allpostStim   
    interval = ((idt-1)*dt+1):(idt*dt); % time interval in samples
    for iChan = good_chans_ipsi
        pwrF_allpostStim(idt,iChan) = SignalPower(F_postStim(interval,iChan),Fs);
    end
end

%% Optional: Change Variable Names to Save Data (Monkey F)

pwrF_allbase_th = pwrF_allbase;
pwrF_allduringPT_th = pwrF_allduringPT;
pwrF_allpostPT_th = pwrF_allpostPT;
pwrF_allstim_th = pwrF_allstim;
pwrF_allpostStim_th = pwrF_allpostStim;

%% Optional: Plot Time-Series Power (dt = 10s) for Monkey F

%{
figure;
subplot(5,1,1)
plot(pwrF_allbase)

subplot(5,1,2)
plot(pwrF_allduringPT)

subplot(5,1,3)
plot(pwrF_allpostPT)

subplot(5,1,4)
plot(pwrF_allstim)

subplot(5,1,5)
plot(pwrF_allpostStim)
%}
