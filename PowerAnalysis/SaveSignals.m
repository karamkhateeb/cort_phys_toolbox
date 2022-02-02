%% Load, Downsample, and Filter PT4 and PT5 uECoG Signals

% The purpose of this script is to load, downsample, and filter the PT4 and
% PT5 ECoG signals across baseline, during illumination, and postPT
% periods.
%% 1. Initialize Variables

animal = 2; % 1 will be Monkey D, 2 will be Monkey E
period = 3; % which time period to get: 1 is baseline, 2 is duringIllumination, 3 is postPT

switch animal
    case 1
%         filepath = 'C:\Users\kkhateeb\Documents\MATLAB\PT_Neurophys\PT4_NeurophysiologyData';
    case 2
%         filepath = 'C:\Users\kkhateeb\Documents\MATLAB\PT_Neurophys\PT5_NeurophysiologyData';
end

% PT4:
% left is IPSILESIONAL: Ports A and B
% right is CONTRALESIONAL: Ports C and D
% 
% PT5:
% right is IPSILESIONAL: Ports A and B
% left is CONTRALESIONAL: Ports C and D

ports = 1:4; % Nomad port A = 1, B = 2, C = 3, D = 4
channels = 1:16; % 16 channels per port, but have 32 per hemisphere

fs = 30000; % sampling frequency in Hz (raw signal sampled at 30kS/s
downsampleFactor = 30; % Downsample by a factor of 30

disp('Initialized Variable')
%% 2. Load and Downsample Recorded Signals

switch animal
    case 1
        name = 'PT4';
        switch period
            case 1
                baselineDS.cont = cell(1,32);
                baselineDS.ipsi = cell(1,32);
            case 2
                duringPTDS.cont = cell(1,32); 
                duringPTDS.ipsi = cell(1,32);
            case 3
                postPTDS.cont = cell(1,32); 
                postPTDS.ipsi = cell(1,32); % for PT4 this will be postillumination3PT4 files
            case 4
                postPTDS1.cont = cell(1,32);
                postPTDS1.ipsi = cell(1,32); % for PT4 this will be postilluminationPT4 files
            case 5
                postPTDS2.cont = cell(1,32);
                postPTDS2.cont = cell(1,32); % for PT4 this will be postillumination2PT4 files
        end   
        
    case 2
        name = 'PT5';
        switch period
            case 1
                baselineDS.cont = cell(1,32);
                baselineDS.ipsi = cell(1,32);
            case 2
                duringPTDS.cont = cell(1,32); 
                duringPTDS.ipsi = cell(1,32);
            case 3
                postPTDS.cont = cell(1,32);
                postPTDS.ipsi = cell(1,32);
        end
end
        
tic
for iPort = ports
    
    disp(['Port ' num2str(iPort)]);
    
    for iChan = channels
        
        % Load the LFPs using Devon/Ripple's function
        % example: [x,y]=extract_all_data_function('Raw',4,true,'C:\Users\user\Trellis\dataFiles\20191113\datafile0001.ns5');
        
        dataChannelNum =128*(iPort - 1) + iChan;
        
        disp(['Channel ' num2str(iChan)]);
                
        if iPort == 1
            switch period
                case 1
                    disp(['Baseline Ipsilateral']);
                    [tempBaseRaw(:,1), tempBaseRaw(:,2)] = extract_all_data_function('Raw', dataChannelNum, false, [filepath '\baseline' name '.ns5']); % only get the first 30 minutes
                    baselineDS.ipsi{iChan} = decimate(tempBaseRaw(:,2), downsampleFactor); % Downsample the raw signal
                    toc
                case 2
                    disp(['DuringPT Ipsilateral']);
                    [tempPostRaw(:,1), tempPostRaw(:,2)] = extract_all_data_function('Raw', dataChannelNum, false, [filepath '\duringillumination' name '.ns5']); % only get the last 30 minutes
                    duringPTDS.ipsi{iChan} = decimate(tempPostRaw(:,2), downsampleFactor); % Downsample the raw signal
                    toc
                case 3
                    disp(['PostPT Ipsilateral']);
                    [tempPostRaw(:,1), tempPostRaw(:,2)] = extract_all_data_function('Raw', dataChannelNum, false, [filepath '\postillumination3' name '.ns5']); % only get the last 30 minutes
                    postPTDS.ipsi{iChan} = decimate(tempPostRaw(:,2), downsampleFactor); % Downsample the raw signal
                    toc
                case 4
                    disp(['PostPT1 Ipsilateral']);
                    [tempPostRaw(:,1), tempPostRaw(:,2)] = extract_all_data_function('Raw', dataChannelNum, false, [filepath '\postillumination' name '.ns5']); % only get the last 30 minutes
                    postPTDS1.ipsi{iChan} = decimate(tempPostRaw(:,2), downsampleFactor); % Downsample the raw signal
                    toc
                case 5
                    disp(['PostPT2 Ipsilateral']);
                    [tempPostRaw(:,1), tempPostRaw(:,2)] = extract_all_data_function('Raw', dataChannelNum, false, [filepath '\postillumination2' name '.ns5']); % only get the last 30 minutes
                    postPTDS2.ipsi{iChan} = decimate(tempPostRaw(:,2), downsampleFactor); % Downsample the raw signal
                    toc
            end
            
        elseif iPort == 2
            switch period
                case 1
                    disp(['Baseline Ipsilateral']);
                    [tempBaseRaw(:,1), tempBaseRaw(:,2)] = extract_all_data_function('Raw', dataChannelNum, false, [filepath '\baseline' name '.ns5']); % only get the first 30 minutes
                    baselineDS.ipsi{iChan+16} = decimate(tempBaseRaw(:,2), downsampleFactor);
                    toc
                case 2
                    disp(['duringPT Ipsilateral']);
                    [tempPostRaw(:,1), tempPostRaw(:,2)] = extract_all_data_function('Raw', dataChannelNum, false, [filepath '\duringillumination' name '.ns5']); % only get the last 30 minutes
                    duringPTDS.ipsi{iChan+16} = decimate(tempPostRaw(:,2), downsampleFactor);
                    toc
                case 3
                    disp(['PostPT Ipsilateral']);
                    [tempPostRaw(:,1), tempPostRaw(:,2)] = extract_all_data_function('Raw', dataChannelNum, false, [filepath '\postillumination3' name '.ns5']); % only get the last 30 minutes
                    postPTDS.ipsi{iChan+16} = decimate(tempPostRaw(:,2), downsampleFactor);
                    toc
                case 4
                    disp(['PostPT1 Ipsilateral']);
                    [tempPostRaw(:,1), tempPostRaw(:,2)] = extract_all_data_function('Raw', dataChannelNum, false, [filepath '\postillumination' name '.ns5']); % only get the last 30 minutes
                    postPTDS1.ipsi{iChan+16} = decimate(tempPostRaw(:,2), downsampleFactor);
                    toc
                case 5
                    disp(['PostPT2 Ipsilateral']);
                    [tempPostRaw(:,1), tempPostRaw(:,2)] = extract_all_data_function('Raw', dataChannelNum, false, [filepath '\postillumination2' name '.ns5']); % only get the last 30 minutes
                    postPTDS2.ipsi{iChan+16} = decimate(tempPostRaw(:,2), downsampleFactor);
                    toc
            end
        elseif iPort == 3
            switch period
                case 1
                    disp(['Baseline Contralateral']);
                    [tempBaseRaw(:,1), tempBaseRaw(:,2)] = extract_all_data_function('Raw', dataChannelNum, false, [filepath '\baseline' name '.ns5']); % only get the first 30 minutes
                    baselineDS.cont{iChan} = decimate(tempBaseRaw(:,2), downsampleFactor);
                    toc
                case 2
                    disp(['DuringPT Contralateral']);
                    [tempPostRaw(:,1), tempPostRaw(:,2)] = extract_all_data_function('Raw', dataChannelNum, false, [filepath '\duringillumination' name '.ns5']); % only get the last 30 minutes
                    duringPTDS.cont{iChan} = decimate(tempPostRaw(:,2), downsampleFactor);
                    toc
                case 3
                    disp(['PostPT Contralateral']);
                    [tempPostRaw(:,1), tempPostRaw(:,2)] = extract_all_data_function('Raw', dataChannelNum, false, [filepath '\postillumination3' name '.ns5']); % only get the last 30 minutes
                    postPTDS.cont{iChan} = decimate(tempPostRaw(:,2), downsampleFactor);
                    toc
                case 4
                    disp(['PostPT1 Contralateral']);
                    [tempPostRaw(:,1), tempPostRaw(:,2)] = extract_all_data_function('Raw', dataChannelNum, false, [filepath '\postillumination' name '.ns5']); % only get the last 30 minutes
                    postPTDS1.cont{iChan} = decimate(tempPostRaw(:,2), downsampleFactor);
                    toc
                case 5
                    disp(['PostPT2 Contralateral']);
                    [tempPostRaw(:,1), tempPostRaw(:,2)] = extract_all_data_function('Raw', dataChannelNum, false, [filepath '\postillumination2' name '.ns5']); % only get the last 30 minutes
                    postPTDS2.cont{iChan} = decimate(tempPostRaw(:,2), downsampleFactor);
                    toc 
            end
            
        elseif iPort == 4
            switch period
                case 1
                    disp(['Baseline Contralateral']);
                    [tempBaseRaw(:,1), tempBaseRaw(:,2)] = extract_all_data_function('Raw', dataChannelNum, false, [filepath '\baseline' name '.ns5']); % only get the first 30 minutes
                    baselineDS.cont{iChan+16} = decimate(tempBaseRaw(:,2), downsampleFactor);
                    toc
                case 2
                    disp(['duringPT Contralateral']);
                    [tempPostRaw(:,1), tempPostRaw(:,2)] = extract_all_data_function('Raw', dataChannelNum, false, [filepath '\duringillumination' name '.ns5']); % only get the last 30 minutes
                    duringPTDS.cont{iChan+16} = decimate(tempPostRaw(:,2), downsampleFactor);
                    toc
                case 3
                    disp(['PostPT Contralateral']);
                    [tempPostRaw(:,1), tempPostRaw(:,2)] = extract_all_data_function('Raw', dataChannelNum, false, [filepath '\postillumination3' name '.ns5']); % only get the last 30 minutes
                    postPTDS.cont{iChan+16} = decimate(tempPostRaw(:,2), downsampleFactor);
                    toc
                case 4
                    disp(['PostPT1 Contralateral']);
                    [tempPostRaw(:,1), tempPostRaw(:,2)] = extract_all_data_function('Raw', dataChannelNum, false, [filepath '\postillumination' name '.ns5']); % only get the last 30 minutes
                    postPTDS1.cont{iChan+16} = decimate(tempPostRaw(:,2), downsampleFactor);
                    toc
                case 5
                    disp(['PostP2 Contralateral']);
                    [tempPostRaw(:,1), tempPostRaw(:,2)] = extract_all_data_function('Raw', dataChannelNum, false, [filepath '\postillumination2' name '.ns5']); % only get the last 30 minutes
                    postPTDS2.cont{iChan+16} = decimate(tempPostRaw(:,2), downsampleFactor);
                    toc
            end
        end
        
        clearvars tempBaseRaw; clearvars tempPostRaw;
        
    end
end

% Save the cells
switch animal
    case 1, filename = 'DownsampledSignalsPT4';
    case 2, filename = 'DownsampledSignalsPT5';
end

switch period
    case 1, variable = 'baselineDS';
    case 2, variable = 'duringPTDS';
    case 3, variable = 'postPTDS';
    case 4, variable = 'postPTDS1';
    case 5, variable = 'postPTDS2';
end

disp(['Saving ' variable])

save([filepath filesep filename], variable, '-v7.3', '-nocompression','-append');

disp('Saved')
%% 3. Filter the Signal

% determine the desired frequency bad
FHighPass = 1; % Low cutoff in Hz
FLowPass = 3; % High cutoff in Hz

band = [num2str(FHighPass) '_' num2str(FLowPass)]; % for variable naming
switch period
    case 1, filtVar = 'baselineBPF';
    case 2, filtVar = 'duringPTBPF';
    case 3, filtVar = 'postPTBPF';
    case 4, filtVar = 'postPT1BPF';
    case 5, filtVar = 'postPT2BPF';
end

% Design notch filters
d60 = designfilt('bandstopiir', 'FilterOrder', 2, 'HalfPowerFrequency1',...
    59.75, 'HalfPowerFrequency2', 60.25, 'DesignMethod', 'butter',...
    'SampleRate', fs/downsampleFactor);
d120 = designfilt('bandstopiir', 'FilterOrder', 2, 'HalfPowerFrequency1',...
    119.75, 'HalfPowerFrequency2', 120.25, 'DesignMethod', 'butter',...
    'SampleRate', fs/downsampleFactor);
d180 = designfilt('bandstopiir', 'FilterOrder', 2, 'HalfPowerFrequency1',...
    179.75, 'HalfPowerFrequency2', 180.25, 'DesignMethod', 'butter',...
    'SampleRate', fs/downsampleFactor);
d240 = designfilt('bandstopiir', 'FilterOrder', 2, 'HalfPowerFrequency1',...
    239.75, 'HalfPowerFrequency2', 240.25, 'DesignMethod', 'butter',...
    'SampleRate', fs/downsampleFactor);
        
                
tic
hemiChannels = 1:32; % 32 channels per hemisphere
for iChan = hemiChannels        
        
        disp(['Channel ' num2str(iChan)]);
        
        switch period
            case 1
                % Notch filter for 60 Hz
                baselineNotch.ipsi{iChan} = filtfilt(d60, baselineDS.ipsi{iChan});
                baselineNotch.cont{iChan} = filtfilt(d60, baselineDS.cont{iChan});
                
                % Notch filter for 120 Hz
                baselineNotch.ipsi{iChan} = filtfilt(d120, baselineNotch.ipsi{iChan});
                baselineNotch.cont{iChan} = filtfilt(d120, baselineNotch.cont{iChan});
                
                % Notch filter for 180 Hz
                baselineNotch.ipsi{iChan} = filtfilt(d180, baselineNotch.ipsi{iChan});
                baselineNotch.cont{iChan} = filtfilt(d180, baselineNotch.cont{iChan});
                
                % Notch filter for 240 Hz
                baselineNotch.ipsi{iChan} = filtfilt(d240, baselineNotch.ipsi{iChan});
                baselineNotch.cont{iChan} = filtfilt(d240, baselineNotch.cont{iChan});
                
                % Filter the signal to the desired frequency band
                baselineBPF1_3.ipsi{iChan} = bandpass(baselineNotch.ipsi{iChan},...
                    [FHighPass FLowPass], fs/downsampleFactor,...
                    'ImpulseResponse', 'iir', 'Steepness', 0.95);
                baselineBPF1_3.cont{iChan} = bandpass(baselineNotch.cont{iChan},...
                    [FHighPass FLowPass], fs/downsampleFactor,...
                    'ImpulseResponse', 'iir', 'Steepness', 0.95);
            case 2
                % Notch filter for 60 Hz
                duringPTNotch.ipsi{iChan} = filtfilt(d60, duringPTDS.ipsi{iChan});
                duringPTNotch.cont{iChan} = filtfilt(d60, duringPTDS.cont{iChan});
                
                % Notch filter for 120 Hz
                duringPTNotch.ipsi{iChan} = filtfilt(d120, duringPTNotch.ipsi{iChan});
                duringPTNotch.cont{iChan} = filtfilt(d120, duringPTNotch.cont{iChan});
                
                % Notch filter for 180 Hz
                duringPTNotch.ipsi{iChan} = filtfilt(d180, duringPTNotch.ipsi{iChan});
                duringPTNotch.cont{iChan} = filtfilt(d180, duringPTNotch.cont{iChan});
                
                % Notch filter for 240 Hz
                duringPTNotch.ipsi{iChan} = filtfilt(d240, duringPTNotch.ipsi{iChan});
                duringPTNotch.cont{iChan} = filtfilt(d240, duringPTNotch.cont{iChan});
                
                % Filter the signal to the desired frequency band
                duringPTBPF1_3.ipsi{iChan} = bandpass(duringPTNotch.ipsi{iChan},...
                    [FHighPass FLowPass], fs/downsampleFactor,...
                    'ImpulseResponse', 'iir', 'Steepness', 0.95);
                duringPTBPF1_3.cont{iChan} = bandpass(duringPTNotch.cont{iChan},...
                    [FHighPass FLowPass], fs/downsampleFactor,...
                    'ImpulseResponse', 'iir', 'Steepness', 0.95);
            case 3
                % Notch filter for 60 Hz
                postPTNotch.ipsi{iChan} = filtfilt(d60, postPTDS.ipsi{iChan});
                postPTNotch.cont{iChan} = filtfilt(d60, postPTDS.cont{iChan});
                
                % Notch filter for 120 Hz
                postPTNotch.ipsi{iChan} = filtfilt(d120, postPTNotch.ipsi{iChan});
                postPTNotch.cont{iChan} = filtfilt(d120, postPTNotch.cont{iChan});
                
                % Notch filter for 180 Hz
                postPTNotch.ipsi{iChan} = filtfilt(d180, postPTNotch.ipsi{iChan});
                postPTNotch.cont{iChan} = filtfilt(d180, postPTNotch.cont{iChan});
                
                % Notch filter for 240 Hz
                postPTNotch.ipsi{iChan} = filtfilt(d240, postPTNotch.ipsi{iChan});
                postPTNotch.cont{iChan} = filtfilt(d240, postPTNotch.cont{iChan});
                
                % Filter the signal to the desired frequency band
                postPTBPF1_3.ipsi{iChan} = bandpass(postPTNotch.ipsi{iChan},...
                    [FHighPass FLowPass], fs/downsampleFactor,...
                    'ImpulseResponse', 'iir', 'Steepness', 0.95);
                postPTBPF1_3.cont{iChan} = bandpass(postPTNotch.cont{iChan},...
                    [FHighPass FLowPass], fs/downsampleFactor,...
                    'ImpulseResponse', 'iir', 'Steepness', 0.95);
            case 4
                % Notch filter for 60 Hz
                postPT1Notch.ipsi{iChan} = filtfilt(d60, postPTDS1.ipsi{iChan});
                postPT1Notch.cont{iChan} = filtfilt(d60, postPTDS1.cont{iChan});
                
                % Notch filter for 120 Hz
                postPT1Notch.ipsi{iChan} = filtfilt(d120, postPT1Notch.ipsi{iChan});
                postPT1Notch.cont{iChan} = filtfilt(d120, postPT1Notch.cont{iChan});
                
                % Notch filter for 180 Hz
                postPT1Notch.ipsi{iChan} = filtfilt(d180, postPT1Notch.ipsi{iChan});
                postPT1Notch.cont{iChan} = filtfilt(d180, postPT1Notch.cont{iChan});
                
                % Notch filter for 240 Hz
                postPT1Notch.ipsi{iChan} = filtfilt(d240, postPT1Notch.ipsi{iChan});
                postPT1Notch.cont{iChan} = filtfilt(d240, postPT1Notch.cont{iChan});
                
                % Filter the signal to the desired frequency band
                postPT1BPF1_3.ipsi{iChan} = bandpass(postPT1Notch.ipsi{iChan},...
                    [FHighPass FLowPass], fs/downsampleFactor,...
                    'ImpulseResponse', 'iir', 'Steepness', 0.95);
                postPT1BPF1_3.cont{iChan} = bandpass(postPT1Notch.cont{iChan},...
                    [FHighPass FLowPass], fs/downsampleFactor,...
                    'ImpulseResponse', 'iir', 'Steepness', 0.95);
            case 5
                % Notch filter for 60 Hz
                postPT2Notch.ipsi{iChan} = filtfilt(d60, postPTDS2.ipsi{iChan});
                postPT2Notch.cont{iChan} = filtfilt(d60, postPTDS2.cont{iChan});
                
                % Notch filter for 120 Hz
                postPT2Notch.ipsi{iChan} = filtfilt(d120, postPT2Notch.ipsi{iChan});
                postPT2Notch.cont{iChan} = filtfilt(d120, postPT2Notch.cont{iChan});
                
                % Notch filter for 180 Hz
                postPT2Notch.ipsi{iChan} = filtfilt(d180, postPT2Notch.ipsi{iChan});
                postPT2Notch.cont{iChan} = filtfilt(d180, postPT2Notch.cont{iChan});
                
                % Notch filter for 240 Hz
                postPT2Notch.ipsi{iChan} = filtfilt(d240, postPT2Notch.ipsi{iChan});
                postPT2Notch.cont{iChan} = filtfilt(d240, postPT2Notch.cont{iChan});
                
                % Filter the signal to the desired frequency band
                postPT2BPF1_3.ipsi{iChan} = bandpass(postPT2Notch.ipsi{iChan},...
                    [FHighPass FLowPass], fs/downsampleFactor,...
                    'ImpulseResponse', 'iir', 'Steepness', 0.95);
                postPT2BPF1_3.cont{iChan} = bandpass(postPT2Notch.cont{iChan},...
                    [FHighPass FLowPass], fs/downsampleFactor,...
                    'ImpulseResponse', 'iir', 'Steepness', 0.95);
        end
        toc
end

% Save filtered signals for later use
disp('Saving Filtered Signal')
% save([filepath filesep 'FilteredSignals' name '.mat'], [filtVar band],...
%     '-v7.3', '-nocompression','-append');
save(['FilteredSignals' name '.mat'], [filtVar band],...
    '-v7.3', '-nocompression','-append');
disp('Saved Filtered Signals')