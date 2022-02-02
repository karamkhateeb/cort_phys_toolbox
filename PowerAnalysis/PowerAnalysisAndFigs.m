%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clustering Neural Activity Following Stroke to Infer Lesion Location in
% Macaque Sensorimotor Cortex
% Karam Khateeb
% March 31, 2020 - March 25, 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Initialize Important Variables

Fs = 1000; % downsampled frequency
limit30 = 30*60*Fs; % number of samples in 30 minutes
channels = 1:32; % there are 32 channels

% Gate to stop from accidentally loading all signals
getAllSigs = false;

% What do I want to plot?
plotTS = true; % Plot time-series signal before and after artifact removal
plotPwr = true; % Plot time-averaged signal power over time

% X and Y locations of electrodes in array
x = [4 4 3 2 2 2 1 1 1 1 2 2 3 3 3 4 4 4 5 6 6 7 8 7 7 7 6 6 6 5 5 5]*0.75;
y = [3 1 2 1 5 3 4 2 8 6 9 7 6 4 8 5 7 9 8 9 7 8 5 6 2 4 3 5 1 6 4 2]*0.75;

% Number of Clusters
numClusters = 3;

% Choose a colormap for plotting all of the figures
% cmp = colormap(parula(3)); close gcf
cmp = [248 117 117;64 78 77;8 178 227]/255;
map = [linspace(248,64,64) linspace(64,8,64);
    linspace(117,78,64) linspace(78,178,64);
    linspace(117,77,64) linspace(77,227,64)]' / 255;

disp('Initialized Important Variables')
%% 2. Get PT4 and PT5 baseline and postPT Ipsilesional Recordings

% load filtered PT4 and PT5 signals saved from PTNeurophysPower.m or
% SaveSignals.m
tic
% PT4 (monkey D)
load('C:\Users\kkhateeb\Documents\MATLAB\PT_Neurophys\PT4_NeurophysiologyData\FilteredSignalsPT4.mat');

signalD_base_hg = zeros(limit30,length(channels));
signalD_base_lg = zeros(limit30,length(channels));
signalD_base_bt = zeros(limit30,length(channels));
signalD_base_th = zeros(limit30,length(channels));
signalD_post_hg = zeros(limit30,length(channels));
signalD_post_lg = zeros(limit30,length(channels));
signalD_post_bt = zeros(limit30,length(channels));
signalD_post_th = zeros(limit30,length(channels));

for iChan = channels
    
    disp(['Channel ' num2str(iChan)])
    
    % Baseline Signals - first 30 minutes
    signalD_base_hg(:,iChan) = baselineBPF60_150.ipsi{1,iChan}(1:limit30); % high gamma signal
    signalD_base_lg(:,iChan) = baselineBPF30_59.ipsi{1,iChan}(1:limit30); % low gamma signal
    signalD_base_bt(:,iChan) = baselineBPF12_29.ipsi{1,iChan}(1:limit30); % beta signal
    signalD_base_th(:,iChan) = baselineBPF4_7.ipsi{1,iChan}(1:limit30); % theta signal

    % Post PT Signals - last 30 minutes
    signalD_post_hg(:,iChan) = postPTBPF60_150.ipsi{1,iChan}((end-limit30+1):end); % high gamma signal
    signalD_post_lg(:,iChan) = postPTBPF30_59.ipsi{1,iChan}((end-limit30+1):end); % low gamma signal
    signalD_post_bt(:,iChan) = postPTBPF12_29.ipsi{1,iChan}((end-limit30+1):end); % beta signal
    signalD_post_th(:,iChan) = postPTBPF4_7.ipsi{1,iChan}((end-limit30+1):end); % theta signal
    toc
end

disp('Loaded PT4 Signals')
%%
% PT5 (monkey E)
load('C:\Users\kkhateeb\Documents\MATLAB\PT_Neurophys\PT5_NeurophysiologyData\FilteredSignalsPT5.mat');

signalE_base_hg = zeros(limit30,length(channels));
signalE_base_lg = zeros(limit30,length(channels));
signalE_base_bt = zeros(limit30,length(channels));
signalE_base_th = zeros(limit30,length(channels));
signalE_post_hg = zeros(limit30,length(channels));
signalE_post_lg = zeros(limit30,length(channels));
signalE_post_bt = zeros(limit30,length(channels));
signalE_post_th = zeros(limit30,length(channels));

for iChan = channels
    
    disp(['Channel ' num2str(iChan)])
    
    % Baseline Signals - first 30 minutes
    signalE_base_hg(:,iChan) = baselineBPF60_150.ipsi{1,iChan}(1:limit30); % high gamma signal
    signalE_base_lg(:,iChan) = baselineBPF30_59.ipsi{1,iChan}(1:limit30); % low gamma signal
    signalE_base_bt(:,iChan) = baselineBPF12_29.ipsi{1,iChan}(1:limit30); % beta signal
    signalE_base_th(:,iChan) = baselineBPF4_7.ipsi{1,iChan}(1:limit30); % theta signal

    % Post PT Signals - last 30 minutes
    signalE_post_hg(:,iChan) = postPTBPF60_150.ipsi{1,iChan}((end-limit30+1):end); % high gamma signal
    signalE_post_lg(:,iChan) = postPTBPF30_59.ipsi{1,iChan}((end-limit30+1):end); % low gamma signal
    signalE_post_bt(:,iChan) = postPTBPF12_29.ipsi{1,iChan}((end-limit30+1):end); % beta signal
    signalE_post_th(:,iChan) = postPTBPF4_7.ipsi{1,iChan}((end-limit30+1):end); % theta signal
    toc
end

disp('Loaded PT5 Signals')

clearvars -REGEXP baseline* postPT* during*

%% 3. Plot the Time Series Signal

if plotTS
    figure, set(gcf, 'color', 'w', 'units', 'normalized', 'outerposition', [0 0 1 1])
    subplot(4,2,1), plot(signalD_base_hg), box off, ylabel('High Gamma'), title('Baseline')
    subplot(4,2,2), plot(signalD_post_hg), box off, title('Post PT')
    subplot(4,2,3), plot(signalD_base_lg), box off, ylabel('Low Gamma')
    subplot(4,2,4), plot(signalD_post_lg), box off
    subplot(4,2,5), plot(signalD_base_bt), box off, ylabel('Beta')
    subplot(4,2,6), plot(signalD_post_bt), box off
    subplot(4,2,7), plot(signalD_base_th), box off, ylabel('Theta')
    subplot(4,2,8), plot(signalD_post_th), box off
    [~,hx] = suplabel('Time', 'x');
    [~,ht] = suplabel('Monkey D Signal', 't');
    
    figure, set(gcf, 'color', 'w', 'units', 'normalized', 'outerposition', [0 0 1 1])
    subplot(4,2,1), plot(signalE_base_hg), box off, ylabel('High Gamma'), title('Baseline')
    subplot(4,2,2), plot(signalE_post_hg), box off, title('Post PT')
    subplot(4,2,3), plot(signalE_base_lg), box off, ylabel('Low Gamma')
    subplot(4,2,4), plot(signalE_post_lg), box off
    subplot(4,2,5), plot(signalE_base_bt), box off, ylabel('Beta')
    subplot(4,2,6), plot(signalE_post_bt), box off
    subplot(4,2,7), plot(signalE_base_th), box off, ylabel('Theta')
    subplot(4,2,8), plot(signalE_post_th), box off
    [~,hx] = suplabel('Time', 'x');
    [~,ht] = suplabel('Monkey E Signal', 't');
end

%% 4. Clean up the Signal (Remove Artifacts)

% Find ranges of noisy data to be removed using the findArtifacts.m script
noiseD_base = [
    345600, 346100;
    1536000, 1538000; 
    1, 700];

noiseD_post = [354300, 354700];

noiseE_base = [22130, 22960;
    704200, 705300;
    916500, 917100;
    1674000, 1675000;
    1723000, 1724000;
    262700, 263900;
    1, 1150];

noiseE_post = [1, 4900;
    42540, 45380;
    130800, 158600;
    170700, 177400;
    434600, 438000;
    442600, 444000;
    500600, 528300;
    732800, 734900;
    1485000, 1488000;
    1502000, 1508500;
    1640000, 1645000;
    1655000, 1656500;
    1659000, 1668000;
    1677000, 1682000;
    1689000, 1691000;
    1711150, 1714000;
    1715500, 1717500;
    1718000, 1721000;
    1722500, 1727500;
    1745000, 1753000];

% Save artifact indices for later use
% save('ArtifactIndices.mat','noiseD_base','noiseD_post','noiseE_base','noiseE_post');

cleanD_base_hg = signalD_base_hg;
cleanD_base_lg = signalD_base_lg;
cleanD_base_bt = signalD_base_bt;
cleanD_base_th = signalD_base_th;
cleanD_post_hg = signalD_post_hg;
cleanD_post_lg = signalD_post_lg;
cleanD_post_bt = signalD_post_bt;
cleanD_post_th = signalD_post_th;

cleanE_base_hg = signalE_base_hg;
cleanE_base_lg = signalE_base_lg;
cleanE_base_bt = signalE_base_bt;
cleanE_base_th = signalE_base_th;
cleanE_post_hg = signalE_post_hg;
cleanE_post_lg = signalE_post_lg;
cleanE_post_bt = signalE_post_bt;
cleanE_post_th = signalE_post_th;

% replace noisey data with NaN
for i = 1:size(noiseD_base,1)
    range = noiseD_base(i,1):noiseD_base(i,2);
    numSamps = 1+noiseD_base(i,2)-noiseD_base(i,1);
    cleanD_base_hg(range,:) = NaN*ones(numSamps,length(channels));
    cleanD_base_lg(range,:) = NaN*ones(numSamps,length(channels));
    cleanD_base_bt(range,:) = NaN*ones(numSamps,length(channels));
    cleanD_base_th(range,:) = NaN*ones(numSamps,length(channels));
    
end

for i = 1:size(noiseD_post,1)
    range = noiseD_post(i,1):noiseD_post(i,2);
    numSamps = 1+noiseD_post(i,2)-noiseD_post(i,1);
    cleanD_post_hg(range,:) = NaN*ones(numSamps,length(channels));
    cleanD_post_lg(range,:) = NaN*ones(numSamps,length(channels));
    cleanD_post_bt(range,:) = NaN*ones(numSamps,length(channels));
    cleanD_post_th(range,:) = NaN*ones(numSamps,length(channels));
end

for i = 1:size(noiseE_base,1)
    range = noiseE_base(i,1):noiseE_base(i,2);
    numSamps = 1+noiseE_base(i,2)-noiseE_base(i,1);
    cleanE_base_hg(range,:) = NaN*ones(numSamps,length(channels));
    cleanE_base_lg(range,:) = NaN*ones(numSamps,length(channels));
    cleanE_base_bt(range,:) = NaN*ones(numSamps,length(channels));
    cleanE_base_th(range,:) = NaN*ones(numSamps,length(channels));
end

for i = 1:size(noiseE_post,1)
    range = noiseE_post(i,1):noiseE_post(i,2);
    numSamps = 1+noiseE_post(i,2)-noiseE_post(i,1);
    cleanE_post_hg(range,:) = NaN*ones(numSamps,length(channels));
    cleanE_post_lg(range,:) = NaN*ones(numSamps,length(channels));
    cleanE_post_bt(range,:) = NaN*ones(numSamps,length(channels));
    cleanE_post_th(range,:) = NaN*ones(numSamps,length(channels));
end

% save cleaned signals for later use
% save('CleanedSignals.mat','cleanD_base_hg','cleanD_base_lg',...
%     'cleanD_base_bt','cleanD_base_th','cleanD_post_hg','cleanD_post_lg',...
%     'cleanD_post_bt','cleanD_post_th','cleanE_base_hg','cleanE_base_lg',...
%     'cleanE_base_bt','cleanE_base_th','cleanE_post_hg','cleanE_post_lg',...
%     'cleanE_post_bt','cleanE_post_th');


disp('Cleaned up Signal')
%% 5. Plot the Time Series Signal

if plotTS
    figure, set(gcf, 'color', 'w', 'units', 'normalized', 'outerposition', [0 0 1 1])
    subplot(4,2,1), plot(cleanD_base_hg), box off, title('Baseline'), ylabel('High Gamma')
    subplot(4,2,2), plot(cleanD_post_hg), box off, title('Post PT')
    subplot(4,2,3), plot(cleanD_base_lg), box off, ylabel('Low Gamma')
    subplot(4,2,4), plot(cleanD_post_lg), box off
    subplot(4,2,5), plot(cleanD_base_bt), box off, ylabel('Beta')
    subplot(4,2,6), plot(cleanD_post_bt), box off
    subplot(4,2,7), plot(cleanD_base_th), box off, ylabel('Theta')
    subplot(4,2,8), plot(cleanD_post_th), box off
    [~,hx] = suplabel('Time', 'x');
    [~,ht] = suplabel('Monkey D Signal', 't');
    
    figure, set(gcf, 'color', 'w', 'units', 'normalized', 'outerposition', [0 0 1 1])
    subplot(4,2,1), plot(cleanE_base_hg), box off, title('Baseline'), ylabel('High Gamma')
    subplot(4,2,2), plot(cleanE_post_hg), box off, title('Post PT')
    subplot(4,2,3), plot(cleanE_base_lg), box off, ylabel('Low Gamma')
    subplot(4,2,4), plot(cleanE_post_lg), box off
    subplot(4,2,5), plot(cleanE_base_bt), box off, ylabel('Beta')
    subplot(4,2,6), plot(cleanE_post_bt), box off
    subplot(4,2,7), plot(cleanE_base_th), box off, ylabel('Theta')
    subplot(4,2,8), plot(cleanE_post_th), box off
    [~,hx] = suplabel('Time', 'x');
    [~,ht] = suplabel('Monkey E Signal', 't');
end

%% 6. Take Signal Power Across Entire 30 minute Time Period

totPwrD_base_hg = zeros(1,length(channels));
totPwrD_base_lg = zeros(1,length(channels));
totPwrD_base_bt = zeros(1,length(channels));
totPwrD_base_th = zeros(1,length(channels));
totPwrD_post_hg = zeros(1,length(channels));
totPwrD_post_lg = zeros(1,length(channels));
totPwrD_post_bt = zeros(1,length(channels));
totPwrD_post_th = zeros(1,length(channels));

totPwrE_base_hg = zeros(1,length(channels));
totPwrE_base_lg = zeros(1,length(channels));
totPwrE_base_bt = zeros(1,length(channels));
totPwrE_base_th = zeros(1,length(channels));
totPwrE_post_hg = zeros(1,length(channels));
totPwrE_post_lg = zeros(1,length(channels));
totPwrE_post_bt = zeros(1,length(channels));
totPwrE_post_th = zeros(1,length(channels));

for iChan = channels
        
    totPwrD_base_hg(iChan) = SignalPower(cleanD_base_hg(:,iChan),Fs);
    totPwrD_base_lg(iChan) = SignalPower(cleanD_base_lg(:,iChan),Fs);
    totPwrD_base_bt(iChan) = SignalPower(cleanD_base_bt(:,iChan),Fs);
    totPwrD_base_th(iChan) = SignalPower(cleanD_base_th(:,iChan),Fs);
    totPwrD_post_hg(iChan) = SignalPower(cleanD_post_hg(:,iChan),Fs);
    totPwrD_post_lg(iChan) = SignalPower(cleanD_post_lg(:,iChan),Fs);
    totPwrD_post_bt(iChan) = SignalPower(cleanD_post_bt(:,iChan),Fs);
    totPwrD_post_th(iChan) = SignalPower(cleanD_post_th(:,iChan),Fs);

    totPwrE_base_hg(iChan) = SignalPower(cleanE_base_hg(:,iChan),Fs);
    totPwrE_base_lg(iChan) = SignalPower(cleanE_base_lg(:,iChan),Fs);
    totPwrE_base_bt(iChan) = SignalPower(cleanE_base_bt(:,iChan),Fs);
    totPwrE_base_th(iChan) = SignalPower(cleanE_base_th(:,iChan),Fs);
    totPwrE_post_hg(iChan) = SignalPower(cleanE_post_hg(:,iChan),Fs);
    totPwrE_post_lg(iChan) = SignalPower(cleanE_post_lg(:,iChan),Fs);
    totPwrE_post_bt(iChan) = SignalPower(cleanE_post_bt(:,iChan),Fs);
    totPwrE_post_th(iChan) = SignalPower(cleanE_post_th(:,iChan),Fs);
end

% Remove Channel 8 from Monkey D analysis and channels 8 and 24 for monk E
% Bad channels determined using VisBadChannels.m script
totPwrD_base_hg(8) = nan; totPwrD_base_lg(8) = nan;
totPwrD_base_bt(8) = nan; totPwrD_base_th(8) = nan;
totPwrD_post_hg(8) = nan; totPwrD_post_lg(8) = nan;
totPwrD_post_bt(8) = nan; totPwrD_post_th(8) = nan;

totPwrE_base_hg(8) = nan; totPwrE_base_lg(8) = nan;
totPwrE_base_bt(8) = nan; totPwrE_base_th(8) = nan;
totPwrE_post_hg(8) = nan; totPwrE_post_lg(8) = nan;
totPwrE_post_bt(8) = nan; totPwrE_post_th(8) = nan;

totPwrE_base_hg(24) = nan; totPwrE_base_lg(24) = nan;
totPwrE_base_bt(24) = nan; totPwrE_base_th(24) = nan;
totPwrE_post_hg(24) = nan; totPwrE_post_lg(24) = nan;
totPwrE_post_bt(24) = nan; totPwrE_post_th(24) = nan;

totPwrE_base_hg(25) = nan; totPwrE_base_lg(25) = nan;
totPwrE_base_bt(25) = nan; totPwrE_base_th(25) = nan;
totPwrE_post_hg(25) = nan; totPwrE_post_lg(25) = nan;
totPwrE_post_bt(25) = nan; totPwrE_post_th(25) = nan;

minD_hg = min([totPwrD_base_hg,totPwrD_post_hg]);
minD_lg = min([totPwrD_base_lg,totPwrD_post_lg]);
minD_bt = min([totPwrD_base_bt,totPwrD_post_bt]);
minD_th = min([totPwrD_base_th,totPwrD_post_th]);
maxD_hg = max([totPwrD_base_hg,totPwrD_post_hg]);
maxD_lg = max([totPwrD_base_lg,totPwrD_post_lg]);
maxD_bt = max([totPwrD_base_bt,totPwrD_post_bt]);
maxD_th = max([totPwrD_base_th,totPwrD_post_th]);

minE_hg = min([totPwrE_base_hg,totPwrE_post_hg]);
minE_lg = min([totPwrE_base_lg,totPwrE_post_lg]);
minE_bt = min([totPwrE_base_bt,totPwrE_post_bt]);
minE_th = min([totPwrE_base_th,totPwrE_post_th]);
maxE_hg = max([totPwrE_base_hg,totPwrE_post_hg]);
maxE_lg = max([totPwrE_base_lg,totPwrE_post_lg]);
maxE_bt = max([totPwrE_base_bt,totPwrE_post_bt]);
maxE_th = max([totPwrE_base_th,totPwrE_post_th]);

% Show Monkey D
figure, set(gcf, 'color', 'w', 'units', 'normalize', 'outerposition', [.3 0 .3 1])
subplot(4,2,1), mapArray(totPwrD_base_hg)
title('Baseline'), ylabel('High Gamma', 'fontweight', 'bold');
caxis([minD_hg maxD_hg]);

subplot(4,2,2), mapArray(totPwrD_post_hg), axis off
title('Post PT', 'fontweight', 'bold')
caxis([minD_hg maxD_hg]);

subplot(4,2,3), mapArray(totPwrD_base_lg)
ylabel('Low Gamma', 'fontweight', 'bold');
caxis([minD_lg maxD_lg]); 

subplot(4,2,4), mapArray(totPwrD_post_lg), axis off
caxis([minD_lg maxD_lg]);

subplot(4,2,5), mapArray(totPwrD_base_bt)
ylabel('Beta', 'fontweight', 'bold');
caxis([minD_bt maxD_bt]);

subplot(4,2,6), mapArray(totPwrD_post_bt), axis off
caxis([minD_bt maxD_bt]);

subplot(4,2,7), mapArray(totPwrD_base_th)
ylabel('Theta', 'fontweight', 'bold');
caxis([minD_th maxD_th]);

subplot(4,2,8), mapArray(totPwrD_post_th), axis off
caxis([minD_th maxD_th]);

[~,ht] = suplabel('Monkey D Total Power','t',[0.08 .125 .84 .84]);
set(ht,'fontsize',18);

% Show Monkey E
figure, set(gcf, 'color', 'w', 'units', 'normalize', 'outerposition', [.3 0 .3 1])
subplot(4,2,1), mapArray(totPwrE_base_hg)
title('Baseline', 'fontweight', 'bold')
ylabel('High Gamma', 'fontweight', 'bold');
caxis([minE_hg maxE_hg]);

subplot(4,2,2), mapArray(totPwrE_post_hg), axis off
title('Post PT', 'fontweight', 'bold')
caxis([minE_hg maxE_hg]);

subplot(4,2,3), mapArray(totPwrE_base_lg)
ylabel('Low Gamma', 'fontweight', 'bold');
caxis([minE_lg maxE_lg]);

subplot(4,2,4), mapArray(totPwrE_post_lg), axis off
caxis([minE_lg maxE_lg]);

subplot(4,2,5), mapArray(totPwrE_base_bt)
ylabel('Beta', 'fontweight', 'bold');
caxis([minE_bt maxE_bt]);

subplot(4,2,6), mapArray(totPwrE_post_bt), axis off
caxis([minE_bt maxE_bt]);

subplot(4,2,7), mapArray(totPwrE_base_th)
ylabel('Theta', 'fontweight', 'bold');
caxis([minE_th maxE_th]);

subplot(4,2,8), mapArray(totPwrE_post_th), axis off
caxis([minE_th maxE_th]);

[~,ht] = suplabel('Monkey E Total Power','t',[0.08 .125 .84 .84]);
set(ht,'fontsize',18);

%% 7. Calculate the Changes in Power (dt = 30 min)

% Subtract baseline powers from postPT powers
subTotD_hg = totPwrD_post_hg - totPwrD_base_hg;
subTotD_lg = totPwrD_post_lg - totPwrD_base_lg;
subTotD_bt = totPwrD_post_bt - totPwrD_base_bt;
subTotD_th = totPwrD_post_th - totPwrD_base_th;

subTotE_hg = totPwrE_post_hg - totPwrE_base_hg;
subTotE_lg = totPwrE_post_lg - totPwrE_base_lg;
subTotE_bt = totPwrE_post_bt - totPwrE_base_bt;
subTotE_th = totPwrE_post_th - totPwrE_base_th;

disp('Calculated Power Changes')

%% 8. Visualize Power (dt = 30 min) Changes

% Visualize subtracted power values
figure('color','w','units','normalized','outerposition',[.3 0 .3 1])
subplot(4,2,1), mapArray(subTotD_hg), title('Monkey D'), ylabel('High Gamma')
subplot(4,2,2), mapArray(subTotE_hg), title('Monkey E'), axis off
subplot(4,2,3), mapArray(subTotD_lg), ylabel('Low Gamma')
subplot(4,2,4), mapArray(subTotE_lg), axis off
subplot(4,2,5), mapArray(subTotD_bt), ylabel('Beta')
subplot(4,2,6), mapArray(subTotE_bt), axis off
subplot(4,2,7), mapArray(subTotD_th), ylabel('Theta')
subplot(4,2,8), mapArray(subTotE_th), axis off

[~,ht] = suplabel('Subtracted Power Values','t',[.08 .125 .84 .84]);
set(ht,'fontsize',18)

%% 9. Load all Filtered Signals for PT4

% Don't want to accidentally run this
if getAllSigs
    % load filtered PT4 and PT5 signals saved from SaveSignals.m to make 
    % one continuous time series
    tic
    % PT4 (monkey D)
    load('C:\Users\kkhateeb\Documents\MATLAB\PT_Neurophys\PT4_NeurophysiologyData\FilteredSignalsPT4.mat');

    lengthD = length(baselineBPF60_150.ipsi{1,1}) + ...
        length(duringPTBPF60_150.ipsi{1,1}) + ...
        length(postPT1BPF60_150.ipsi{1,1}) + ...
        length(postPT2BPF60_150.ipsi{1,1}) + length(postPTBPF60_150.ipsi{1,1});
    lengthD_base = length(baselineBPF60_150.ipsi{1,1});
    lengthD_drng = length(duringPTBPF60_150.ipsi{1,1});
    lengthD_post = length(postPT1BPF60_150.ipsi{1,1}) + ...
        length(postPT2BPF60_150.ipsi{1,1}) + length(postPTBPF60_150.ipsi{1,1});

    signalD_allbase_hg = zeros(lengthD_base,length(channels));
    signalD_allbase_lg = zeros(lengthD_base,length(channels));
    signalD_allbase_bt = zeros(lengthD_base,length(channels));
    signalD_allbase_th = zeros(lengthD_base,length(channels));

    signalD_drngill_hg = zeros(lengthD_drng,length(channels));
    signalD_drngill_lg = zeros(lengthD_drng,length(channels));
    signalD_drngill_bt = zeros(lengthD_drng,length(channels));
    signalD_drngill_th = zeros(lengthD_drng,length(channels));

    signalD_allpost_hg = zeros(lengthD_post,length(channels));
    signalD_allpost_lg = zeros(lengthD_post,length(channels));
    signalD_allpost_bt = zeros(lengthD_post,length(channels));
    signalD_allpost_th = zeros(lengthD_post,length(channels));
    % 
    for iChan = channels

        disp(['Channel ' num2str(iChan)])

        % All signals split into frequency bands
        signalD_allbase_hg(:,iChan) = [baselineBPF60_150.ipsi{1,iChan}]; % high gamma signal
        signalD_allbase_lg(:,iChan) = [baselineBPF30_59.ipsi{1,iChan}]; % low gamma signal
        signalD_allbase_bt(:,iChan) = [baselineBPF12_29.ipsi{1,iChan}]; % beta signal
        signalD_allbase_th(:,iChan) = [baselineBPF4_7.ipsi{1,iChan}]; % theta signal

        signalD_drngill_hg(:,iChan) = [duringPTBPF60_150.ipsi{1,iChan}]; % high gamma signal
        signalD_drngill_lg(:,iChan) = [duringPTBPF30_59.ipsi{1,iChan}]; % low gamma signal
        signalD_drngill_bt(:,iChan) = [duringPTBPF12_29.ipsi{1,iChan}]; % beta signal
        signalD_drngill_th(:,iChan) = [duringPTBPF4_7.ipsi{1,iChan}]; % theta signal

        signalD_allpost_hg(:,iChan) = [postPT1BPF60_150.ipsi{1,iChan};...
            postPT2BPF60_150.ipsi{1,iChan};postPTBPF60_150.ipsi{1,iChan}]; % high gamma signal
        signalD_allpost_lg(:,iChan) = [postPT1BPF30_59.ipsi{1,iChan};...
            postPT2BPF30_59.ipsi{1,iChan};postPTBPF30_59.ipsi{1,iChan}]; % low gamma signal
        signalD_allpost_bt(:,iChan) = [postPT1BPF12_29.ipsi{1,iChan};...
            postPT2BPF12_29.ipsi{1,iChan};postPTBPF12_29.ipsi{1,iChan}]; % beta signal
        signalD_allpost_th(:,iChan) = [postPT1BPF4_7.ipsi{1,iChan};...
            postPT2BPF4_7.ipsi{1,iChan};postPTBPF4_7.ipsi{1,iChan}]; % theta signal
        toc
    end
    % PT4 baseline artifacts same as before, no new ones ID'd
    disp('Loaded PT4 Signals')

    %% 10. Clean up all Signals for PT4 - Monkey D

    % as before, artifact indices determined using findArtifacts function
    noiseD_allbase = noiseD_base; % no change in baseline noise

    noiseD_drngill = [1, 1700;
        1583000, 1587500;
        1382000, 1383000;
        1384500, 1385550];

    noiseD_allpost = [5250320, 5250720;
        1, 1000;
        1618550, 1619500;
        2773250, 2773525;
        4101750, 4103250;
        4768500, 4769500;
        1550500, 1551250;];

    cleanD_allbase_hg = signalD_allbase_hg;
    cleanD_allbase_lg = signalD_allbase_lg;
    cleanD_allbase_bt = signalD_allbase_bt;
    cleanD_allbase_th = signalD_allbase_th;

    cleanD_drngill_hg = signalD_drngill_hg;
    cleanD_drngill_lg = signalD_drngill_lg;
    cleanD_drngill_bt = signalD_drngill_bt;
    cleanD_drngill_th = signalD_drngill_th;

    cleanD_allpost_hg = signalD_allpost_hg;
    cleanD_allpost_lg = signalD_allpost_lg;
    cleanD_allpost_bt = signalD_allpost_bt;
    cleanD_allpost_th = signalD_allpost_th;

    % replace noisey data with NaN
    for i = 1:size(noiseD_allbase,1)
        range = noiseD_allbase(i,1):noiseD_allbase(i,2);
        numSamps = 1+noiseD_allbase(i,2)-noiseD_allbase(i,1);
        cleanD_allbase_hg(range,:) = NaN*ones(numSamps,length(channels));
        cleanD_allbase_lg(range,:) = NaN*ones(numSamps,length(channels));
        cleanD_allbase_bt(range,:) = NaN*ones(numSamps,length(channels));
        cleanD_allbase_th(range,:) = NaN*ones(numSamps,length(channels));

    end

    for i = 1:size(noiseD_drngill,1)
        range = noiseD_drngill(i,1):noiseD_drngill(i,2);
        numSamps = 1+noiseD_drngill(i,2)-noiseD_drngill(i,1);
        cleanD_drngill_hg(range,:) = NaN*ones(numSamps,length(channels));
        cleanD_drngill_lg(range,:) = NaN*ones(numSamps,length(channels));
        cleanD_drngill_bt(range,:) = NaN*ones(numSamps,length(channels));
        cleanD_drngill_th(range,:) = NaN*ones(numSamps,length(channels));
    end

    for i = 1:size(noiseD_allpost,1)
        range = noiseD_allpost(i,1):noiseD_allpost(i,2);
        numSamps = 1+noiseD_allpost(i,2)-noiseD_allpost(i,1);
        cleanD_allpost_hg(range,:) = NaN*ones(numSamps,length(channels));
        cleanD_allpost_lg(range,:) = NaN*ones(numSamps,length(channels));
        cleanD_allpost_bt(range,:) = NaN*ones(numSamps,length(channels));
        cleanD_allpost_th(range,:) = NaN*ones(numSamps,length(channels));
    end

    % save signals for later use
    save('CleanedSignalsPT4.mat','cleanD_allbase_hg','cleanD_allbase_lg',...
        'cleanD_allbase_bt','cleanD_allbase_th','cleanD_allpost_hg','cleanD_allpost_lg',...
        'cleanD_allpost_bt','cleanD_allpost_th','cleanD_drngill_hg','cleanD_drngill_lg',...
        'cleanD_drngill_bt','cleanD_drngill_th', '-v7.3');

    disp('Cleaned up PT4 Signal')

    %% 11. Load all Filtered Signals for PT5

    % PT5 (monkey E)
    load('C:\Users\kkhateeb\Documents\MATLAB\PT_Neurophys\PT5_NeurophysiologyData\FilteredSignalsPT5.mat');

    lengthE = length(baselineBPF60_150.ipsi{1,1}) + ...
        length(duringPTBPF60_150.ipsi{1,1}) + length(postPTBPF60_150.ipsi{1,1});
    lengthE_base = length(baselineBPF60_150.ipsi{1,1});
    lengthE_drng = length(duringPTBPF60_150.ipsi{1,1});
    lengthE_post = length(postPTBPF60_150.ipsi{1,1});

    signalE_allbase_hg = zeros(lengthE_base,length(channels));
    signalE_allbase_lg = zeros(lengthE_base,length(channels));
    signalE_allbase_bt = zeros(lengthE_base,length(channels));
    signalE_allbase_th = zeros(lengthE_base,length(channels));

    signalE_drngill_hg = zeros(lengthE_drng,length(channels));
    signalE_drngill_lg = zeros(lengthE_drng,length(channels));
    signalE_drngill_bt = zeros(lengthE_drng,length(channels));
    signalE_drngill_th = zeros(lengthE_drng,length(channels));

    signalE_allpost_hg = zeros(lengthE_post,length(channels));
    signalE_allpost_lg = zeros(lengthE_post,length(channels));
    signalE_allpost_bt = zeros(lengthE_post,length(channels));
    signalE_allpost_th = zeros(lengthE_post,length(channels));
    tic
    for iChan = channels

        disp(['Channel ' num2str(iChan)])

        % All signals split into frequency bands
        signalE_allbase_hg(:,iChan) = [baselineBPF60_150.ipsi{1,iChan}]; % high gamma signal
        signalE_allbase_lg(:,iChan) = [baselineBPF30_59.ipsi{1,iChan}]; % low gamma signal
        signalE_allbase_bt(:,iChan) = [baselineBPF12_29.ipsi{1,iChan}]; % beta signal
        signalE_allbase_th(:,iChan) = [baselineBPF4_7.ipsi{1,iChan}]; % theta signal

        signalE_drngill_hg(:,iChan) = [duringPTBPF60_150.ipsi{1,iChan}]; % high gamma signal
        signalE_drngill_lg(:,iChan) = [duringPTBPF30_59.ipsi{1,iChan}]; % low gamma signal
        signalE_drngill_bt(:,iChan) = [duringPTBPF12_29.ipsi{1,iChan}]; % beta signal
        signalE_drngill_th(:,iChan) = [duringPTBPF4_7.ipsi{1,iChan}]; % theta signal

        signalE_allpost_hg(:,iChan) = [postPTBPF60_150.ipsi{1,iChan}]; % high gamma signal
        signalE_allpost_lg(:,iChan) = [postPTBPF30_59.ipsi{1,iChan}]; % low gamma signal
        signalE_allpost_bt(:,iChan) = [postPTBPF12_29.ipsi{1,iChan}]; % beta signal
        signalE_allpost_th(:,iChan) = [postPTBPF4_7.ipsi{1,iChan}]; % theta signal
        toc
    end

    disp('Loaded PT5 Signals')

    clearvars -REGEXP baseline* postPT* during*

    %% 12. Clean up all Signals for PT5 - Monkey E

    % as before, artifact indices determined using findArtifacts function
    noiseE_allbase = noiseE_base;

    noiseE_drngill = [1, 1160;
        107200, 109100;
        110400, 113800;
        127600, 131700;
        550300, 556700;
        1491000, lengthE_drng];

    noiseE_allpost = [6620, 7814;
        46700, 47430;
        53220, 57870;
        1088000, 1092700;
        1094000, 1099000;
        1106000, 1116500;
        1111000, 1111125;
        1121000, 1128000;
        1136000, 1138000;
        1143800, 1146000;
        1157000, 1158700;
        1164000, 1166000;
        1689750, 1690600;
        1684700, 1685250;
        1700000, 1702500;
        1123500, 1124500;
        2959000, 2860500;
        3545000, 3546500;                                
        2977500, 2978500;
        2519000, 2520000;
        2520500, 2521000;
        2426000, 2427000;
        2445000, 2446500;
        2447100, 2447500;
        2459000, 2460500;
        2529300, 2530000;
        2959500, 2960300;
        3868000, 3871000;
        3873000, 3874000;
        3874700, 3875300;
        3876300, 3879700;
        3882000, 3886500;
        3898200, 3898500;
        3959500, 3961500;
        3960400, 3960700;
        3966000, 3975000;
        4067600, 4068000;
        5413900, 5415300;
        5418800, 5419300;
        6499000, 6503500;
        6764500, 6767500;
        7817000, 7849000;
        7823000, 7823300;
        7873000, 7879000;
        7966000, 7970000;
        7993000, 7998000;
        8005000, 8011000;
        8294000, 8297000;
        8304000, 8309000;
        8341250, 8342000;
        8953000, 8958000;
        9008500, 9011500;
        9101000, 9110000;
        9112000, 9116000;
        9120000, 9129000;
        9132000, 9136000;
        9147000, 9151000;
        9155000, 9159000;
        lengthE_post-limit30+noiseE_post];

    cleanE_allbase_hg = signalE_allbase_hg;
    cleanE_allbase_lg = signalE_allbase_lg;
    cleanE_allbase_bt = signalE_allbase_bt;
    cleanE_allbase_th = signalE_allbase_th;

    cleanE_drngill_hg = signalE_drngill_hg;
    cleanE_drngill_lg = signalE_drngill_lg;
    cleanE_drngill_bt = signalE_drngill_bt;
    cleanE_drngill_th = signalE_drngill_th;

    cleanE_allpost_hg = signalE_allpost_hg;
    cleanE_allpost_lg = signalE_allpost_lg;
    cleanE_allpost_bt = signalE_allpost_bt;
    cleanE_allpost_th = signalE_allpost_th;

    tic
    % replace noisey data with NaN
    for i = 1:size(noiseE_allbase,1)
        range = noiseE_allbase(i,1):noiseE_allbase(i,2);
        numSamps = 1+noiseE_allbase(i,2)-noiseE_allbase(i,1);
        cleanE_allbase_hg(range,:) = NaN*ones(numSamps,length(channels));
        cleanE_allbase_lg(range,:) = NaN*ones(numSamps,length(channels));
        cleanE_allbase_bt(range,:) = NaN*ones(numSamps,length(channels));
        cleanE_allbase_th(range,:) = NaN*ones(numSamps,length(channels));
        toc
    end

    for i = 1:size(noiseE_drngill,1)
        range = noiseE_drngill(i,1):noiseE_drngill(i,2);
        numSamps = 1+noiseE_drngill(i,2)-noiseE_drngill(i,1);
        cleanE_drngill_hg(range,:) = NaN*ones(numSamps,length(channels));
        cleanE_drngill_lg(range,:) = NaN*ones(numSamps,length(channels));
        cleanE_drngill_bt(range,:) = NaN*ones(numSamps,length(channels));
        cleanE_drngill_th(range,:) = NaN*ones(numSamps,length(channels));
        toc
    end

    for i = 1:size(noiseE_allpost,1)
        range = noiseE_allpost(i,1):noiseE_allpost(i,2);
        numSamps = 1+noiseE_allpost(i,2)-noiseE_allpost(i,1);
        cleanE_allpost_hg(range,:) = NaN*ones(numSamps,length(channels));
        cleanE_allpost_lg(range,:) = NaN*ones(numSamps,length(channels));
        cleanE_allpost_bt(range,:) = NaN*ones(numSamps,length(channels));
        cleanE_allpost_th(range,:) = NaN*ones(numSamps,length(channels));
        toc
    end

    % save cleaned signals for later use
    save('CleanedSignalsPT5.mat','cleanE_allbase_hg','cleanE_allbase_lg',...
        'cleanE_allbase_bt','cleanE_allbase_th','cleanE_allpost_hg','cleanE_allpost_lg',...
        'cleanE_allpost_bt','cleanE_allpost_th','cleanE_drngill_hg','cleanE_drngill_lg',...
        'cleanE_drngill_bt','cleanE_drngill_th', '-v7.3', '-append');

    disp('Cleaned up PT5 Signal')
    
    % save artifact indices for later use
    save('AllArtifactIndices.mat','noiseD_allbase','noiseD_drngill',...
        'noiseD_allpost','noiseE_allbase','noiseE_drngill','noiseE_allpost');

    clearvars -REGEXP signalD_all* signalD_drngill* signalE_all* signalE_drngill*
end
%% 13. Load Monkey D All Time-Series

% if not already loaded, load the cleaned up signals
if ~getAllSigs
    tic
    load('CleanedSignalsPT4.mat');
    disp('Loading Cleaned PT4 Signals')
    toc
end

%% 14. Calculate Power (dt = 10s) for Monkey D All Time-Series

dt = 10*Fs; % number of samples every 10 seconds

lngthD_allbase = size(cleanD_allbase_hg,1);
lngthD_drngill = size(cleanD_drngill_hg,1);
lngthD_allpost = size(cleanD_allpost_hg,1);

NdtD_allbase = floor(lngthD_allbase/dt); % number of 10 second intervals
NdtD_drngill = floor(lngthD_drngill/dt);
NdtD_allpost = floor(lngthD_allpost/dt);

pwrD_allbase_hg = zeros(NdtD_allbase,length(channels));
pwrD_allbase_lg = zeros(NdtD_allbase,length(channels));
pwrD_allbase_bt = zeros(NdtD_allbase,length(channels));
pwrD_allbase_th = zeros(NdtD_allbase,length(channels));

pwrD_drngill_hg = zeros(NdtD_drngill,length(channels));
pwrD_drngill_lg = zeros(NdtD_drngill,length(channels));
pwrD_drngill_bt = zeros(NdtD_drngill,length(channels));
pwrD_drngill_th = zeros(NdtD_drngill,length(channels));

pwrD_allpost_hg = zeros(NdtD_allpost,length(channels));
pwrD_allpost_lg = zeros(NdtD_allpost,length(channels));
pwrD_allpost_bt = zeros(NdtD_allpost,length(channels));
pwrD_allpost_th = zeros(NdtD_allpost,length(channels));

tic
for idt = 1:NdtD_allbase
    
    interval = ((idt-1)*dt+1):(idt*dt); % time interval in samples
    
    for iChan = channels
        if iChan ~= 8
            pwrD_allbase_hg(idt,iChan) = SignalPower(cleanD_allbase_hg(interval,iChan),Fs);
            pwrD_allbase_lg(idt,iChan) = SignalPower(cleanD_allbase_lg(interval,iChan),Fs);       
            pwrD_allbase_bt(idt,iChan) = SignalPower(cleanD_allbase_bt(interval,iChan),Fs);       
            pwrD_allbase_th(idt,iChan) = SignalPower(cleanD_allbase_th(interval,iChan),Fs);
        end        
    end
end
toc
disp('Calculated All Baseline Power')

for idt = 1:NdtD_drngill
    
    interval = ((idt-1)*dt+1):(idt*dt); % time interval in samples
    
    for iChan = channels
        if iChan ~= 8
            pwrD_drngill_hg(idt,iChan) = SignalPower(cleanD_drngill_hg(interval,iChan),Fs);
            pwrD_drngill_lg(idt,iChan) = SignalPower(cleanD_drngill_lg(interval,iChan),Fs);       
            pwrD_drngill_bt(idt,iChan) = SignalPower(cleanD_drngill_bt(interval,iChan),Fs);       
            pwrD_drngill_th(idt,iChan) = SignalPower(cleanD_drngill_th(interval,iChan),Fs);
        end        
    end
end
toc
disp('Calculated All During Illumination Power')

for idt = 1:NdtD_allpost
    
    interval = ((idt-1)*dt+1):(idt*dt); % time interval in samples
    
    for iChan = channels
        if iChan ~= 8
            pwrD_allpost_hg(idt,iChan) = SignalPower(cleanD_allpost_hg(interval,iChan),Fs);
            pwrD_allpost_lg(idt,iChan) = SignalPower(cleanD_allpost_lg(interval,iChan),Fs);       
            pwrD_allpost_bt(idt,iChan) = SignalPower(cleanD_allpost_bt(interval,iChan),Fs);      
            pwrD_allpost_th(idt,iChan) = SignalPower(cleanD_allpost_th(interval,iChan),Fs);
        end        
    end
end
toc
disp('Calculated All PostPT Power')

% Remove channel 8 for Monkey D
pwrD_allbase_hg(:,8) = nan*ones(size(pwrD_allbase_hg(:,8))); 
pwrD_allbase_lg(:,8) = nan*ones(size(pwrD_allbase_lg(:,8)));
pwrD_allbase_bt(:,8) = nan*ones(size(pwrD_allbase_bt(:,8)));
pwrD_allbase_th(:,8) = nan*ones(size(pwrD_allbase_th(:,8)));

pwrD_drngill_hg(:,8) = nan*ones(size(pwrD_drngill_hg(:,8))); 
pwrD_drngill_lg(:,8) = nan*ones(size(pwrD_drngill_lg(:,8)));
pwrD_drngill_bt(:,8) = nan*ones(size(pwrD_drngill_bt(:,8)));
pwrD_drngill_th(:,8) = nan*ones(size(pwrD_drngill_th(:,8)));

pwrD_allpost_hg(:,8) = nan*ones(size(pwrD_allpost_hg(:,8)));
pwrD_allpost_lg(:,8) = nan*ones(size(pwrD_allpost_lg(:,8)));
pwrD_allpost_bt(:,8) = nan*ones(size(pwrD_allpost_bt(:,8)));
pwrD_allpost_th(:,8) = nan*ones(size(pwrD_allpost_th(:,8)));

disp('Finished Calculated Power')

pwrD_all_hg = [pwrD_allbase_hg; pwrD_drngill_hg; pwrD_allpost_hg];
pwrD_all_lg = [pwrD_allbase_lg; pwrD_drngill_lg; pwrD_allpost_lg];
pwrD_all_bt = [pwrD_allbase_bt; pwrD_drngill_bt; pwrD_allpost_bt];
pwrD_all_th = [pwrD_allbase_th; pwrD_drngill_th; pwrD_allpost_th];

%% 15. Load Monkey E All Time-Series

% if not loaded, load the cleaned up signals
if ~getAllSigs
    tic
    load('CleanedSignalsPT5.mat');
    disp('Loading Cleaned PT5 Signals')
    toc
end

%% 15. Calculate Power (dt = 10s) for Monkey E All Time-Series

dt = 10*Fs; % number of samples every 10 seconds

lngthE_allbase = size(cleanE_allbase_hg,1);
lngthE_drngill = size(cleanE_drngill_hg,1);
lngthE_allpost = size(cleanE_allpost_hg,1);

NdtE_allbase = floor(lngthE_allbase/dt); % number of 10 second intervals
NdtE_drngill = floor(lngthE_drngill/dt);
NdtE_allpost = floor(lngthE_allpost/dt);

pwrE_allbase_hg = zeros(NdtE_allbase,length(channels));
pwrE_allbase_lg = zeros(NdtE_allbase,length(channels));
pwrE_allbase_bt = zeros(NdtE_allbase,length(channels));
pwrE_allbase_th = zeros(NdtE_allbase,length(channels));

pwrE_drngill_hg = zeros(NdtE_drngill,length(channels));
pwrE_drngill_lg = zeros(NdtE_drngill,length(channels));
pwrE_drngill_bt = zeros(NdtE_drngill,length(channels));
pwrE_drngill_th = zeros(NdtE_drngill,length(channels));

pwrE_allpost_hg = zeros(NdtE_allpost,length(channels));
pwrE_allpost_lg = zeros(NdtE_allpost,length(channels));
pwrE_allpost_bt = zeros(NdtE_allpost,length(channels));
pwrE_allpost_th = zeros(NdtE_allpost,length(channels));

tic
for idt = 1:NdtE_allbase
    
    interval = ((idt-1)*dt+1):(idt*dt); % time interval in samples
    
    for iChan = channels
        if iChan ~= 8
            pwrE_allbase_hg(idt,iChan) = SignalPower(cleanE_allbase_hg(interval,iChan),Fs);
            pwrE_allbase_lg(idt,iChan) = SignalPower(cleanE_allbase_lg(interval,iChan),Fs);       
            pwrE_allbase_bt(idt,iChan) = SignalPower(cleanE_allbase_bt(interval,iChan),Fs);       
            pwrE_allbase_th(idt,iChan) = SignalPower(cleanE_allbase_th(interval,iChan),Fs);
        end        
    end
end
toc
disp('Calculated All Baseline Power')

for idt = 1:NdtE_drngill
    
    interval = ((idt-1)*dt+1):(idt*dt); % time interval in samples
    
    for iChan = channels
        if iChan ~= 8
            pwrE_drngill_hg(idt,iChan) = SignalPower(cleanE_drngill_hg(interval,iChan),Fs);
            pwrE_drngill_lg(idt,iChan) = SignalPower(cleanE_drngill_lg(interval,iChan),Fs);       
            pwrE_drngill_bt(idt,iChan) = SignalPower(cleanE_drngill_bt(interval,iChan),Fs);       
            pwrE_drngill_th(idt,iChan) = SignalPower(cleanE_drngill_th(interval,iChan),Fs);
        end        
    end
end
toc
disp('Calculated All During Illumination Power')

for idt = 1:NdtE_allpost
    
    interval = ((idt-1)*dt+1):(idt*dt); % time interval in samples
    
    for iChan = channels
        if iChan ~= 8 && iChan ~= 24 && iChan ~= 25
            pwrE_allpost_hg(idt,iChan) = SignalPower(cleanE_allpost_hg(interval,iChan),Fs);
            pwrE_allpost_lg(idt,iChan) = SignalPower(cleanE_allpost_lg(interval,iChan),Fs);       
            pwrE_allpost_bt(idt,iChan) = SignalPower(cleanE_allpost_bt(interval,iChan),Fs);      
            pwrE_allpost_th(idt,iChan) = SignalPower(cleanE_allpost_th(interval,iChan),Fs);
                    
        end
    end
end
toc
disp('Calculated All PostPT Power')

% Remove channel 8, 24, and 25 for Monkey E
for iChan = [8, 24, 25]
    pwrE_allbase_hg(:,iChan) = nan*ones(size(pwrE_allbase_hg(:,8))); 
    pwrE_allbase_lg(:,iChan) = nan*ones(size(pwrE_allbase_lg(:,8)));
    pwrE_allbase_bt(:,iChan) = nan*ones(size(pwrE_allbase_bt(:,8)));
    pwrE_allbase_th(:,iChan) = nan*ones(size(pwrE_allbase_th(:,8)));

    pwrE_drngill_hg(:,iChan) = nan*ones(size(pwrE_drngill_hg(:,8))); 
    pwrE_drngill_lg(:,iChan) = nan*ones(size(pwrE_drngill_lg(:,8)));
    pwrE_drngill_bt(:,iChan) = nan*ones(size(pwrE_drngill_bt(:,8)));
    pwrE_drngill_th(:,iChan) = nan*ones(size(pwrE_drngill_th(:,8)));

    pwrE_allpost_hg(:,iChan) = nan*ones(size(pwrE_allpost_hg(:,8)));
    pwrE_allpost_lg(:,iChan) = nan*ones(size(pwrE_allpost_lg(:,8)));
    pwrE_allpost_bt(:,iChan) = nan*ones(size(pwrE_allpost_bt(:,8)));
    pwrE_allpost_th(:,iChan) = nan*ones(size(pwrE_allpost_th(:,8)));
end

disp('Finished Calculated Power')

pwrE_all_hg = [pwrE_allbase_hg; pwrE_drngill_hg; pwrE_allpost_hg];
pwrE_all_lg = [pwrE_allbase_lg; pwrE_drngill_lg; pwrE_allpost_lg];
pwrE_all_bt = [pwrE_allbase_bt; pwrE_drngill_bt; pwrE_allpost_bt];
pwrE_all_th = [pwrE_allbase_th; pwrE_drngill_th; pwrE_allpost_th];

% Save power time course for potential future use
save('PowerTimeCourse.mat',...
    'pwrD_allbase_hg','pwrD_drngill_hg','pwrD_allpost_hg',...
    'pwrD_allbase_lg','pwrD_drngill_lg','pwrD_allpost_lg',...
    'pwrD_allbase_bt','pwrD_drngill_bt','pwrD_allpost_bt',...
    'pwrD_allbase_th','pwrD_drngill_th','pwrD_allpost_th',...
    'pwrD_all_hg','pwrD_all_lg','pwrD_all_bt','pwrD_all_th',...
    'pwrE_allbase_hg','pwrE_drngill_hg','pwrE_allpost_hg',...
    'pwrE_allbase_lg','pwrE_drngill_lg','pwrE_allpost_lg',...
    'pwrE_allbase_bt','pwrE_drngill_bt','pwrE_allpost_bt',...
    'pwrE_allbase_th','pwrE_drngill_th','pwrE_allpost_th',...
    'pwrE_all_hg','pwrE_all_lg','pwrE_all_bt','pwrE_all_th');
    
%% 16. Smooth Power vs Time for Both Monkeys

smoFactor = 0.4; % value from 0 to 1 (0.25 default)

smPwrD_all_hg = smoothdata(pwrD_all_hg,'SmoothingFactor',smoFactor);
smPwrD_all_lg = smoothdata(pwrD_all_lg,'SmoothingFactor',smoFactor);
smPwrD_all_bt = smoothdata(pwrD_all_bt,'SmoothingFactor',smoFactor);
smPwrD_all_th = smoothdata(pwrD_all_th,'SmoothingFactor',smoFactor);

smoFactor = 0.6; % value from 0 to 1 (0.25 default)
smPwrE_all_hg = smoothdata(pwrE_all_hg,'SmoothingFactor',smoFactor);
smPwrE_all_lg = smoothdata(pwrE_all_lg,'SmoothingFactor',smoFactor);
smPwrE_all_bt = smoothdata(pwrE_all_bt,'SmoothingFactor',smoFactor);
smPwrE_all_th = smoothdata(pwrE_all_th,'SmoothingFactor',smoFactor);

%% 17. Make Nice Heat Map Figures for the Manuscript

% Get the same color scale for each frequency band
min_hg = min([totPwrD_base_hg,totPwrD_post_hg]);
max_hg = max([totPwrD_base_hg,totPwrD_post_hg]);

min_lg = min([totPwrD_base_lg,totPwrD_post_lg]);
max_lg = max([totPwrD_base_lg,totPwrD_post_lg]);

min_bt = min([totPwrD_base_bt,totPwrD_post_bt]);
max_bt = max([totPwrD_base_bt,totPwrD_post_bt]);

min_th = min([totPwrD_base_th,totPwrD_post_th]);
max_th = max([totPwrD_base_th,totPwrD_post_th]);

% set x,y,w,h for subplots of array heat maps
xPos = [.04 .35 .65]; yPos = [.76 .52 .28 .04]; w = .2; h = .2;

% Show Baseline, PostpT, and Change in Power Heat Maps for Monkey D only
figure('color','w','units','normalized','outerposition',[.3 0 .3 1])
subplot('position',[xPos(1) yPos(1) w h]), mapArray(totPwrD_base_hg), ...
    caxis([min_hg,max_hg]), title('Baseline Power','fontsize',14), ...
    ylabel('High Gamma','fontsize',14)
subplot('position',[xPos(2) yPos(1) w h]), mapArray(totPwrD_post_hg), ...
    caxis([min_hg,max_hg]), title('Post-PT Power','fontsize',14), ...
    axis off, cbarset(1,1,xPos(1),yPos(1));
subplot('position',[xPos(3) yPos(1) w h]), mapArray(subTotD_hg), ...
    axis off, title('Power Change','fontsize',14), ...
    cbarset(2,1,xPos(3),yPos(1));

subplot('position',[xPos(1) yPos(2) w h]), mapArray(totPwrD_base_lg), ...
    caxis([min_lg,max_lg]), ylabel('Low Gamma','fontsize',14)
subplot('position',[xPos(2) yPos(2) w h]), mapArray(totPwrD_post_lg), ...
    caxis([min_lg,max_lg]), axis off, cbarset(1,2,xPos(1),yPos(2));
subplot('position',[xPos(3) yPos(2) w h]), mapArray(subTotD_lg), ...
    axis off, cbarset(2,2,xPos(3),yPos(2));

subplot('position',[xPos(1) yPos(3) w h]), mapArray(totPwrD_base_bt), ...
    caxis([min_bt,max_bt]), ylabel('Beta','fontsize',14)
subplot('position',[xPos(2) yPos(3) w h]), mapArray(totPwrD_post_bt), ...
    caxis([min_bt,max_bt]), axis off, cbarset(1,3,xPos(1),yPos(3));
subplot('position',[xPos(3) yPos(3) w h]), mapArray(subTotD_bt), ...
    axis off, cbarset(2,3,xPos(3),yPos(3));

subplot('position',[xPos(1) yPos(4) w h]), mapArray(totPwrD_base_th), ...
    caxis([min_th,max_th]), ylabel('Theta','fontsize',14)
subplot('position',[xPos(2) yPos(4) w h]), mapArray(totPwrD_post_th), ...
    caxis([min_th,max_th]), axis off, cbarset(1,4,xPos(1),yPos(4));
subplot('position',[xPos(3) yPos(4) w h]), mapArray(subTotD_th), ...
    axis off, cbarset(2,4,xPos(3),yPos(4));

[~,ht] = suplabel('Changes Network Power','t',[.08 .147 .84 .84]);
set(ht,'fontweight','bold','fontsize',18);


%% 19. Group Channels into Lesion and Non-Lesioned Groups

% Lesioned group: post shows significant difference from baseline
% Non-lesioned group: failed to show a difference

dt = 60*Fs; % number of samples every N seconds * Fs
Ndt = limit30/dt; % number of dt/Fs second intervals in 30 minutes

dtpwrD_base_hg = zeros(Ndt,length(channels));
dtpwrD_base_lg = zeros(Ndt,length(channels));
dtpwrD_base_bt = zeros(Ndt,length(channels));
dtpwrD_base_th = zeros(Ndt,length(channels));
dtpwrD_post_hg = zeros(Ndt,length(channels));
dtpwrD_post_lg = zeros(Ndt,length(channels));
dtpwrD_post_bt = zeros(Ndt,length(channels));
dtpwrD_post_th = zeros(Ndt,length(channels));

dtpwrE_base_hg = zeros(Ndt,length(channels));
dtpwrE_base_lg = zeros(Ndt,length(channels));
dtpwrE_base_bt = zeros(Ndt,length(channels));
dtpwrE_base_th = zeros(Ndt,length(channels));
dtpwrE_post_hg = zeros(Ndt,length(channels));
dtpwrE_post_lg = zeros(Ndt,length(channels));
dtpwrE_post_bt = zeros(Ndt,length(channels));
dtpwrE_post_th = zeros(Ndt,length(channels));

for idt = 1:Ndt
    
    interval = ((idt-1)*dt+1):(idt*dt); % time interval in samples
    
    for iChan = channels
        dtpwrD_base_hg(idt,iChan) = SignalPower(cleanD_base_hg(interval,iChan),Fs);
        dtpwrD_base_lg(idt,iChan) = SignalPower(cleanD_base_lg(interval,iChan),Fs);       
        dtpwrD_base_bt(idt,iChan) = SignalPower(cleanD_base_bt(interval,iChan),Fs);       
        dtpwrD_base_th(idt,iChan) = SignalPower(cleanD_base_th(interval,iChan),Fs);
        dtpwrD_post_hg(idt,iChan) = SignalPower(cleanD_post_hg(interval,iChan),Fs);
        dtpwrD_post_lg(idt,iChan) = SignalPower(cleanD_post_lg(interval,iChan),Fs);       
        dtpwrD_post_bt(idt,iChan) = SignalPower(cleanD_post_bt(interval,iChan),Fs);      
        dtpwrD_post_th(idt,iChan) = SignalPower(cleanD_post_th(interval,iChan),Fs);
        
        dtpwrE_base_hg(idt,iChan) = SignalPower(cleanE_base_hg(interval,iChan),Fs);
        dtpwrE_base_lg(idt,iChan) = SignalPower(cleanE_base_lg(interval,iChan),Fs);       
        dtpwrE_base_bt(idt,iChan) = SignalPower(cleanE_base_bt(interval,iChan),Fs);       
        dtpwrE_base_th(idt,iChan) = SignalPower(cleanE_base_th(interval,iChan),Fs);
        dtpwrE_post_hg(idt,iChan) = SignalPower(cleanE_post_hg(interval,iChan),Fs);
        dtpwrE_post_lg(idt,iChan) = SignalPower(cleanE_post_lg(interval,iChan),Fs);       
        dtpwrE_post_bt(idt,iChan) = SignalPower(cleanE_post_bt(interval,iChan),Fs);       
        dtpwrE_post_th(idt,iChan) = SignalPower(cleanE_post_th(interval,iChan),Fs);
    end
end

% Remove channel 8 for Monkey D and channels 8, 24 and 25 for monkey E
dtpwrD_base_hg(:,8) = nan*ones(size(dtpwrD_base_hg(:,8))); 
dtpwrD_base_lg(:,8) = nan*ones(size(dtpwrD_base_lg(:,8)));
dtpwrD_base_bt(:,8) = nan*ones(size(dtpwrD_base_bt(:,8)));
dtpwrD_base_th(:,8) = nan*ones(size(dtpwrD_base_th(:,8)));
dtpwrD_post_hg(:,8) = nan*ones(size(dtpwrD_post_hg(:,8)));
dtpwrD_post_lg(:,8) = nan*ones(size(dtpwrD_post_lg(:,8)));
dtpwrD_post_bt(:,8) = nan*ones(size(dtpwrD_post_bt(:,8)));
dtpwrD_post_th(:,8) = nan*ones(size(dtpwrD_post_th(:,8)));

dtpwrE_base_hg(:,8) = nan*ones(size(dtpwrE_base_hg(:,8))); 
dtpwrE_base_lg(:,8) = nan*ones(size(dtpwrE_base_lg(:,8)));
dtpwrE_base_bt(:,8) = nan*ones(size(dtpwrE_base_bt(:,8)));
dtpwrE_base_th(:,8) = nan*ones(size(dtpwrE_base_th(:,8)));
dtpwrE_post_hg(:,8) = nan*ones(size(dtpwrE_post_hg(:,8)));
dtpwrE_post_lg(:,8) = nan*ones(size(dtpwrE_post_lg(:,8)));
dtpwrE_post_bt(:,8) = nan*ones(size(dtpwrE_post_bt(:,8)));
dtpwrE_post_th(:,8) = nan*ones(size(dtpwrE_post_th(:,8)));

dtpwrE_base_hg(:,24) = nan*ones(size(dtpwrE_base_hg(:,24))); 
dtpwrE_base_lg(:,24) = nan*ones(size(dtpwrE_base_lg(:,24)));
dtpwrE_base_bt(:,24) = nan*ones(size(dtpwrE_base_bt(:,24)));
dtpwrE_base_th(:,24) = nan*ones(size(dtpwrE_base_th(:,24)));
dtpwrE_post_hg(:,24) = nan*ones(size(dtpwrE_post_hg(:,24)));
dtpwrE_post_lg(:,24) = nan*ones(size(dtpwrE_post_lg(:,24)));
dtpwrE_post_bt(:,24) = nan*ones(size(dtpwrE_post_bt(:,24)));
dtpwrE_post_th(:,24) = nan*ones(size(dtpwrE_post_th(:,24)));

dtpwrE_base_hg(:,25) = nan*ones(size(dtpwrE_base_hg(:,25))); 
dtpwrE_base_lg(:,25) = nan*ones(size(dtpwrE_base_lg(:,25)));
dtpwrE_base_bt(:,25) = nan*ones(size(dtpwrE_base_bt(:,25)));
dtpwrE_base_th(:,25) = nan*ones(size(dtpwrE_base_th(:,25)));
dtpwrE_post_hg(:,25) = nan*ones(size(dtpwrE_post_hg(:,25)));
dtpwrE_post_lg(:,25) = nan*ones(size(dtpwrE_post_lg(:,25)));
dtpwrE_post_bt(:,25) = nan*ones(size(dtpwrE_post_bt(:,25)));
dtpwrE_post_th(:,25) = nan*ones(size(dtpwrE_post_th(:,25)));

disp('Calculated Power')

m = length(channels)*4*2; % number of comparisons: numChan*freqBands*numAnimals
FWER = 0.001; % family-wise error rate to use
alpha = FWER/m;

[hlD_hg,plD_hg,sortD_hg] = SortChannelsNeg(dtpwrD_post_hg,dtpwrD_base_hg,alpha);
[hlD_lg,plD_lg,sortD_lg] = SortChannelsNeg(dtpwrD_post_lg,dtpwrD_base_lg,alpha);
[hlD_bt,plD_bt,sortD_bt] = SortChannelsNeg(dtpwrD_post_bt,dtpwrD_base_bt,alpha);
[hlD_th,plD_th,sortD_th] = SortChannelsNeg(dtpwrD_post_th,dtpwrD_base_th,alpha);

[hlE_hg,plE_hg,sortE_hg] = SortChannelsNeg(dtpwrE_post_hg,dtpwrE_base_hg,alpha);
[hlE_lg,plE_lg,sortE_lg] = SortChannelsNeg(dtpwrE_post_lg,dtpwrE_base_lg,alpha);
[hlE_bt,plE_bt,sortE_bt] = SortChannelsNeg(dtpwrE_post_bt,dtpwrE_base_bt,alpha);
[hlE_th,plE_th,sortE_th] = SortChannelsNeg(dtpwrE_post_th,dtpwrE_base_th,alpha);

disp('Sorted the Channels')

figure('color','w','units','normalized','outerposition',[.3 0 .4 1])
subplot(4,2,1), clustArray(sortD_hg,cmp), caxis([1 3]), title('Monkey D'), ylabel('High Gamma')
subplot(4,2,2), clustArray(sortE_hg,cmp), caxis([1 3]), title('Monkey E'), axis off
subplot(4,2,3), clustArray(sortD_lg,cmp), caxis([1 3]), ylabel('Low Gamma')
subplot(4,2,4), clustArray(sortE_lg,cmp), caxis([1 3]), axis off
subplot(4,2,5), clustArray(sortD_bt,cmp), caxis([1 3]), ylabel('Beta')
subplot(4,2,6), clustArray(sortE_bt,cmp), caxis([1 3]), axis off
subplot(4,2,7), clustArray(sortD_th,cmp), caxis([1 3]), ylabel('Theta')
subplot(4,2,8), clustArray(sortE_th,cmp), caxis([1 3]), axis off
[~,ht] = suplabel(['Sorted Ipsi Channels, FWER = ' num2str(FWER) ', dt = ' num2str(dt/Fs) ' seconds'],'t',[.08 .125 .84 .84]);
set(ht,'fontsize',18);

% Save channel grouping for later use
save(['Groupings_' num2str(dt/Fs) 'dt_' num2str(FWER) 'FWER.mat'],...
    'sortD_hg','sortD_lg','sortD_bt','sortD_th',...
    'sortE_hg','sortE_lg','sortE_bt','sortE_th');

% Save figure
save(['LowGammaIpsiGroupings_' num2str(dt/Fs) 'dt_' num2str(FWER) 'FWER.mat'],'sortD_lg','sortE_lg');
%% 20. Make Sorted Array Figure for Manuscript
min1 = min([min_lg minE_lg]); max1 = max([max_lg maxE_lg]);
min2 = min([subTotD_lg subTotE_lg],[],'all');
max2 = max([subTotD_lg subTotE_lg],[],'all');
MAX = max1; MIN = min2;
% MAX = max1; MIN = -max1;

figure('color','w','units','normalize','outerposition',[0 0 .4 1])
xPos = [0 .25 .5 .75]; yPos = [.6 .35]; w = .25; h = .25;
ax = subplot('position',[xPos(1) yPos(1) w h],'color','none'),...
    mapArray(totPwrD_base_lg);...
    title('Baseline','fontsize',14),...
    ylabel('Monkey D','fontsize',14),...
    caxis([MIN MAX]);
ax = subplot('position',[xPos(2) yPos(1) w h],'color','none'),...
    mapArray(totPwrD_post_lg);...
    title('Post-PT','fontsize',14),...
    caxis([MIN MAX]); axis off
ax = subplot('position',[xPos(3) yPos(1) w h],'color','none'),...
    mapArray(subTotD_lg);...
    title('Change','fontsize',14),...
    caxis([MIN MAX]); axis off
ax = subplot('position',[xPos(4) yPos(1) w h],'color','none'),...
    clustArray(sortD_lg,cmp(1:2,:),gca);...
    title('Sorted','fontsize',14); axis off, hold off

ax = subplot('position',[xPos(1) yPos(2) w h],'color','none'),...
    mapArray(totPwrE_base_lg,NaN,ax);...
    ylabel('Monkey E','fontsize',14),...
    caxis([MIN MAX]);
ax = subplot('position',[xPos(2) yPos(2) w h],'color','none'),...
    mapArray(totPwrE_post_lg,NaN,ax);...
    axis off,...
    caxis([MIN MAX]);
ax = subplot('position',[xPos(3) yPos(2) w h],'color','none'),...
    mapArray(subTotE_lg,NaN,ax);...
    axis off,...
    caxis([MIN MAX]);
c = colorbar('horiz'); c.TickDirection = 'in';...
    c.Position = [xPos(1)+.0816 yPos(2)-.07 2.34*w .025];...
    set(c, 'fontsize',14,'xtick',[-3e5 -2e5 -1e5 0 1e5 2e5 3e5],'TickLength',[.005 .025]); 
subplot('position',[xPos(4) yPos(2) w h],'color','none'),...
    clustArray(sortE_lg,cmp(1:2,:),gca);...
    axis off

print -dpdf -painters SortedPwrray

%% 21. Show Power vs Time for Sorted Channels

dt = 10*Fs;

figure('color','w','units','normalize','outerposition',[0 0 1 1])
xPos3 = 0.05; yPos3 = [.76 .54 .32 .1]; w = .9; h = .2;
subplot('position',[xPos3(1) yPos3(2) w h]),...
     rectangle('position',[NdtD_allbase 100 (NdtD_drngill) 5e5],'facecolor',[[255 255 92]/255 0.7],'edgecolor','none'); hold on,...
     clustPlot(smPwrD_all_lg,sortD_lg,cmp), box off, ...
     ylabel('Low Gamma','fontsize',14), ylim([0 5e5]),...
     xline(NdtD_allbase,'--','linewidth',3);...
     xline(NdtD_drngill+NdtD_allbase,'--','linewidth',3);...
     xline(NdtD_allpost+NdtD_drngill+NdtD_allbase,'--','linewidth',3);...
     xline(NdtD_allpost+NdtD_drngill+NdtD_allbase-(limit30/dt)+1,'--','linewidth',3);...
     xlim([0 NdtD_allpost+NdtD_drngill+NdtD_allbase]);...
     xticks([180 360 540 720 900]); xticklabels([]);
     xl = xlabel('Time (min)','fontsize',14);...
%      xl.Position(2) = xl.Position(2)+.03;...
     xticks([180 360 540 720 900]); xticklabels({30,60,90,120,180});

print -dpdf -painters sortedPwrTPlot

%% 22. Show bar graph of sorted low gamma groups

mColor = [0.75 0.75 0.75]; mAlpha = 0.5;
eColor = [0.5 0.5 0.5];

% seperate channels in each group
negGroupD = subTotD_lg(sortD_lg == 1);
regGroupD = subTotD_lg(sortD_lg == 2);

negGroupE = subTotE_lg(sortE_lg == 1);
regGroupE = subTotE_lg(sortE_lg == 2);

% store mean power change for each group into single array
groupsD = [mean(negGroupD) mean(regGroupD)];
groupsE = [mean(negGroupE) mean(regGroupE)];

% store the standard error of each group
stdGroupsD = [std(negGroupD)./sqrt(length(negGroupD)) std(regGroupD)./sqrt(length(regGroupD))];
stdGroupsE = [std(negGroupE)./sqrt(length(negGroupE)) std(regGroupE)./sqrt(length(regGroupE))];

groupedC(1,1,:) = cmp(1,:); groupedC(2,1,:) = cmp(2,:);
xpos = getGroupPos(groupsD(1,:)'); xpos = xpos([2 3]);
group = {[xpos(1),xpos(2)]};

% Do a one-sided sign test for groups 1 and 2
pGroupSTD1 = signtest(negGroupD,0,'tail','left'); % median clust1 negative?
pGroupSTD2 = signtest(regGroupD,0,'tail','left'); % median clust2 negative?

pGroupSTE1 = signtest(negGroupE,0,'tail','left'); % median clust1 negative?
pGroupSTE2 = signtest(regGroupE,0,'tail','left'); % median clust2 negative?

% compare negGroup and regGroup for each monkey
[p_groupD,t_groupD] = clustStats(negGroupD,regGroupD,[]);

[p_groupE,t_groupE] = clustStats(negGroupE,regGroupE,[]);

% don't use multCompare p value since only 2 comparisons
mpval_groupD = [p_groupD];
mpval_groupE = [p_groupE];

figure('color','w','units','normalize','outerposition',[0 0 1 1])
xPos2 = [.3 .425]+.03; yPos2 = [.55]; w = .1; h = .3;
subplot('position',[xPos2(1) yPos2(1) w h],'color','none'),...
    superbar(groupsD(1,:),'E', stdGroupsD(1,:),'BarFaceColor',groupedC,'PStarShowGT', false,...
    'PStarIcon',char(167),'PStarOffset',.2*abs(max(groupsD(1,:)))); box off, xlim([.5 2.5]),...
    set(gca,'xticklabel',{[]},'fontsize',14); H = sigstar(group,mpval_groupD,0,1); ...
    title('Monkey D','fontsize',14); ...
    set(H(:,2),'fontsize',14); ylim([-20e4 13e4]);...
    ylabel('Power Change (\muV^2/s)'), hold on, ...
    scatter(ones(size(negGroupD)),negGroupD,'filled','MarkerEdgeColor',eColor,'MarkerFaceColor',mColor,'MarkerFaceAlpha',mAlpha), ...
    scatter(2*ones(size(regGroupD)),regGroupD,'filled','MarkerEdgeColor',eColor,'MarkerFaceColor',mColor,'MarkerFaceAlpha',mAlpha)
subplot('position',[xPos2(2) yPos2(1) w h],'color','none'),...
    superbar(groupsE(1,:),'E', stdGroupsE(1,:),'BarFaceColor',groupedC,'PStarShowGT', false,...
    'PStarIcon',char(167),'PStarOffset',.2*abs(max(groupsE(1,:)))); box off, xlim([.5 2.5]),...
    set(gca,'xticklabel',{[]},'fontsize',14); H = sigstar(group,mpval_groupE,0,1); ...
    title('Monkey E','fontsize',14); ...
    set(H(:,2),'fontsize',14); ylim([-20e4 13e4]); hold on, ...
    scatter(ones(size(negGroupE)),negGroupE,'filled','MarkerEdgeColor',eColor,'MarkerFaceColor',mColor,'MarkerFaceAlpha',mAlpha), ...
    scatter(2*ones(size(regGroupE)),regGroupE,'filled','MarkerEdgeColor',eColor,'MarkerFaceColor',mColor,'MarkerFaceAlpha',mAlpha)

% print -dpdf -painters groupedBarGraph


%% Local Functions That I Don't Need to Be External
% Since these are super specific to this script

function cbarset(setID,yID,xPos,yPos)
if setID == 1
    c = colorbar('horiz'); c.TickDirection = 'in';
    % c.Location = 'manual'; %c.Position = [0.87 yPosition 0.025 0.2];
    c.Position = [xPos yPos+.007 0.5 0.02];
    set(c, 'fontsize', 14);
    switch yID
        case 1, set(c,'ytick',[1e5 2e5]);
        case 2, set(c,'ytick',[1e5 3e5]);
        case 3, set(c,'ytick',[2e5 10e5]);
        case 4, set(c,'ytick',[1e6 3e6]);
    end
    % colorTitleHandle = get(c, 'ylabel'); titleString = 'Power';
    % set(colorTitleHandle, 'String', titleString, 'rotation', -90, 'fontsize', 14);
    % pos = get(colorTitleHandle, 'position'); pos(1) = pos(1) + 1;
    % set(colorTitleHandle, 'position', pos);
elseif setID == 2
    c = colorbar('horiz'); c.TickDirection = 'in';
    % c.Location = 'manual'; %c.Position = [0.87 yPosition 0.025 0.2];
    c.Position = [xPos yPos+.007 0.2 0.02];
    set(c, 'fontsize', 14);
    switch yID
        case 1, set(c,'ytick',[-5e4 0 5e4]);
        case 2, set(c,'ytick',[-1e5 0 1e5]);
        case 3, set(c,'ytick',[-8e5 0]);
        case 4, set(c,'ytick',[-10e5 0]);
    end
end
end

function cbarset2(setID,xID,xPos,yPos)
if setID == 1
    c = colorbar; c.TickDirection = 'in';
    c.Location = 'manual'; %c.Position = [0.87 yPosition 0.025 0.2];
    c.Position = [xPos+.09 yPos+.015 0.0055 0.27];
    set(c, 'fontsize', 14);
    switch xID
        case 1, set(c,'ytick',[1e5 2e5]);
        case 2, set(c,'ytick',[1e5 3e5]);
        case 3, set(c,'ytick',[2e5 10e5]);
        case 4, set(c,'ytick',[1e6 3e6]);
    end
    % colorTitleHandle = get(c, 'ylabel'); titleString = 'Power';
    % set(colorTitleHandle, 'String', titleString, 'rotation', -90, 'fontsize', 14);
    % pos = get(colorTitleHandle, 'position'); pos(1) = pos(1) + 1;
    % set(colorTitleHandle, 'position', pos);
elseif setID == 2
    c = colorbar; c.TickDirection = 'in';
    % c.Location = 'manual'; %c.Position = [0.87 yPosition 0.025 0.2];
    c.Position = [xPos+.09 yPos+.0075 0.0055 0.135];
    set(c, 'fontsize', 14);
    switch xID
        case 1, set(c,'ytick',[-5e4 5e4]);
        case 2, set(c,'ytick',[-1e5 0 1e5]);
        case 3, set(c,'ytick',[-8e5 0]);
        case 4, set(c,'ytick',[-10e5 0]);
    end
end
end

function cbarset3(setID,xID,xPos,yPos,w,h)
if setID == 1
    c = colorbar; c.TickDirection = 'in';
    c.Location = 'manual'; %c.Position = [0.87 yPosition 0.025 0.2];
    c.Position = [xPos+(.76*w) yPos+(.1*h) .01 2.7*h];
    set(c, 'fontsize', 14);
    switch xID
        case 1, set(c,'ytick',[1e5 2e5]);
        case 2, set(c,'ytick',[1e5 2e5]);
        case 3, set(c,'ytick',[2e5 10e5]);
        case 4, set(c,'ytick',[1e6 3e6]);
    end
    % colorTitleHandle = get(c, 'ylabel'); titleString = 'Power';
    % set(colorTitleHandle, 'String', titleString, 'rotation', -90, 'fontsize', 14);
    % pos = get(colorTitleHandle, 'position'); pos(1) = pos(1) + 1;
    % set(colorTitleHandle, 'position', pos);
elseif setID == 2
    c = colorbar; c.TickDirection = 'in';
    % c.Location = 'manual'; %c.Position = [0.87 yPosition 0.025 0.2];
    c.Position = [xPos+(.76*w) yPos+(.15*h) .01 .9*3*h];
    set(c, 'fontsize', 14);
    switch xID
        case 1, set(c,'ytick',[-1e5 3e5]);
        case 2, set(c,'ytick',[-1e5 0 1e5]);
        case 3, set(c,'ytick',[-8e5 0]);
        case 4, set(c,'ytick',[-10e5 0]);
    end
end
end

function cbarset4(xID,xPos,yPos)
    c = colorbar; c.TickDirection = 'in';
    c.Location = 'manual'; %c.Position = [0.87 yPosition 0.025 0.2];
    c.Position = [xPos+.09 yPos+.015 0.0055 0.42];
    set(c, 'fontsize', 14);
    switch xID
        case 1, set(c,'ytick',[-5e4 0 5e4 10e4 20e4]);
        case 2, set(c,'ytick',[0e5 1e5 2e5 3e5]);
        case 3, set(c,'ytick',[-5e5 0 5e5 10e5]);
        case 4, set(c,'ytick',[0e6 1e6 2e6]);
    end
    % colorTitleHandle = get(c, 'ylabel'); titleString = 'Power';
    % set(colorTitleHandle, 'String', titleString, 'rotation', -90, 'fontsize', 14);
    % pos = get(colorTitleHandle, 'position'); pos(1) = pos(1) + 1;
    % set(colorTitleHandle, 'position', pos);

end
