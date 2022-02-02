% Find time points that have very large power for each channel and see if
% they are consistent across all channels to remove artifacts

% define signal here
signal = signalE_allpost_bt;

chans = 1:size(signal,2);

% First calculate the power at 1 second intervals
dt = 0.5*Fs; % number of samples every 1 second
Ndt = size(signal,1)/dt; % number of 1 second intervals in 30 minutes

%% Find the number of very high power time points for each channel

% First normalize the power values for each channel;
sigN = normalize(signal);

% determine threshold defining artifacts (#standard deviations)
thresh = 25;

figure, set(gcf,'color','w','units','normalize','outerposition',[0 0 1 1])
plot(sigN), hold on, plot(thresh*ones(size(sigN,1),1)), drawnow

highVals = cell(1,length(chans));

for iChan = chans
    highVals{1,iChan} = find(abs(sigN(:,iChan)) > thresh);    
end

%% In How Many Channels is Each Artifact Repeated?

% This is to make sure artifacts are happening in all channels

tic
figure, set(gcf,'color','w','units','normalize','outerposition',[0 0 1 1])
for iChan = chans
    plot(iChan*ones(length(highVals{1,iChan})),highVals{1,iChan},'o')
    xlim([0 32]), %ylim([0 2e6]) 
    hold on, ylabel('Beta'), box off
    toc
end
xlabel('Channels')
[~,hy] = suplabel('Sample Index','y');
