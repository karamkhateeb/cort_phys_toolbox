%% Determine Which Channels are Bad Channels

%% Monkey D

figure('color','w','units','normalize','outerposition',[0 0 1 1])
channels = 1:32;

for iChan = channels
    
    subplot(1,2,1)
    plot(pwrD_allbase_th,'k','handlevisibility','off'), hold on
    plot(pwrD_allbase_th(:,iChan),'r','linewidth',4)
    box off, %ylim([0 1e6])
    title('Monkey E Baseline')
    ylabel('Low Gamma Power')
    hold off
    
    subplot(1,2,2)
    plot(pwrD_allpost_th,'k'), hold on
    plot(pwrD_allpost_th(:,iChan),'r','linewidth',4)
    box off, %ylim([0 2e6])
    title('Monkey E Post PT')
    hold off
    [~,ht] = suplabel(['Channel ' num2str(iChan)],'t');
    set(ht,'fontsize',18,'fontweight', 'bold'); 
    drawnow, pause(1)
end

%% Monkey E

figure('color','w','units','normalize','outerposition',[0 0 1 1])
channels = 1:32;

for iChan = channels
    
    subplot(1,2,1)
    plot(pwrD_base_hg,'k','handlevisibility','off'), hold on
    plot(pwrD_base_hg(:,iChan),'r','linewidth',4)
    box off
    title('Monkey D Baseline')
    ylabel('High Gamma Power')
    hold off
    
    subplot(1,2,2)
    plot(pwrD_post_hg,'k'), hold on
    plot(pwrD_post_hg(:,iChan),'r','linewidth',4)
    box off
    title('Monkey D Post PT')
    hold off
    [~,ht] = suplabel(['Channel ' num2str(iChan)],'t');
    set(ht,'fontsize',18,'fontweight', 'bold');
    drawnow, pause(1)
end
