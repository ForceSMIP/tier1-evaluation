
% Used to make Fig. 12 of Wills et al. 2025
% to be run after ForceSMIP_index_skill.m

index = 1:30; 
examples = [6 7 21 24 25];

switch length(index_ref)
    case 73
        time = 1950:2022;
    case 64
        time = 1954.5:2017.5;
    case 63 
        time = 1955:2017;
        if max(index_ref)>50
            index_ref = index_ref./100;
            index_all = index_all./100;
        end
end
index_estimate = index(RMSs(3:end)<RMSs(2)&STDs(3:end)>0.5&CORs(3:end)>0.1);
q50 = quantile(index_all(:,index_estimate),0.50,2);
q17 = quantile(index_all(:,index_estimate),0.17,2);
q83 = quantile(index_all(:,index_estimate),0.83,2);
figure; shadedErrorBar(time,q50,[q83-q50, q50-q17]');
hold on; plot(time,index_ref,'k','linewidth',2)
hold on; plot(time,index_all(:,26),'k--','linewidth',2)
hold on; plot(time,index_all(:,examples(1)),'linewidth',2)
hold on; plot(time,index_all(:,examples(2)),'linewidth',2)
hold on; plot(time,index_all(:,examples(3)),'linewidth',2)
hold on; plot(time,index_all(:,examples(4)),'linewidth',2)
hold on; plot(time,index_all(:,examples(5)),'linewidth',2)
set(gca,'ygrid','on')

% switch variable
%     case 'tos'
%         pretty_figure(700,300,'Year','Niño3.4 Anomaly (°C)','none','none',16); set(gca,'xlim',[1954.5 2017.5])
%         %pretty_figure(700,300,'Year','NASSTI AMV (°C)','none','none',16); set(gca,'xlim',[1954.5 2017.5])
%     case 'tas'
%         pretty_figure(700,300,'Year','GMST Anomaly (°C)','none','none',16); set(gca,'xlim',[1950 2022])
%     case 'pr'
%         pretty_figure(700,300,'Year','Sahel Precip. Anomaly (mm day^{-1})','none','none',16); set(gca,'xlim',[1954.5 2017.5])
%     case 'psl'
%         pretty_figure(700,300,'Year','Aleutian Low Anomaly (hPa)','none','none',16); set(gca,'xlim',[1955 2017])
%     case 'monmaxtasmax'
%         pretty_figure(700,300,'Year',{'Continental Europe','TXx Anomaly (°C)'},'none','none',16); set(gca,'xlim',[1950 2022])
% end