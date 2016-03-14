%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Wei's practice on Chapter 9%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% compute ERP at each electrode
load sampleEEGdata.mat

for i=1:EEG.nbchan
erp(i,:)=[double(squeeze(mean(EEG.data(i,:,:),3)))]; 
end

%% select five time points at which to show topographical plots

times2plot = dsearchn(EEG.times',(0:100:400)');
color_limit = 10; 

figure
set(gcf,'Name',[ num2str(length(times2plot)) ' topographical plots '])

for i=1:length(times2plot)
    
    % figure out how many subplots we need
    subplot(ceil(length(times2plot)/ceil(sqrt(length(times2plot)))),ceil(sqrt(length(times2plot))),i)
    
    % extract EEG data 
    topoplot(double(squeeze(mean(EEG.data(:,times2plot(i),:),3))),EEG.chanlocs,'maplimits',[-color_limit color_limit]);
    title([ num2str(round(EEG.times(times2plot(i)))) ' ms' ])
end

%% Show average of 40 ms at each each time point
% option 1
average_points = ceil(20/(1/EEG.srate*1000))

figure
set(gcf,'Name',[ num2str(length(times2plot)) ' topographical plots at average of 50 ms'],'Number','off')

for i=1:length(times2plot)
    
    % figure out how many subplots we need
    subplot(ceil(length(times2plot)/ceil(sqrt(length(times2plot)))),ceil(sqrt(length(times2plot))),i)
    
    % extract EEG data 
    eegdata2plot = mean(double(squeeze(mean(EEG.data(:,(times2plot(i)-average_points:times2plot(i)+average_points),:),2))),2);
    
    topoplot(eegdata2plot,EEG.chanlocs,'maplimits',[-color_limit color_limit]);
    title([ num2str(round(EEG.times(times2plot(i)))) ' ms' ])
end

% option 2

times=[0:100:400];

figure
set(gcf,'Name',[ num2str(length(times2plot)) ' topographical plots at average of 50 ms'],'Number','off')

for i=1:length(times2plot)
    
    % figure out how many subplots we need
    subplot(ceil(length(times2plot)/ceil(sqrt(length(times2plot)))),ceil(sqrt(length(times2plot))),i)
    
    % extract EEG data 
    eegdata2plot = mean (double(squeeze(mean(EEG.data(:,(EEG.times> times(i)-20 & EEG.times < times(i)+20),:),3)))')';
    
    topoplot(eegdata2plot,EEG.chanlocs,'maplimits',[-color_limit color_limit]);
    title([ num2str(round(EEG.times(times2plot(i)))) ' ms' ])
end



%% find peak times for each electrodes

for i=1:size(erp,1)
    
    [m,n]=max(erp(i,:));
    time_peaks(i)=EEG.times(n);
    
end

time_peaks = time_peaks(time_peaks>100 & time_peaks <400)


%% topo for peak times

time_peaks = unique(time_peaks)


figure (1)
set(gcf,'Name',[  ' Topographical plots at ' num2str(length(time_peaks)) ' peak times'],'Number','off')

for i=1:length(time_peaks)
    
    % figure out how many subplots we need
    subplot(ceil(length(time_peaks)/ceil(sqrt(length(time_peaks)))),ceil(sqrt(length(time_peaks))),i)
    
    % extract EEG data 
    topoplot(double(squeeze(mean(EEG.data(:,(EEG.times==time_peaks(i)),:),3))),EEG.chanlocs,'maplimits',[-color_limit color_limit]);
    title([ num2str(round(time_peaks(i))) ' ms' ])
    % colorbar 
end
%set(gcf,'position',[100 100 1100 650]) % sets figure size
% fig1Pos = get(figure(1),'position');
hb = colorbar('location','eastoutside');
% set(figure(1),'Units','normalized', 'position', [0 0 0.9 0.9]);
set(hb,'Units','normalized', 'position', [0.81 0.1 0.03 0.22]);





