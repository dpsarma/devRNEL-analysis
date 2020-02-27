%% R01  Resbub Instramuscular EMG - Part 2
% Overlay Movement Periods for selected Movements and Muscles

%% Wrist - Pro/Sup - Neutral Position - Flat

load('R:\data_raw\human\emg_decoding\EMGC14\20170424_s15\Exports\EMG\RNEL_Raw_30khz\RNEL_EMGC14_EMG_DATA_Raw_30kHz_trial_0053.mat');
data = emg_30khz';
tRef = linspace(0, length(data)/30000, length(data));

trials = readtable('R:\data_raw\human\emg_decoding\EMGC14\20170424_s15\Exports\metaTable\EMG\EMG_30kHz\emgraw_metaTrial.csv','ReadRowNames',true);
trials(9,:) = [];

load('C:\Users\dsarma\Box Sync\R-01 Resubmit Figs\Kinematics_Organized.mat');

muscleName = {'APL', 'FPL', 'SUP', 'PTER', 'ECU', 'FCR'};
taskName = trials.sTrialType(end);
mIdx = find(contains(chanelsLabel, muscleName));

Fsmooth = 10; %Lowpass cutoff frequency for smoothing the rectified signal 10 vs 50
Norder =  4;
Fs = 30e3;
dt = 1/Fs;
[b,a] = butter(Norder,Fsmooth/(0.5*Fs));
% for i = 1:size(data,1)
    emgF = filtfilt(b,a,abs(data)')';
% end

emgF_task = emgF(mIdx,:);

cycleStarts = KinematicsIdx_Structured.Sbj(14).Wrist.PS.Neut.Flat.Idx_25'-4500;
cycleStops = cycleStarts+Fs;%KinematicsIdx_Structured.Sbj(14).Wrist.PS.Neut.Flat.Idx_100+3000;
if length(cycleStops) ~= length(cycleStarts)
    cycleStarts = cycleStarts(1:length(cycleStops));
end
cycleLengths = (cycleStops - cycleStarts +1);


for c = 1:length(cycleStops)
    testEpochs{c} = emgF_task(:, cycleStarts(c):cycleStops(c));
end

%% Plot Version 1 (Overlap)
for tIdx = 1:length(cycleLengths)
    t_tmp = linspace(0,cycleLengths(tIdx)/Fs,cycleLengths(tIdx)); %
%     t_tmp = linspace(cycleStarts(tIdx)/Fs,cycleStops(tIdx)/Fs,cycleLengths(tIdx));
    epochTime{tIdx} = t_tmp; %linspace(cycleStarts(tIdx)/Fs,cycleStops(tIdx)/Fs,cycleLengths(tIdx));
end

figure; 
% maximize;%
% suptitle(['Wrist - Pro/Sup - Neutral - Flat: 10 Cycles in Series']) %
suptitle({'W - Pro/Sup - Neu - Flt:', '10 Cycles Overlapped'}); %

for idx=1:length(mIdx)
    subplot(length(mIdx),1,idx)
% %     figure; maximize;
    hold on;
    for eIdx = 1:length(testEpochs)
        plot(epochTime{eIdx},testEpochs{eIdx}(idx,:));
    end
    set(gca, 'YScale', 'log')
    ylabel([muscleName(idx) ' (uV)']);
    vline([cycleLengths/Fs]); %
   
    xlim([0 max(cycleLengths/Fs)]);%
    ylim([0 max(max(emgF_task))]);
end
xlabel('Epoch Time (sec)');

%% Plot Version 2 (Series)

for tIdx = 1:length(cycleLengths)
%     t_tmp = linspace(0,cycleLengths(tIdx)/Fs,cycleLengths(tIdx)); %
    t_tmp = linspace(cycleStarts(tIdx)/Fs,cycleStops(tIdx)/Fs,cycleLengths(tIdx));
    epochTime{tIdx} = t_tmp; %linspace(cycleStarts(tIdx)/Fs,cycleStops(tIdx)/Fs,cycleLengths(tIdx));
end

figure; 
maximize;%
suptitle(['Wrist - Pro/Sup - Neutral - Flat: 10 Cycles in Series']) %
% suptitle({'W - Pro/Sup - Neu - Flt:', '10 Cycles Overlapped'}); %

for idx=1:length(mIdx)
    subplot(length(mIdx),1,idx)
% %     figure; maximize;
    hold on;
    for eIdx = 1:length(testEpochs)
        plot(epochTime{eIdx},testEpochs{eIdx}(idx,:));
    end
    ylabel([muscleName(idx) ' (uV)']);
    vline([cycleStarts/Fs],'g'); vline([cycleStops/Fs],'r:'); %
   
%     xlim([0 max(cycleLengths/Fs)]);%
    ylim([0 max(max(emgF_task))]);
end
xlabel('Epoch Time (sec)');



%% Semilog Y
minTimeRange = min(cycleLengths);
epochTimeLog = linspace(0,minTimeRange/Fs,minTimeRange);

figure; 
% maximize;%
% suptitle(['Wrist - Pro/Sup - Neutral - Flat: 10 Cycles in Series']) %
suptitle({'(log) W - Pro/Sup - Neu - Flt:', '10 Cycles Overlapped'}); %

for muscle=1:length(mIdx)
    subplot(length(mIdx),1,muscle)
% %     figure; maximize;
    for rep = 1:length(testEpochs)
        remapEpoch(rep,:,muscle) = abs(testEpochs{rep}(muscle,:));
    end
    semilogy(epochTimeLog, remapEpoch(:,:,muscle));
    
    ylabel([muscleName(muscle) ' (log uV)']);
%     vline([cycleStarts/Fs],'g'); vline([cycleStops/Fs],'r:'); %
   
    xlim([0 max(cycleLengths/Fs)]);%
%     ylim([0 max(max(emgF_task))]);
end
xlabel('Epoch Time (sec)');

%% Semilog Y 2

% h = figure; %Toggle Comment for Subplot vs Figure.
% suptitle('W - Pro/Sup - Neu - Flt (log-scale)'); % %Toggle Comment for Subplot vs Figure.
% maximize;% %Toggle Comment for Subplot vs Figure.
% suptitle(['Wrist - Pro/Sup - Neutral - Flat: 10 Cycles in Series']) % %Toggle Comment for Subplot vs Figure.
for muscle2=1:length(mIdx)
%     subplot(length(mIdx),1,muscle2) %Toggle Comment for Subplot vs Figure.
    figure; %maximize; %Toggle Comment for Subplot vs Figure.
    for rep2 = 1:length(testEpochs)
        remapEpoch2(rep2,:,muscle2) = abs(testEpochs{rep2}(muscle2,:));
    end
    meanEpochs(muscle2,:) = mean(remapEpoch2(:,:,muscle2));
    semilogy(epochTimeLog, remapEpoch2(:,:,muscle2));
    hold on;
    semilogy(epochTimeLog, meanEpochs(muscle2,:), 'Color','k','LineWidth',1.2);
    
    ylabel([muscleName(muscle2)]);
%     vline([cycleStarts/Fs],'g'); vline([cycleStops/Fs],'r:'); %
    xlabel('Epoch Time (sec)');%Toggle Comment for Subplot vs Figure.
    
    xlim([0 max(cycleLengths/Fs)]);%
    ylim([0 10e4]);
    xticks([0 .25 .5 .75 1])
    yticks([0 10e0 10e2 10e4])
end
% xlabel('Epoch Time (sec)');  %Toggle Comment for Subplot vs Figure.
% h.Position = [1204 49 296 947]; %Toggle Comment for Subplot vs Figure.

