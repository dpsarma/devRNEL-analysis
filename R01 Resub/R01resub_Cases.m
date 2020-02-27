%% R01 Resub - Subset Cases for IM comparison with Sleeve Data
% This script plots data for comparison with sleeve data for a subset of
% previously determined muscle groups and actions from the HAPTIX EMG
% Decoding Experiments.  It relies on Data from Subject 14 (EMG14C) and the
% resultant Kinematic Indexes created by Carl Beringer.  This script is
% based off of iterative scripts by Dev Sarma. 9/6/2018

dataPath = 'R:\data_raw\human\emg_decoding\EMGC14\20170424_s15\Exports\EMG\RNEL_Raw_30khz\';
savPath = 'D:\Figures\R01-Resub\';

%% Setup Trial Table and Kinematics Index File
load('C:\Users\dsarma\Box Sync\IM EMG for R01\Kinematics_Organized.mat'); %Carl's Kinematics based Movement Indexes

trials = readtable('R:\data_raw\human\emg_decoding\EMGC14\20170424_s15\Exports\metaTable\EMG\EMG_30kHz\emgraw_metaTrial.csv','ReadRowNames',true); %Trials table form HAPTIX
trials(9,:) = []; %Remove Dead File/Trial


%% Filtering Parameters for Smoothing
Fsmooth = 10; %Lowpass cutoff frequency for smoothing the rectified signal 10 vs 50
Norder =  4;
Fs = 30e3;
dt = 1/Fs;
[b,a] = butter(Norder,Fsmooth/(0.5*Fs));

%% Establish Which Subset to Run for Plotting
actionCase = input(['Enter a number:\n' ...
    '1, Wrist - Pro/Sup - Neutral Position - Flat\n' ...
    '2, Wrist - Flx/Ext - Wrist Pronated - Neutral\n' ...
    '3, D1 - Flx/Ext - Wrist Neutral\n' ...
    '4, D1 - Flx/Ext - Wrist Pronated\n' ...
    '5, D2 - Flx/Ext - Wrist Neutral\n' ...
    '6, D2 - Flx/Ext - Wrist Supinated\n' ...
    '7, D4 - Flx/Ext - Wrist Neutral\n' ...
    '8, D4 - Flx/Ext - Wrist Pronated\n'])

switch actionCase
    case 1 %% Wrist - Pro/Sup - Neutral Position - Flat
        trialFile = 'RNEL_EMGC14_EMG_DATA_Raw_30kHz_trial_0053.mat'
        muscleName = {'APL', 'FPL', 'SUP', 'PTER', 'ECU', 'FCR'};
        cycleStarts = KinematicsIdx_Structured.Sbj(14).Wrist.PS.Neut.Flat.Idx_25'-4500;  
        
    case 2 %% Wrist - Flx/Ext - Wrist Pronated - Neutral
        trialFile = 'RNEL_EMGC14_EMG_DATA_Raw_30kHz_trial_0044.mat'
        muscleName = {'EPL', 'APL', 'FPL', 'SUP', 'PTER', 'ED3', 'FCR'};
        cycleStarts = KinematicsIdx_Structured.Sbj(14).Wrist.FE.Pro.Flat.Idx_25'-4500;
        
    case 3 %% D1 - Flx/Ext - Wrist Neutral
        trialFile = 'RNEL_EMGC14_EMG_DATA_Raw_30kHz_trial_0013.mat'
        muscleName = {'APL', 'FPL', 'SUP', 'PTER', 'ECU', 'FCR', 'EPL', 'FDP2'};
        cycleStarts = KinematicsIdx_Structured.Sbj(14).Digits(1).Neut.Idx_25'-4500;
        
    case 4 %% D1 - Flx/Ext - Wrist Pronated 
        trialFile = 'RNEL_EMGC14_EMG_DATA_Raw_30kHz_trial_0032.mat'
        muscleName = {'APL', 'FPL', 'SUP', 'PTER', 'ECU', 'FDP2'};
        cycleStarts = KinematicsIdx_Structured.Sbj(14).Digits(1).Pro.Idx_25'-4500;
        
    case 5 %% D2 - Flx/Ext - Wrist Neutral
        trialFile = 'RNEL_EMGC14_EMG_DATA_Raw_30kHz_trial_0009.mat'
        muscleName = {'FDP2', 'FCU', 'SUP', 'PTER', 'ECU', 'FCR', 'ED2', 'EPL'};
        cycleStarts = KinematicsIdx_Structured.Sbj(14).Digits(2).Neut.Idx_25'-4500;

    case 6 %% D2 - Flx/Ext - Wrist Supinated
        trialFile = 'RNEL_EMGC14_EMG_DATA_Raw_30kHz_trial_0034.mat'
        muscleName = {'FDP2', 'FDS2', 'ED3', 'ED2','SUP', 'PTER', 'FCR'};
        cycleStarts = KinematicsIdx_Structured.Sbj(14).Digits(2).Sup.Idx_25'-4500;
        
    case 7 %% D4 - Flx/Ext - Wrist Neutral
        trialFile = 'RNEL_EMGC14_EMG_DATA_Raw_30kHz_trial_0011.mat'
        muscleName = {'FDS2', 'FDS4', 'ED3', 'EPL','SUP', 'PTER', 'ED4', 'FCR'};
        cycleStarts = KinematicsIdx_Structured.Sbj(14).Digits(4).Neut.Idx_25'-4500;
        
    case 8 %% D4 - Flx/Ext - Wrist Pronated 
        trialFile = 'RNEL_EMGC14_EMG_DATA_Raw_30kHz_trial_0030.mat'
        muscleName = {'FDS2', 'FDS4', 'ED3', 'EPL','SUP', 'PTER', 'ED4', 'FCR'};
        cycleStarts = KinematicsIdx_Structured.Sbj(14).Digits(4).Pro.Idx_25'-4500;
        
    otherwise
        disp('Please enter a correct option or make a new action case!');
        return;
end


%% Task Name & Path
taskName = replace(trials({erase(trialFile, 'RNEL_')},:).sTrialType, '_', '-');
taskPath = fullfile(savPath,taskName);
mkdir(char(taskPath));

%% Load and Process Data:
load([dataPath trialFile]);
data = emg_30khz'; % Channels x Samples
tRef = linspace(0, length(data)/30000, length(data)); % Time Vector


emgF = filtfilt(b,a,abs(data)')'; % Filter Data

for count = 1:length(muscleName)
    mIdx(count) = find(contains(chanelsLabel, muscleName(count))); % Index of Desired EMG Channels
    emgF_task(count,:) = emgF(mIdx(count),:);%Extract Relevant Channels
    testNames(count,:) = chanelsLabel(mIdx(count));
end


% Identify Epochs based upon Case
cycleStops = cycleStarts+Fs; %End Period for Data Epochs
if length(cycleStops) ~= length(cycleStarts)
    cycleStarts = cycleStarts(1:length(cycleStops));
end
cycleLengths = (cycleStops - cycleStarts +1);

% Extract Epochs
for c = 1:length(cycleStops)
    testEpochs{c} = emgF_task(:, cycleStarts(c):cycleStops(c));
end

%% Plot Version 1 (Epochs in Series)

for tIdx = 1:length(cycleLengths)
%     t_tmp = linspace(0,cycleLengths(tIdx)/Fs,cycleLengths(tIdx)); %
    t_tmp = linspace(cycleStarts(tIdx)/Fs,cycleStops(tIdx)/Fs,cycleLengths(tIdx));
    epochTime{tIdx} = t_tmp; %linspace(cycleStarts(tIdx)/Fs,cycleStops(tIdx)/Fs,cycleLengths(tIdx));
end

hSeries = figure; 
maximize;%
suptitle([taskName '10 Cycles in Series']) %
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
% 
saveas(hSeries,[char(taskPath) '\' char(taskName) '_InSeries.png']);
savefig(hSeries,[char(taskPath) '\' char(taskName) '_InSeries']);


%% Semilog Y
minTimeRange = min(cycleLengths);
epochTimeLog = linspace(0,minTimeRange/Fs,minTimeRange);

hSLogY = figure; %Toggle Comment for Subplot vs Figure.
% tmptit = suptitle({char(taskName), '10 Cycles Overlapped','(log scale)'}); % %Toggle Comment for Subplot vs Figure.
% maximize;% %Toggle Comment for Subplot vs Figure.
% set(tmptit,'fontsize',10)

for muscle2=1:length(mIdx)
    hSY(muscle2) = subplot(length(mIdx),1,muscle2); %Toggle Comment for Subplot vs Figure.
%     hSY(muscle2) = figure; %maximize; %Toggle Comment for Subplot vs Figure.
    for rep2 = 1:length(testEpochs)
        remapEpoch2(rep2,:,muscle2) = abs(testEpochs{rep2}(muscle2,:));
    end
    meanEpochs2(muscle2,:) = mean(remapEpoch2(:,:,muscle2));
    semilogy(epochTimeLog, remapEpoch2(:,:,muscle2));
    hold on;
    semilogy(epochTimeLog, meanEpochs2(muscle2,:), 'Color','k','LineWidth',1.2);
    
    ylabel([muscleName(muscle2) ' (log uV)']);
%     vline([cycleStarts/Fs],'g'); vline([cycleStops/Fs],'r:'); %
%     xlabel('Epoch Time (sec)');%Toggle Comment for Subplot vs Figure.
    
    xlim([0 max(cycleLengths/Fs)]);%
    ylim([0 10e4]);
    xticks([0 .25 .5 .75 1])
    yticks([0 10e0 10e2 10e4])
%     suptitle({'(log) Wrist - Pro/Sup - Neu - Flt: ',[char(muscleName(muscle2)) ', 10 Cycles Overlapped']});
end
xlabel('Epoch Time (sec)'); %Toggle Comment for Subplot vs Figure.
linkaxes([hSY],'xy');
hSLogY.Position = [1204 49 296 947]; %Toggle Comment for Subplot vs Figure.
pause(0.1)
tmptit = suptitle({char(taskName), '10 Cycles Overlapped','(log scale)'});

saveas(hSLogY,[char(taskPath) '\' char(taskName) '_LogScale.png']);
savefig(hSLogY,[char(taskPath) '\' char(taskName) '_LogScale']);


%% Plot Version 2 (Epochs Overlapped)
for tIdx = 1:length(cycleLengths)
    t_tmp = linspace(0,cycleLengths(tIdx)/Fs,cycleLengths(tIdx)); %
%     t_tmp = linspace(cycleStarts(tIdx)/Fs,cycleStops(tIdx)/Fs,cycleLengths(tIdx));
    epochTime{tIdx} = t_tmp; %linspace(cycleStarts(tIdx)/Fs,cycleStops(tIdx)/Fs,cycleLengths(tIdx));
end

hOverlap = figure; 
% maximize;%
% suptitle(['Wrist - Pro/Sup - Neutral - Flat: 10 Cycles in Series']) %
% suptitle([taskName '10 Cycles in Series']) %

for idx=1:length(mIdx)
    hOL(idx) = subplot(length(mIdx),1,idx);
% %     figure; maximize;
    hold on;
    for rep = 1:length(testEpochs)
        remapEpoch(rep,:,idx) = abs(testEpochs{rep}(idx,:));
    end
    meanEpochs(idx,:) = mean(remapEpoch(:,:,idx));
    
    for eIdx = 1:length(testEpochs)
        plot(epochTime{eIdx},testEpochs{eIdx}(idx,:));
    end
    plot(epochTimeLog, meanEpochs(idx,:), 'Color','k','LineWidth',1.2);
    
%     set(gca, 'YScale', 'log')
    ylabel([muscleName(idx) ' (uV)']);
    vline([cycleLengths/Fs]); %
   
    xlim([0 max(cycleLengths/Fs)]);%
    ylim([0 max(max(emgF_task))]);
end
xlabel('Epoch Time (sec)');
linkaxes([hOL],'xy');
hOverlap.Position = [1204 49 296 947];
pause(0.1)
tmptit = suptitle({char(taskName), '10 Cycles Overlapped'});

saveas(hOverlap,[char(taskPath) '\' char(taskName) '_Overlap.png']);
savefig(hOverlap,[char(taskPath) '\' char(taskName) '_Overlap']);


pause(1); close all; clear all;

