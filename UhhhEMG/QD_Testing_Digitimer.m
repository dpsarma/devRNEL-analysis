%% UH3 EMG - Recruitment Curves - Quick Load Trial Stats - Testing Version %%
% The Purpose of this script is to quickly generate summary plots for the
% UH3 Aim 2: Reflexive EMG from stimulation Project.
% The Script utilizes the Stim Parameter .mat files (from Sensory Rig 1)
% and the trellis data.mat files (from Sensory Rig 2).  It will batch
% process 'n' number of files at once.  This should primarily be used for
% each step of the Recruitment Curve Process.  Please ensure that the
% target files are in the proper initial paths (RNELshare for stim params &
% the local datatank folder for trellis files).  Plots and Text Log will
% save in the designated Report Path for remote users to peruse.

%%Things to ADD: Subject Name Parameter to Automate. Remove Visibility of
%%Plots


Muscle_ID = {'Right VM', 'Right RF', 'Right VL', 'Right MG', 'Right SO med', ...
    'Right TA', 'Right LG', 'Right SO lat'};


% stimAmp = [1 2 3 4 4.5 5 5.5 6];

reportPath = 'R:\data_generated\human\uh3_stim\reports\testing\Digitimer\';

%% Identify files to load
disp('Please Select the Ripple Data Folder');
[emgFilenames, emgPathname] = uigetfile('R:\data_raw\human\uh3_stim\testing\*.nev','Pick files','MultiSelect', 'on');
disp(['User selected ', fullfile(emgPathname, emgFilenames)]);

% disp('Please Select the Folder with the Stim Trial Info:');
trialPathname = uigetdir(reportPath, 'OpenLoop Stim Trial Info');

for f = 1:length(emgFilenames)
    trialFilenames{f} = [emgFilenames{f}(1:end-4) '.mat']; %2014a friendly, no erase function
end

%% Initialize Analysis Parameters

% Filtering Settings:
fs = 30e3;
lowCut = 10; %lowest frequency to pass
highCut = 7500; %highest frequency to pass
Norder = 4;
Wp = [lowCut, highCut]/(.5*fs);
[b,a]=butter(Norder, Wp);

%Smoothing Settings:
Fsmooth = 10; %Lowpass cutoff frequency for smoothing the rectified signal 10 vs 50
NorderS =  4;
dt = 1/fs;
[bS,aS] = butter(NorderS,Fsmooth/(0.5*fs));

% Set Epoch Window for Stim-triggered Averaging
dtPre = 20e-3; % msec before stim onset 
dtPost = 80e-3; % msec after stim onset 


%% Begin Parsing each file:
for fnum = 7%1:length(emgFilenames)

%     trialName(fnum,:) = erase(trialFilenames{fnum}, '.mat'); %only for
%     2016 and beyond
    trial(fnum).Name = trialFilenames{fnum}(1:end-4);
    % Load Trial Stim Parameters
    if exist(fullfile(trialPathname,trialFilenames{fnum}), 'file')==0
        disp(['This Trial.mat DOES NOT EXIST - ' trial(fnum).Name])
        fileExist{fnum}.trial = 0;
        fileExist{fnum}.name = trial(fnum).Name;
        continue;
    else
        disp(['Trial: ' trialFilenames{fnum}(1:end-4)]);
        fileExist{fnum}.trial = 1;
        fileExist{fnum}.name = trial(fnum).Name;
    end
       
    trialInfo = load(fullfile(trialPathname,trialFilenames{fnum}));
    
    if isempty(trialInfo.Stim_params)
        disp(['This Trial.mat IS EMPTY - ' trial(fnum).Name]);
        fileExist(fnum).trial = 2; %Where 2 equals exists but empty
        continue;
    end
    
    %% Identify Stim Events Channel
    spinalElec(fnum) = cell2mat(trialInfo.Stim_params.SpinalElecs);
    nevStimCh(fnum) = trialInfo.Stim_params.NSChannels{1,1}(1);%find_stimevtchannel(spinalElec);
    
    % Set Pulse Width, Stim Duration, & Stim Frequency
    stimAmp(fnum) = cell2mat(trialInfo.Stim_params.SpinalElecAmps);
    pulseWidth(fnum) = cell2mat(trialInfo.Stim_params.PulseWidth); %this is provided in ms
    stimDuration(fnum) = cell2mat(trialInfo.Stim_params.Duration)/1000; %this is provided in ms
    stimFrequency(fnum) = cell2mat(trialInfo.Stim_params.Frequency); %this is provided in Hz
    numStims(fnum) = stimDuration(fnum)*stimFrequency(fnum);
    nSampsPre = floor(dtPre*fs);
    nSampsPost = floor(dtPost*fs);

    if contains(trial(fnum).Name,'Set008')
        stimAmp(fnum) = 5900; %Correct for weird '76' fault
    end
    
    if contains(trial(fnum).Name,'Set012')
        stimAmp(fnum) = 5000 + cell2mat(trialInfo.Stim_params.SpinalElecAmps); %To go beyond 12mA
        if stimAmp(fnum) == 10200
            stimAmp(fnum) = 12000;
        elseif stimAmp(fnum) == 10400
            stimAmp(fnum) = 14000;
        elseif stimAmp(fnum) == 10600
            stimAmp(fnum) = 16000;
        elseif stimAmp(fnum) == 10800
            stimAmp(fnum) = 18000;
        elseif stimAmp(fnum) == 10800
            stimAmp(fnum) = 18000;
        elseif stimAmp(fnum) == 11000
            stimAmp(fnum) = 20000; 
        end
    end
    
          
    
    stimLength(fnum) = (pulseWidth(fnum)/1000)*fs; %Pulse Width in ms, interstim interval is 2 ticks of clock (2/fs)
    trial(fnum).epochTimeVec = linspace(-dtPre, (stimLength(fnum)/fs) + dtPost, nSampsPre + nSampsPost + stimLength(fnum));
    
    [stimEvts] = read_digitalEvents([emgPathname emgFilenames{fnum}],1);
    trial(fnum).stims = floor(cell2mat(stimEvts.timeStamp)*fs);
    
  
    [analogData,timeVec] = read_continuousData([emgPathname emgFilenames{fnum}], 'raw' , 256+[1:32]); %128 vs 256
    
    %Monopolar Recordings to Bipolar EMG
    for i = 1:2:size(analogData,1)
        epochRMS = analogData(i:i+1,:);
        emg(round(i/2),:) = diff(flipud(epochRMS));
    end
    
    emgF = filtfilt(b,a,emg');
    emgF = emgF';
    
    emgS = filtfilt(bS,aS,emgF')';
    % Epoch Windows 
    for i = 1:size(emg,1)
        asd = emgF(i,:);
        emgStim{i} = cell2mat(arrayfun(@(x) asd((trial(fnum).stims(x)-(nSampsPre)):(trial(fnum).stims(x)+(stimLength(fnum)+nSampsPost-1))), 1:length(trial(fnum).stims), 'UniformOutput',false)');
        qwe = emgS(i,:);
        emgStim_Sm{i} = cell2mat(arrayfun(@(x) asd((trial(fnum).stims(x)-(nSampsPre)):(trial(fnum).stims(x)+(stimLength(fnum)+nSampsPost-1))), 1:length(trial(fnum).stims), 'UniformOutput',false)');
         
    end
    
    %Plot Stim averages
    hF = figure;
    for c = 1:length(Muscle_ID)
        hFs(c) = subplot(2,round(length(Muscle_ID)/2),c);
        
     %Plot /Channel
     trial(fnum).meanTrace{c} = mean(emgStim{1,c});
     plot(trial(fnum).epochTimeVec,emgStim{1,c});hold on;
     plot(trial(fnum).epochTimeVec,trial(fnum).meanTrace{c},'k','LineWidth',1.2);
     vline([0 stimLength(fnum)/fs],'r-');
     title([Muscle_ID{c}]);
              
    end
    suplabel('Amplitude(uV)','y');
    suplabel('Time(sec)','x');
    
    Pix_SS = get(0,'screensize');
    hF.Position = Pix_SS;
    tmptit = suptitle([trial(fnum).Name ' Stim Traces']);
    set(tmptit, 'Interpreter', 'none'); 
    linkaxes([hFs],'xy');
    hFs(1).XLim = [trial(fnum).epochTimeVec(1) trial(fnum).epochTimeVec(end)];

    disp('Making Trial Folder')
    mkdir(reportPath,trial(fnum).Name)
    savefig(hF,[fullfile(reportPath,trial(fnum).Name) '\EpochedStimTraces_R']);
    saveas(hF,[fullfile(reportPath,trial(fnum).Name) '\EpochedStimTraces_R.png']);

    %% Get RMS acros epochs
    trial(fnum).epochRMS = cellfun(@rms,emgStim,'UniformOutput', false);
%     figure;cellfun(@plot,epochRMS);
    StimEnd = nSampsPre + stimLength(fnum);
    
    hF2 = figure;
    for c = 1:length(Muscle_ID)
              
     %%Plot /Channel
     plot(trial(fnum).epochTimeVec,trial(fnum).epochRMS{c});hold on;
          vline([0 stimLength/fs],'r-'); %vline(stimEnd/fs,'g-');
%      trial(fnum).responseRMS(c,:) = mean(trial(fnum).epochRMS{c}(StimEnd:end));
                  
    end
    xlabel('RMS Amplitude (uV)');
    ylabel('Time (Sec)');
    title([trial(fnum).Name ' - Epoched RMS'],'Interpreter', 'none');
    vline(stimEnd/fs);
    legend(Muscle_ID);
    
    saveas(hF2,[fullfile(reportPath,trial(fnum).Name) '\Epoched_RMS_R.png']);
   
    
    
    
    close all
    
    
    clearvars -except emgFilenames emgPathname reportPath Muscle_ID fs b a bS aS ...
                dtPre dtPost fnum spinalElec nevStimCh stimAmp pulseWidth...
                stimDuration stimFrequency responseRMS trialName ...
                trialPathname trialFilenames trial
end

%% Plot Recruitment Curves

responseRMS = [trial(:).responseRMS]';

hF3 = figure;
C = linspecer(length(Muscle_ID),'qualitative');

% % % subplot(2,3,1)
% % for m = 1:length(Muscle_ID)
% %     yyS = smooth(stimAmp/1000,responseRMS(:,m),'sgolay');
% %     plot(stimAmp/1000,yyS,'Color',C(m,:),'LineWidth',3);
% %     hold on;
% % end
% %     ylabel('(post-Stim) Mean RMS  (uV)')
% %     xlabel('Stimulation Amplitude (mA)')
% %     title('Recruitment Curve - Amplitude - 1to6mA, Standing');
    
% subplot(2,3,2)
% % for m = 1:length(Muscle_ID)/2
% % %     scatter(pulseWidth,responseRMS(:,m)','Color',C(m,:))
% %     yyP = smooth(pulseWidth,responseRMS(:,m),'sgolay');
% %     plot(pulseWidth,yyP,'Color',C(m,:),'LineWidth',3);
% %     hold on;
% % end
% % legend([Muscle_ID(1:length(Muscle_ID))],'Location','northoutside', 'Orientation','horizontal');
% %     ylabel('(post-Stim) Mean RMS  (uV)')
% %     xlabel('Stimulation PulseWidth (us)')
% %     title('Recruitment Curve - PulseWidth');
%     
% subplot(2,3,3)
for m = 1:length(Muscle_ID)/2
%     scatter(stimFrequency,responseRMS(:,m)')
    yyF = smooth(stimFrequency(1:6),responseRMS(:,m),'sgolay');
    plot(stimFrequency(1:6),yyF,'Color',C(m,:),'LineWidth',3);
    hold on;
end
    ylabel('(post-Stim) Mean RMS  (uV)')
    xlabel('Stimulation Frequency (Hz)')
    title('Recruitment Curve - Frequency');
%     
%     subplot(2,3,4)
% for m = (length(Muscle_ID)/2):length(Muscle_ID)
%     plot(stimAmp/1000,responseRMS(:,m)')
%     hold on;
% end
%     ylabel('(post-Stim) Mean RMS  (uV)')
%     xlabel('Stimulation Amplitude (mA)')
%     title('Recruitment Curve - Amplitude');
%     
% subplot(2,3,5)
% for m = (length(Muscle_ID)/2):length(Muscle_ID)
%     scatter(pulseWidth,responseRMS(:,m)')
%     hold on;
% end
% legend([Muscle_ID((length(Muscle_ID)/2)+1:length(Muscle_ID))],'Location','southoutside', 'Orientation','horizontal');
%     ylabel('(post-Stim) Mean RMS  (uV)')
%     xlabel('Stimulation PulseWidth (ms)')
%     title('Recruitment Curve - PulseWidth');
%     
% subplot(2,3,6)
% for m = (length(Muscle_ID)/2):length(Muscle_ID)
%     scatter(stimFrequency,responseRMS(:,m)')
%     hold on;
% end
%     ylabel('(post-Stim) Mean RMS  (uV)')
%     xlabel('Stimulation Frequency (Hz)')
%     title('Recruitment Curve - Frequency');
    
Pix_SS = get(0,'screensize');
hF3.Position = Pix_SS;
    
saveas(hF3,[reportPath '\' strrep(datestr(now),':','_') '_FreqRecruitmentCurves.png']);
savefig(hF3,[reportPath '\' strrep(datestr(now),':','_') '_FreqRecruitmentCurves']);
    
    

varNames = {'TrialName','SpinalElectrodes', 'NEV_StimChannel', 'StimAmplitudes_mA', ...
    'PulseWidths_ms', 'StimDurations_s', 'StimFrequencies_Hz'};

TrialNames = [{trial(:).Name}'];

TrialTable = table(TrialNames, spinalElec', nevStimCh', (stimAmp/1000)', pulseWidth',...
    stimDuration', stimFrequency','VariableNames', varNames)

 
writetable(TrialTable,[reportPath '\' strrep(datestr(now),':','_') '_Summary.txt'],'Delimiter','tab')
writetable(TrialTable,[reportPath '\' strrep(datestr(now),':','_') '_Summary.xls'])



for m = 1:length(Muscle_ID)
    for i = 1:length(stimFrequency)
        Z(i,:) = [trial(i).meanTrace{1,m}];
%         Y(i,:) = [trial(i).epochTimeVec];
%         X(i,:) = [stimAmp(i)/1000];
    end
    hWs(m) = figure;
    maximize;
    waterfall(Z)
    title({'Mean EMG Response - Digitimer - 4.5mA, Frequency Modulated, Seated', Muscle_ID{m}});
    zlim([min(min(Z)) inf])
    xLAB = xticks/fs;
    xlabel('Samples (30K)');
    ylabel(['Trial #', 'Modulating Frequency 1-50Hz']);
    zlabel('Voltage(uV)')
    view([-2.7, 29.2]);
    
    saveas(hWs(m),[reportPath 'MeanTraces-Set5_' num2str(m) '.png']);
    savefig(hWs(m),[reportPath 'MeanTraces-Set5_' num2str(m)]);
end
% linkaxes([hWs],'xyz');
% legend('1 mA', '2mA', '3mA', '4mA', '4.5mA', '5mA', '5.5mA', '6mA','Location','southoutside');