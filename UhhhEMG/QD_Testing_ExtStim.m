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

Muscle_ID = {'Right SO (lat)', 'Right RF', 'Right VL', 'Right BF', 'Right SE', ...
    'Right TA', 'Right SO (med)', 'Right LG'};


stimAmp = [1 2 3 4 4.5 5 5.5 6];

reportPath = 'R:\data_generated\human\uh3_stim\reports\testing\StimCompare\';

%% Identify files to load
disp('Please Select the Ripple Data Folder');
[emgFilenames, emgPathname] = uigetfile('C:\DataTanks\2018\*.nev','Pick files','MultiSelect', 'on');
disp(['User selected ', fullfile(emgPathname, emgFilenames)]);

% disp('Please Select the Folder with the Stim Trial Info:');
% trialPathname = uigetdir(reportPath, 'OpenLoop Stim Trial Info');

for f = 1:length(emgFilenames)
    trialFilenames{f} = [emgFilenames{f}(1:end-4) '.mat']; %2014a friendly, no erase function
end

%% Initialize Analysis Parameters

% Filtering Settings:
fs = 30e3;
lowCut = 75; %lowest frequency to pass
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
artifactBuffer = 0.015; % post stim artifact settle lag in s, so 0.01 is 10ms
rmsWindow = 0.05; % post stim RMS Window in s, so 0.05 is 50ms

%% Begin Parsing each file:
for fnum = 1:length(emgFilenames)

%     trialName(fnum,:) = erase(trialFilenames{fnum}, '.mat'); %only for
%     2016 and beyond
    trial(fnum).Name = trialFilenames{fnum}(1:end-4);
    %% Load Trial Stim Parameters
%     if exist(fullfile(trialPathname,trialFilenames{fnum}), 'file')==0
%         disp(['This Trial.mat DOES NOT EXIST - ' trialName(fnum,:)])
%         fileExist{fnum}.trial = 0;
%         fileExist{fnum}.name = trialName(fnum,:);
%         continue;
%     else
%         disp(['Trial: ' trialFilenames{fnum}(1:end-4)]);
%         fileExist{fnum}.trial = 1;
%         fileExist{fnum}.name = trialName(fnum,:);
%     end
       
%     trialInfo = load(fullfile(trialPathname,trialFilenames{fnum}));
%     
%     if isempty(trialInfo.Stim_params)
%         disp(['This Trial.mat IS EMPTY - ' trialName(fnum,:)]);
%         fileExist(fnum).trial = 2; %Where 2 equals exists but empty
%         continue;
%     end
    
    %% Identify Stim Events Channel
%     spinalElec(fnum) = cell2mat(trialInfo.Stim_params.SpinalElecs);
%     nevStimCh(fnum) = trialInfo.Stim_params.NSChannels{1,1}(1);%find_stimevtchannel(spinalElec);
    
    % Set Pulse Width, Stim Duration, & Stim Frequency
%     stimAmp(fnum) = cell2mat(trialInfo.Stim_params.SpinalElecAmps);
    trial(fnum).pulseWidth = 1; %cell2mat(trialInfo.Stim_params.PulseWidth); %this is provided in ms
    trial(fnum).stimDuration = 6000; %cell2mat(trialInfo.Stim_params.Duration)/1000; %this is provided in ms
    trial(fnum).stimFrequency = 1; %cell2mat(trialInfo.Stim_params.Frequency); %this is provided in Hz
%     numStims(fnum) = stimDuration(fum)*stimFrequency(fnum);
    trial(fnum).nSampsPre = floor(dtPre*fs);
    trial(fnum).nSampsPost = floor(dtPost*fs);

    trial(fnum).stimLength = (trial(fnum).pulseWidth/1000)*fs; %Pulse Width in ms, interstim interval is 2 ticks of clock (2/fs)
    trial(fnum).epochTimeVec = linspace(-dtPre, (trial(fnum).stimLength/fs)+dtPost, trial(fnum).nSampsPre+trial(fnum).nSampsPost+trial(fnum).stimLength);
    
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
        emgStim{i} = cell2mat(arrayfun(@(x) asd((trial(fnum).stims(x)-(trial(fnum).nSampsPre)):(trial(fnum).stims(x)+(trial(fnum).stimLength+trial(fnum).nSampsPost-1))), 1:length(trial(fnum).stims), 'UniformOutput',false)');
        qwe = emgS(i,:);
        emgStim_Sm{i} = cell2mat(arrayfun(@(x) asd((trial(fnum).stims(x)-(trial(fnum).nSampsPre)):(trial(fnum).stims(x)+(trial(fnum).stimLength+trial(fnum).nSampsPost-1))), 1:length(trial(fnum).stims), 'UniformOutput',false)');
         
    end
    
%     %Plot Stim averages
%     hF = figure;
%     for c = 1:length(Muscle_ID)
%         hFs(c) = subplot(2,round(length(Muscle_ID)/2),c);
%         
%      %Plot /Channel
%      trial(fnum).meanTrace{c} = mean(emgStim{1,c});
%      plot(trial(fnum).epochTimeVec,emgStim{1,c});hold on;
%      plot(trial(fnum).epochTimeVec,trial(fnum).meanTrace{c},'k','LineWidth',1.2);
%      vline([0 trial(fnum).stimLength/fs],'r-');
%      title([Muscle_ID{c}]);
%               
%     end
%     suplabel('Amplitude(uV)','y');
%     suplabel('Time(sec)','x');
%     
%     Pix_SS = get(0,'screensize');
%     hF.Position = Pix_SS;
%     tmptit = suptitle([trial(fnum).Name ' Stim Traces']);
%     set(tmptit, 'Interpreter', 'none'); 
%     linkaxes([hFs],'xy');
%     hFs(1).XLim = [trial(fnum).epochTimeVec(1) trial(fnum).epochTimeVec(end)];
% 
%     disp('Making Trial Folder')
%     mkdir(reportPath,trial(fnum).Name)
%     savefig(hF,[fullfile(reportPath,trial(fnum).Name) '\EpochedStimTraces_R']);
%     saveas(hF,[fullfile(reportPath,trial(fnum).Name) '\EpochedStimTraces_R.png']);
% 
%     %% Get RMS acros epochs
%     trial(fnum).epochRMS = cellfun(@rms,emgStim,'UniformOutput', false);
% %     figure;cellfun(@plot,epochRMS);
    stimEnd = trial(fnum).nSampsPre + trial(fnum).stimLength;
    
%     hF2 = figure;
%     for c = 1:length(Muscle_ID)
%               
%      %Plot /Channel
%      plot(trial(fnum).epochTimeVec,trial(fnum).epochRMS{c});hold on;
%      %      vline([0 stimLength/fs],'r-');
%      trial(fnum).responseRMS(c,:) = mean(trial(fnum).epochRMS{c}(stimEnd:end));
%                   
%     end
%     xlabel('RMS Amplitude (uV)');
%     ylabel('Time (Sec)');
%     title([trial(fnum).Name ' - Epoched RMS'],'Interpreter', 'none');
%     vline(stimEnd/fs);
%     legend(Muscle_ID);
    
%     saveas(hF2,[fullfile(reportPath,trial(fnum).Name) '\Epoched_RMS_R.png']);
   
     trial(fnum).emgStim = emgStim;
    trial(fnum).emgStim_Sm = emgStim_Sm;
    
    %% Get RMS acros epochs
    trial(fnum).epochRMS = cellfun(@rms,emgStim,'UniformOutput', false);
        
    %% Extract "Mean" Trace & Response
%     stimEnd = nSampsPre + stimLength + artifactBuffer*fs;    
    for m = 1:length(Muscle_ID)
        % Mean Trace
        trial(fnum).meanTrace{m} = mean(emgStim{1,m});
        trial(fnum).meanTrace_Sm{m} = mean(emgStim_Sm{1,m});
        % Response
        trial(fnum).responseRMS(m,:) = mean(trial(fnum).epochRMS{m}(stimEnd:stimEnd+(rmsWindow*fs)));
        trial(fnum).responseP2P(m,:) = mean(trial(fnum).emgStim{m}(stimEnd:stimEnd+(rmsWindow*fs)));
    end

    
    
    close all
    
    
    clearvars -except emgFilenames emgPathname reportPath Muscle_ID fs b a bS aS ...
                dtPre dtPost fnum spinalElec nevStimCh stimAmp pulseWidth...
                stimDuration stimFrequency responseRMS trialName ...
                trialPathname trialFilenames trial rmsWindow
end

%% Plotting Stim Averages
for c = 1:length(Muscle_ID)
    disp('Plotting Stim Averages');
    
    hF(c) = figure; %('visible', 'off');
    tmpsup = suptitle([Muscle_ID{c} ' - External Stim -  Stim Traces']);
    set(tmpsup, 'Interpreter', 'none') 
    maximize;
%     set(hF(c),'Position', [1 -798 1080 1803]);
    subCount = 1;
    for fnum = 1:length(emgFilenames)
        hFs(fnum) = subplot(length(emgFilenames),1,fnum);%round(fnum/subcount)); evens by 2
%         hFs(fnum) = subplot(12,1,subCount);
        %Plot /Channel
        plot(trial(fnum).epochTimeVec,trial(fnum).emgStim{c});hold on;
        plot(trial(fnum).epochTimeVec,trial(fnum).meanTrace{c},'k','LineWidth',1.2);
        vline([0    trial(fnum).stimLength/fs],'r-');
        tmptit = title([trial(fnum).Name ' - ' num2str(stimAmp(fnum)) 'mA, ' ...
            num2str(trial(fnum).pulseWidth) 'us, ' num2str(trial(fnum).stimFrequency) 'Hz, ' num2str(trial(fnum).stimDuration) 's']); 
        set(hFs(fnum),'FontSize',7);
        set(tmptit, 'Interpreter', 'none','FontSize',5);
        hFs(fnum).XLim = [trial(fnum).epochTimeVec(1) trial(fnum).epochTimeVec(end)];
        
        subCount = subCount + 1;
    end
    linkaxes([hFs],'xy');
    suplabel('Amplitude(uV)','y');
    suplabel('Time(sec)','x');
%     hF(c).Position = [317 49 568 1787];  
%     
    disp('Saving to Set Folder')
    savefig(hF(c),[reportPath Muscle_ID{c} '_EpochedStimTraces_ExtStim' ]);
    saveas(hF(c),[reportPath Muscle_ID{c} '_EpochedStimTraces_ExtStim.png']);
end


% % %% Plot Recruitment Curves
% % 
% % responseRMS = [trial(:).responseRMS]';
% % 
% % hF3 = figure;
% % 
% % % subplot(2,3,1)
% % for m = 1:length(Muscle_ID)
% %     plot(stimAmp/1000,responseRMS(:,m))
% %     hold on;
% % end
% %     ylabel('(post-Stim) Mean RMS  (uV)')
% %     xlabel('Stimulation Amplitude (mA)')
% %     title('Recruitment Curve - Amplitude - Ext. Stim');
% % 
% % legend([Muscle_ID(1:length(Muscle_ID))],'Location','northoutside', 'Orientation','horizontal');
% % 
% %     
% % Pix_SS = get(0,'screensize');
% % hF3.Position = Pix_SS;
% %     
% % saveas(hF3,[reportPath '\' strrep(datestr(now),':','_') '_AmpRecruitmentCurvesExt.png']);
% % savefig(hF3,[reportPath '\' strrep(datestr(now),':','_') '_AmpRecruitmentCurvesExt']);
% %     
%     
% 
% varNames = {'TrialName','SpinalElectrodes', 'NEV_StimChannel', 'StimAmplitudes_mA', ...
%     'PulseWidths_ms', 'StimDurations_s', 'StimFrequencies_Hz'};
% 
% 
% TrialTable = table(string(trialName), (stimAmp/1000)', pulseWidth',...
%     stimDuration', stimFrequency','VariableNames', varNames)
% 
%  
% writetable(TrialTable,[reportPath '\' strrep(datestr(now),':','_') '_Summary.txt'],'Delimiter','tab')
% writetable(TrialTable,[reportPath '\' strrep(datestr(now),':','_') '_Summary.xls'])



for m = 1:length(Muscle_ID)
    for i = 1:8
        XY(i,:) = [trial(i).meanTrace{1,m}];
    end
    hWs(m) = figure;
    maximize;
    waterfall(XY)
    title({'Mean EMG Response - Digitimer', Muscle_ID{m}});
    zlim([-1500 inf])
    xLAB = xticks/fs;
    xlabel({'Samples (30K)', join(string(xLAB(1:end-1)))});
    ylabel(['Trial #, Stim (mA):' join(string(stimAmp),',')]);
    zlabel('Voltage(uV)')
    view([-2.7, 29.2]);
    
    saveas(hWs(m),[reportPath 'MeanTraces-ExtStim_' num2str(m) '.png']);
    savefig(hWs(m),[reportPath 'MeanTraces-ExtStim_' num2str(m)]);
end
% linkaxes([hWs],'xyz');
% legend('1 mA', '2mA', '3mA', '4mA', '4.5mA', '5mA', '5.5mA', '6mA','Location','southoutside');