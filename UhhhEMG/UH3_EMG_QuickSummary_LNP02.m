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

%% FOR THE EXPERIMENTER:
% 1) Designate the Subject Name. 
% 2) Change Report Paths & Data Tank Path if needed. ('Stimparam'.mat files)
% 2) Ascribe Set (session/block, etc.).
% 3) Add a description (i.e. 1mA-6mA_Standing).
% 4) Set RMS Window and artifact settle buffer as needed.


addpath(genpath('\\share.files.pitt.edu\RnelShare\data_generated\human\uh3_stim\genSummary\code_tools')); 
path_datatank = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\'; %'C:\DataTanks\2018\'; %'R:\data_raw\human\uh3_stim\';
path_stimparams = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\'; %'R:\data_generated\human\uh3_stim\genSummary\'; %'R:\data_raw\human\uh3_stim\';

%% Experiment Session Parameters %'Limit Check f-2Hz PW-1ms' 'RC 500uA-5mA f-1Hz PW-200us'
% subjectName = 'LSP02b'; setName = 'Elec 1'; setDescrpt = 'Day 15 - 100Hz PW500us SD500ms - Sitting'; %file safe names only
subjectName = 'LNP02'; setName = 'RCs'; setDescrpt = 'Set 1'; %file safe names only
externalStimulator = 'no'; %no = NanoStims
multipolar = 'no';
chanPort = 0; % 128 or 256 based on which Grapevine Port (B or C)
artifactBuffer = 0.005; % post stim artifact settle lag in s, so 0.01 is 10ms
rmsWindow = 0.05; % post stim RMS Window in s, so 0.05 is 50ms


% Muscles Used
Muscle_ID = {'Right VM', 'Right RF', 'Right VL', 'Right BF', 'Right ST', ...
    'Right TA', 'Right SO', 'Right LG','Left VM', 'Left RF', 'Left VL',...
    'Left BF', 'Left ST', 'Left TA', 'Left MG', 'Left LG'};

% Muscle_ID = {'Right VM', 'Right RF', 'Right VL', 'Right BF', 'Right ST', ...
%     'Right TA', 'Right SO', 'Right LG','Left VM', 'Left RF', 'Left VL',...
%     'Left BF', 'Left ST', 'Left TA', 'Left MG', 'Left LG',...
%     'Right Glute', 'Left Glute', 'Right Add', 'Left Add', 'Right TFL',...
%     'Left TFL', 'nan', 'nan','Left SM'};

% Muscle_ID = {'Right Add', 'Right RF', 'Right VL', 'Right BF', 'Right ST', ...
%     'Right TA', 'Right SO', 'Right LG', 'Left Add', 'Left RF', 'Left VL',...
%     'Left BF', 'Left ST', 'Left Ham', 'Left SO', 'Left LG'};

C = linspecer(length(Muscle_ID)/2,'qualitative');
C = [C;C];

reportPath = ['C:\figs\UH3\' subjectName '\emgRecruitment_Summary\']; %['R:\data_generated\human\uh3_stim\' subjectName '\emgRecruitment_Summary\'];
setPath = [reportPath setName '\'];
mkdir(setPath); %mkdir(reportPath);

%% Identify files to load
disp('Please Select the Ripple Data Folder');
[emgFilenames, emgPathname] = uigetfile([path_datatank subjectName '\*.nf6'],'Pick files','MultiSelect', 'on');
disp(['User selected ', fullfile(emgPathname, emgFilenames)]);

% disp('Please Select the Folder with the Stim Trial Info:');
trialPathname = uigetdir([path_stimparams subjectName '\'] , 'OpenLoop Stim Trial Info');

for f = 1:length(emgFilenames)
    trialFilenames{f} = [emgFilenames{f}(1:end-10) '.mat']; %2014a friendly, no erase function (-4 for just file, -10 for _x# thing)
end

%% Initialize Analysis Parameters

% Filtering Settings:
fs = 7500;
lowCut = 75; %lowest frequency to pass
highCut = 750; %highest frequency to pass
Norder = 2;
Wp = [lowCut, highCut]/(.5*fs);
[b,a]=butter(Norder, Wp);

%Smoothing Settings:
Fsmooth = 10; %Lowpass cutoff frequency for smoothing the rectified signal 10 vs 50
NorderS =  2;
dt = 1/fs;
[bS,aS] = butter(NorderS,Fsmooth/(0.5*fs));



%% Begin Parsing each file:
for fnum = 1:length(emgFilenames)
    
    disp(['Loading file ' num2str(fnum) ' of ' num2str(length(emgFilenames))]);
    
    %% Get Basic File Info
    trial(fnum).Name = erase(trialFilenames{fnum}, '.mat'); 
    [ns_status, hFile] = ns_OpenFile([emgPathname emgFilenames{fnum}]); 
    [ns_RESULT, nsFileInfo] = ns_GetFileInfo(hFile);
    ns_CloseFile(hFile);
    
    trial(fnum).Date = join(string([nsFileInfo.Time_Month nsFileInfo.Time_Day nsFileInfo.Time_Year]),'-');
    trial(fnum).Time = join(string([nsFileInfo.Time_Hour nsFileInfo.Time_Min nsFileInfo.Time_Sec nsFileInfo.Time_MilliSec]),':');
    
    %% Load Data
    [analogData,timeVec] = read_continuousData([emgPathname emgFilenames{fnum}], 'hifreq' , chanPort+[1:length(Muscle_ID)*2]); %128 vs 256 | 8 vs. 16
    %% 'Data' -> Bipolar EMG
    for i = 1:2:size(analogData,1)
        bipolar = analogData(i:i+1,:);
        emg(round(i/2),:) = diff(flipud(bipolar));
    end
    %% Filtering (Denoising) & Smoothing
    emgF = filtfilt(b,a,emg')';  
    emgS = filtfilt(bS,aS,emgF')';
    
    trial(fnum).emg = emg;
    trial(fnum).emgF = emgF;
    trial(fnum).timeVec = timeVec;
    disp('Filtering and Epoching....')
    
    %% Load Trial Stim Parameters
    trialInfo = load(fullfile(trialPathname,trialFilenames{fnum}));
        % Identify Stim Events Channel
        switch multipolar
            case 'yes' 
                spinalElec(fnum) = trialInfo.Stim_params(1).SpinalElecs{1}(1);
                anodeElec(fnum,:) = trialInfo.Stim_params(1).SpinalElecs{1}(2:end);
                        stimAmp(fnum) = trialInfo.Stim_params(1).SpinalElecAmps{1}(1);
            otherwise
                spinalElec(fnum) = cell2mat(trialInfo.Stim_params(1).SpinalElecs);
                        stimAmp(fnum) = cell2mat(trialInfo.Stim_params(1).SpinalElecAmps);
                        anodeElec = [];
        end
        nevStimCh(fnum) = trialInfo.Stim_params(1).NSChannels{1,1}(1);   
        % Set Pulse Width, Stim Duration, & Stim Frequency
        pulseWidth(fnum) = cell2mat(trialInfo.Stim_params(1).PulseWidth); %provided in ms
        stimDuration(fnum) = cell2mat(trialInfo.Stim_params(1).Duration)/1000; %in s, provided in ms
        stimFrequency(fnum) = cell2mat(trialInfo.Stim_params(1).Frequency); %this is provided in Hz
        numStims(fnum) = stimDuration(fnum)*stimFrequency(fnum);
    

    %% Load Stim Times
    switch externalStimulator
        case 'yes'
            [stimEvts] = read_digitalEvents([emgPathname emgFilenames{fnum}],1);
            stims = floor(cell2mat(stimEvts.timeStamp)*fs);
        otherwise
            [stimEvts,stchannels] = read_stimEvents([emgPathname emgFilenames{fnum}],nevStimCh(fnum));
            stims = floor(cell2mat(stimEvts)*fs);
    end
    
    stimLength = (pulseWidth(fnum)/1000)*fs + 2; %Pulse Width in ms, interstim interval is 2 ticks of clock (2/fs)
    trial(fnum).stimLength = stimLength;
    
    if length(stims) == 2*numStims(fnum)
        stimBegins = stims(1:2:end);
        stimEnds = stims(2:2:end);
    else
        stimBegins = stims;
    end
    
    trial(fnum).stims = stims;
    %% Build Epoch Time Vector
    
    % Set Epoch Window for Stim-triggered Averaging
% %     if stimFrequency(fnum) < 10
    dtPre = 50e-3; % msec before stim onset 
    dtPost = 100e-3; % msec after stim onset 
    nSampsPre = floor(dtPre*fs);
    nSampsPost = floor(dtPost*fs);
    
    trial(fnum).epochTimeVec = linspace(-dtPre, (stimLength/fs) + dtPost, nSampsPre + nSampsPost + stimLength);
   
    %% Epoch data into Windows 
    for i = 1:size(emg,1)
        asd = emgF(i,:);
        emgStim{i} = cell2mat(arrayfun(@(x) asd((stimBegins(x)-nSampsPre):(stimBegins(x)+(stimLength+nSampsPost-1))), 1:length(stimBegins), 'UniformOutput',false)');
        qwe = emgS(i,:);
        emgStim_Sm{i} = cell2mat(arrayfun(@(x) asd((stimBegins(x)-nSampsPre):(stimBegins(x)+(stimLength+nSampsPost-1))), 1:length(stimBegins), 'UniformOutput',false)');
    end
    
    
    trial(fnum).emgStim = emgStim;
    trial(fnum).emgStim_Sm = emgStim_Sm;
    
    %% Get RMS acros epochs
    trial(fnum).epochRMS = cellfun(@rms,emgStim,'UniformOutput', false);
    trial(fnum).epochP2P = cellfun(@peak2peak,emgStim,'UniformOutput', false);
    %% Extract "Mean" Trace & Response
    stimEnd = nSampsPre + stimLength + artifactBuffer*fs;    
    for m = 1:length(Muscle_ID)
        % Mean Trace
        trial(fnum).meanTrace{m} = mean(emgStim{1,m});
        trial(fnum).meanTrace_Sm{m} = mean(emgStim_Sm{1,m});
        trial(fnum).baseline(m,:) = mean(trial(fnum).epochRMS{m}(1:nSampsPre));
        % Response
        trial(fnum).responseRMS(m,:) = mean(trial(fnum).epochRMS{m}(stimEnd:stimEnd+(rmsWindow*fs)));
        trial(fnum).responseP2P(m,:) = mean(trial(fnum).epochP2P{m}(stimEnd:stimEnd+(rmsWindow*fs)));
    end

% %     %% Plot RMS Window for Verification
% %     disp('Plotting RMS Window for Verification...')
% %     
% %     hF2 = figure('visible', 'off');
% %         for c = 1:length(Muscle_ID)
% % 
% %          %%Plot /Channel
% %          plot(trial(fnum).epochTimeVec,trial(fnum).epochRMS{c},'Color',C(c,:));hold on;
% %               vline([0 stimLength/fs],'r-'); %vline(stimEnd/fs,'g-');
% % 
% %         end
% %         ylabel('RMS Amplitude (uV)');
% %         xlabel('Time (Sec)');
% %         title([trial(fnum).Name ' - Epoched RMS'],'Interpreter', 'none');
% %         vline([artifactBuffer artifactBuffer+rmsWindow], '-g');
% %         legend(Muscle_ID);    
% %     saveas(hF2,[setPath 'RMS_window-' trial(fnum).Name ' - ' num2str(stimAmp(fnum)/1000) 'mA_' ...
% %             num2str(pulseWidth(fnum)) 'us_' num2str(stimFrequency(fnum)) 'Hz_' num2str(stimDuration(fnum)) 's.png']);
% %  

    %% clean variables for later plotting.
    
        clearvars -except emgFilenames emgPathname reportPath Muscle_ID fs b a bS aS ...
                dtPre dtPost fnum spinalElec anodeElec nevStimCh stimAmp pulseWidth...
                stimDuration stimFrequency responseRMS trialName ...
                trialPathname trialFilenames trial subjectName setName ...
                setDescrpt setPath artifactBuffer rmsWindow ...
                externalStimulator nSampsPre nSampsPost chanPort C iter multipolar

    close all
    
    disp('End File Processing');
end
%Save Workspace for Later%

save([reportPath '\' setName '_' setDescrpt '.mat'], '-v7.3');
%% Export Summary Table:
disp('Writing Summary Table');

varNames = {'TrialName','SpinalElectrodes', 'NEV_StimChannel', 'StimAmplitudes_mA', ...
    'PulseWidths_ms', 'StimDurations_s', 'StimFrequencies_Hz', 'Record_Date', 'Record_Time'};

TrialNames = [{trial(:).Name}'];
TrialDates = [{trial(:).Date}'];
TrialTimes = [{trial(:).Time}'];

TrialTable = table(TrialNames, spinalElec', nevStimCh', (stimAmp/1000)', pulseWidth',...
    stimDuration', stimFrequency', TrialDates, TrialTimes, 'VariableNames', varNames)

 
writetable(TrialTable,[reportPath '\' strrep(datestr(now),':','_') '_' setName '_' setDescrpt '_Summary.txt'],'Delimiter','tab');
writetable(TrialTable,[reportPath '\' strrep(datestr(now),':','_') '_' setName '_' setDescrpt '_Summary.xls']);

% % %% Plotting Stim Averages
elecs = unique(spinalElec);
for ii = 1:length(elecs)
    eIdx = find(spinalElec == elecs(ii));
    [sortedStim,I] = sort(stimAmp(eIdx));
    [sortedPulse,I2] = sort(pulseWidth(eIdx));
    [sortedFreq,I3] = sort(stimFrequency(eIdx));

    if length(stimAmp(eIdx)) > 10
        iter = floor(length(stimAmp(eIdx))/10);
    else 
        iter = 2;
    end
    
    for c = 1:length(Muscle_ID)
        disp('Plotting Stim Averages');

        hF(c) = figure; %('visible', 'off');
        if ~isempty(anodeElec)
            tmpsup = suptitle([Muscle_ID{c} ' - E' num2str(elecs(ii)) '::' num2str(anodeElec(eIdx(1))) ' -  Stim Traces']);
        else
            tmpsup = suptitle([Muscle_ID{c} ' - E' num2str(elecs(ii)) ' -  Stim Traces']);
        end
        set(tmpsup, 'Interpreter', 'none') 
        maximize;

        subCount = 1; 
        for fnum = 1:iter:length(eIdx)
            hFs(fnum) = subplot(round(length(eIdx)/iter)+1,1,subCount);
    %         hFs(fnum) = subplot(10,1,subCount);
            %Plot /Channel
            plot(trial(eIdx(I(fnum))).epochTimeVec,trial(eIdx(I(fnum))).emgStim{c});hold on;
            plot(trial(eIdx(I(fnum))).epochTimeVec,trial(eIdx(I(fnum))).meanTrace{c},'k','LineWidth',1.2);
            vline([0    trial(eIdx(I(fnum))).stimLength/fs],'r-');
            tmptit = title([trial(eIdx(I(fnum))).Name ' - ' num2str(stimAmp(eIdx(I(fnum)))/1000) 'mA, ' ...
                num2str(pulseWidth(eIdx(I(fnum)))) 'us, ' num2str(stimFrequency(eIdx(I(fnum)))) 'Hz, ' num2str(stimDuration(eIdx(I(fnum)))) 's']); 
            set(hFs(fnum),'FontSize',7);
            set(tmptit, 'Interpreter', 'none','FontSize',5);
            hFs(fnum).XLim = [trial(eIdx(I(fnum))).epochTimeVec(1) trial(eIdx(I(fnum))).epochTimeVec(end)];
       
            subCount = subCount + 1;
        end
        linkaxes([hFs],'xy');
        suplabel('Amplitude(uV)','y');
        suplabel('Time(sec)','x');
        set(hF(c),'Position', [1343 49 570 947]);

        disp('Saving to Set Folder')
        if ~isempty(anodeElec)
            savefig(hF(c),[setPath Muscle_ID{c} '_EpochedStimTraces_' setName '_ E' num2str(elecs(ii)) '-' num2str(anodeElec(eIdx(1)))]);
            saveas(hF(c),[setPath Muscle_ID{c} '_EpochedStimTraces_' setName '_ E'  num2str(elecs(ii)) '-' num2str(anodeElec(eIdx(1))) '.png']);
        else
            savefig(hF(c),[setPath Muscle_ID{c} '_EpochedStimTraces_' setName '_ E' num2str(elecs(ii))]);
            saveas(hF(c),[setPath Muscle_ID{c} '_EpochedStimTraces_' setName '_ E'  num2str(elecs(ii)) '.png']);
        end
        
    end

    close all
    
    

    %% Plotting Recruitment Curves
% %     response = [trial(:).responseRMS]';
    response = [trial(:).responseP2P]';

    switch length(Muscle_ID)
        case 16
            pRow = 2; pCol = 1;
            mCount = length(Muscle_ID)/2;
            pIdx = 2;
        case 8
            pRow = 1; pCol = 3;
            mCount = length(Muscle_ID);
            pIdx = 1;
        otherwise
            disp('How Would you like to Plot the Recruitment Curves?');
            return;
    end

    hF3 = figure; %maximize;
    for pI = 1:pIdx

        hF3s(pI) = subplot(pRow,pCol, pCol*(pI-1)+1);
        for m = 1+mCount*(pI-1):mCount+mCount*(pI-1)
            yyA = smooth(sortedStim/1000,response(eIdx(I),m),'sgolay'); %/max(response(:,m))
            plot(sortedStim/1000,yyA,'Color',C(m,:),'LineWidth',3);
    %         scatter(stimAmp/1000, response(I,m)');
            hold on;
        end
            ylabel('(post-Stim) Mean P2P  (uV)'); %'(post-Stim) Mean RMS  (uV)'
            xlabel('Stimulation Amplitude (mA)');
            title(['E' num2str(elecs(ii)) ' - By Amplitude']);

% %          hF3s(pI) = subplot(pRow,pCol, pCol*(pI-1)+1);
% %         for m = 1+mCount*(pI-1):mCount+mCount*(pI-1)
% %     %         scatter(pulseWidth,response(:,m)')
% %             yyP = smooth(sortedPulse,response(I2,m),'sgolay');
% %             plot(sortedPulse,yyP,'Color',C(m,:),'LineWidth',3);
% %             hold on;
% %         end
% %             ylabel('(post-Stim) Mean RMS  (uV)')
% %             xlabel('Stimulation PulseWidth (us)')
% %             title('By PulseWidth');
    %     %     
    %     subplot(pRow,pCol, pCol*(pI-1)+3)
    %     for m = 1+mCount*(pI-1):mCount+mCount*(pI-1)
    % %         scatter(stimFrequency,response(:,m)')
    %         yyF = smooth(sortedFreq,response(I3,m),'sgolay');
    %         plot(sortedFreq,yyF,'Color',C(m,:),'LineWidth',3);
    %         hold on;
    %     end
    %         ylabel('(post-Stim) Mean RMS  (uV)')
    %         xlabel('Stimulation Frequency (Hz)')
    %         title('By Frequency');

        if pI == 1
            legend([Muscle_ID(1:mCount)],'Location','northwest', 'Orientation','vertical');
        elseif pI == 2
            legend([Muscle_ID((mCount)+1:length(Muscle_ID))],'Location','northwest', 'Orientation','vertical');
        end
    %     
    end
    if ~isempty(anodeElec)
        tmptit = sgtitle([setName ' - ' setDescrpt ' - e' num2str(elecs(ii)) '::' num2str(anodeElec(eIdx(1))) ' -  Recruitment Curves (P2P)']);
    else
        tmptit = sgtitle([setName ' - ' setDescrpt ' - e' num2str(elecs(ii)) ' -  Recruitment Curves (P2P)']);
    end
    set(tmptit, 'Interpreter', 'none');
    % % Pix_SS = get(0,'screensize');
    hF3.Position = [680 49 572 947];
%     linkaxes([hF3s],'xy');

    if ~isempty(anodeElec)
        saveas(hF3,[reportPath '\' strrep(datestr(now),':','_')  '_RecruitmentCurvesP2P' '_' setName '_ E' num2str(elecs(ii)) '-' num2str(anodeElec(eIdx(1))) '.png']);
        savefig(hF3,[reportPath '\' strrep(datestr(now),':','_') '_RecruitmentCurvesP2P_' setName '_ E' num2str(elecs(ii)) '-' num2str(anodeElec(eIdx(1))) ]);
    else
        saveas(hF3,[reportPath '\' strrep(datestr(now),':','_')  '_RecruitmentCurvesP2P' '_' setName '_ E' num2str(elecs(ii)) '.png']);
        savefig(hF3,[reportPath '\' strrep(datestr(now),':','_') '_RecruitmentCurvesP2P_' setName '_ E' num2str(elecs(ii)) ]);
    end
    
    
    %% Plot Normalized RC
    hF4 = figure; %maximize;
    for pI = 1:pIdx

        hF4s(pI) = subplot(pRow,pCol, pCol*(pI-1)+1);
        for m = 1+mCount*(pI-1):mCount+mCount*(pI-1)
            yyA = smooth(sortedStim/1000,response(eIdx(I),m)/max(response(:,m)),'sgolay'); %/max(response(:,m))
            plot(sortedStim/1000,yyA,'Color',C(m,:),'LineWidth',3);
    %         scatter(stimAmp/1000, response(I,m)');
            hold on;
        end
            ylabel('(post-Stim) Normalized RMS')
            xlabel('Stimulation Amplitude (mA)')
            title(['E' num2str(elecs(ii)) ' - By Amplitude']);

% %         hF3s(pI) = subplot(pRow,pCol, pCol*(pI-1)+1);
% %         for m = 1+mCount*(pI-1):mCount+mCount*(pI-1)
% %     %         scatter(pulseWidth,response(:,m)')
% %             yyP = smooth(sortedPulse,response(I2,m)/max(response(:,m)),'sgolay');
% %             plot(sortedPulse,yyP,'Color',C(m,:),'LineWidth',3);
% %             hold on;
% %         end
% %             ylabel('(post-Stim) Mean RMS  (uV)')
% %             xlabel('Stimulation PulseWidth (us)')
% %             title('By PulseWidth');
    %     %     
    %     subplot(pRow,pCol, pCol*(pI-1)+3)
    %     for m = 1+mCount*(pI-1):mCount+mCount*(pI-1)
    % %         scatter(stimFrequency,response(:,m)')
    %         yyF = smooth(sortedFreq,response(I3,m),'sgolay');
    %         plot(sortedFreq,yyF,'Color',C(m,:),'LineWidth',3);
    %         hold on;
    %     end
    %         ylabel('(post-Stim) Mean RMS  (uV)')
    %         xlabel('Stimulation Frequency (Hz)')
    %         title('By Frequency');

        if pI == 1
            legend([Muscle_ID(1:mCount)],'Location','northwest', 'Orientation','vertical');
        elseif pI == 2
            legend([Muscle_ID((mCount)+1:length(Muscle_ID))],'Location','northwest', 'Orientation','vertical');
        end
    %     
    end
    if ~isempty(anodeElec)
        tmptit = sgtitle([setName ' - ' setDescrpt ' - e' num2str(elecs(ii)) '::' num2str(anodeElec(eIdx(1))) ' - Norm Recruitment Curves (RMS)']);
    else
        tmptit = sgtitle([setName ' - ' setDescrpt ' - e' num2str(elecs(ii)) ' - Norm Recruitment Curves (RMS)']);
    end
    set(tmptit, 'Interpreter', 'none');
    % % Pix_SS = get(0,'screensize');
    hF4.Position = [680 49 572 947];
%     linkaxes([hF3s],'xy');

    if ~isempty(anodeElec)
        saveas(hF4,[reportPath '\' strrep(datestr(now),':','_')  '_NormRecruitmentCurvesRMS' '_' setName '_ E' num2str(elecs(ii)) '-' num2str(anodeElec(eIdx(1))) '.png']);
        savefig(hF4,[reportPath '\' strrep(datestr(now),':','_') '_NormRecruitmentCurvesRMS_' setName '_ E' '-' num2str(anodeElec(eIdx(1))) num2str(elecs(ii)) ]);
    else
        saveas(hF4,[reportPath '\' strrep(datestr(now),':','_')  '_NormRecruitmentCurvesRMS' '_' setName '_ E' num2str(elecs(ii)) '.png']);
        savefig(hF4,[reportPath '\' strrep(datestr(now),':','_') '_NormRecruitmentCurvesRMS_' setName '_ E' num2str(elecs(ii)) ]);
    end
    
    close all;
    
    
% %      hF5 = figure; %maximize;
% %     for pI = 1:pIdx
% %          hF3s(pI) = subplot(pRow,pCol, pCol*(pI-1)+1);
% %         for m = 1+mCount*(pI-1):mCount+mCount*(pI-1)
% %     %         scatter(pulseWidth,response(:,m)')
% %             yyP = smooth(sortedPulse,response(I2,m),'sgolay');
% %             plot(sortedPulse,yyP,'Color',C(m,:),'LineWidth',3);
% %             hold on;
% %         end
% %             ylabel('(post-Stim) Mean P2P  (uV)')
% %             xlabel('Stimulation PulseWidth (us)')
% %             title('RC By PulseWidth');
% %     %     %     
% % 
% %         if pI == 1
% %             legend([Muscle_ID(1:mCount)],'Location','northwest', 'Orientation','vertical');
% %         elseif pI == 2
% %             legend([Muscle_ID((mCount)+1:length(Muscle_ID))],'Location','northwest', 'Orientation','vertical');
% %         end
% %     %     
% %     end
% %     if ~isempty(anodeElec)
% %         tmptit = sgtitle([setName ' - ' setDescrpt ' - e' num2str(elecs(ii)) '::' num2str(anodeElec(eIdx(1))) ' -  Recruitment Curves (P2P)']);
% %     else
% %         tmptit = sgtitle([setName ' - ' setDescrpt ' - e' num2str(elecs(ii)) ' -  Recruitment Curves (P2P)']);
% %     end
% %     set(tmptit, 'Interpreter', 'none');
% %     % % Pix_SS = get(0,'screensize');
% %     hF5.Position = [680 49 572 947];
% % %     linkaxes([hF3s],'xy');
% % 
% %     if ~isempty(anodeElec)
% %         saveas(hF5,[reportPath '\' strrep(datestr(now),':','_')  '_RecruitmentCurvesPW' '_' setName '_ E' num2str(elecs(ii)) '-' num2str(anodeElec(eIdx(1))) '.png']);
% %         savefig(hF5,[reportPath '\' strrep(datestr(now),':','_') '_RecruitmentCurvesPW_' setName '_ E' num2str(elecs(ii)) '-' num2str(anodeElec(eIdx(1))) ]);
% %     else
% %         saveas(hF5,[reportPath '\' strrep(datestr(now),':','_')  '_RecruitmentCurvesPW' '_' setName '_ E' num2str(elecs(ii)) '.png']);
% %         savefig(hF5,[reportPath '\' strrep(datestr(now),':','_') '_RecruitmentCurvesPW_' setName '_ E' num2str(elecs(ii)) ]);
% %     end

%% Plotting baseline subtracted RC

    hF6 = figure; %maximize;
    for pI = 1:pIdx

        hF6s(pI) = subplot(pRow,pCol, pCol*(pI-1)+1);
        for m = 1+mCount*(pI-1):mCount+mCount*(pI-1)
            
            for ll =  1:length(trial)
                X(ll,:) = trial(ll).baseline(m);
            end
            mBase = mean(X);
            
            yyA = smooth(sortedStim/1000,response(eIdx(I),m)-mBase,'sgolay'); %/max(response(:,m))
            plot(sortedStim/1000,yyA,'Color',C(m,:),'LineWidth',3);
    %         scatter(stimAmp/1000, response(I,m)');
            hold on;
        end
            ylabel('(post-Stim) basesub-Mean RMS  (uV)'); %'(post-Stim) Mean RMS  (uV)'
            xlabel('Stimulation Amplitude (mA)');
            title(['E' num2str(elecs(ii)) ' - By Amplitude']);


        if pI == 1
            legend([Muscle_ID(1:mCount)],'Location','northwest', 'Orientation','vertical');
        elseif pI == 2
            legend([Muscle_ID((mCount)+1:length(Muscle_ID))],'Location','northwest', 'Orientation','vertical');
        end
    %     
    end
    if ~isempty(anodeElec)
        tmptit = sgtitle([setName ' - ' setDescrpt ' - e' num2str(elecs(ii)) '::' num2str(anodeElec(eIdx(1))) ' -  Recruitment Curves (RMS)']);
    else
        tmptit = sgtitle([setName ' - ' setDescrpt ' - e' num2str(elecs(ii)) ' -  Recruitment Curves (RMS)']);
    end
    set(tmptit, 'Interpreter', 'none');
    % % Pix_SS = get(0,'screensize');
    hF6.Position = [680 49 572 947];
%     linkaxes([hF3s],'xy');

    if ~isempty(anodeElec)
        saveas(hF6,[reportPath '\' strrep(datestr(now),':','_')  '_bsrRecruitmentCurvesRMS' '_' setName '_ E' num2str(elecs(ii)) '-' num2str(anodeElec(eIdx(1))) '.png']);
        savefig(hF6,[reportPath '\' strrep(datestr(now),':','_') '_bsrRecruitmentCurvesRMS_' setName '_ E' num2str(elecs(ii)) '-' num2str(anodeElec(eIdx(1))) ]);
    else
        saveas(hF6,[reportPath '\' strrep(datestr(now),':','_')  '_bsrRecruitmentCurvesRMS' '_' setName '_ E' num2str(elecs(ii)) '.png']);
        savefig(hF6,[reportPath '\' strrep(datestr(now),':','_') '_bsrRecruitmentCurvesRMS_' setName '_ E' num2str(elecs(ii)) ]);
    end
 
end
close all

 


    