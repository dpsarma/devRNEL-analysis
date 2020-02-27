
%% Experiment Session Parameters

subjectName = 'LSP02b';
elecNames = [1 2 3 4 5 6 7 8 25 26 27 28 29 30 31 32 33 34 35 36];
runType = 'RC'; %or 'RC';

switch runType
    case 'LC'
        setDescrpt = 'Limit Check f-2Hz PW-1ms';
    case 'LC2'
        setDescrpt = 'Limit Check f-2Hz PW-200us';  
    case 'RC'
        setDescrpt = 'RC 500uA-XmA f-1Hz PW-200us';
    otherwise
        disp('Wrong Trial Case');
        return;
end

path_UH3 = 'R:\data_raw\human\uh3_stim\';
addpath('C:\Users\dsarma\Documents\GitHub\mattools\linspecer\');

chanPort = 128; % 128 or 256 based on which Grapevine Port (B or C)
artifactBuffer = 0.005; % post stim artifact settle lag in s, so 0.01 is 10ms
rmsWindow = 0.05; % post stim RMS Window in s, so 0.05 is 50ms
iter = 1;

Muscle_ID = {'Right VM', 'Right RF', 'Right VL', 'Right BF', 'Right ST', ...
    'Right TA', 'Right SO', 'Right LG', 'Left VM', 'Left RF', 'Left VL',...
    'Left BF', 'Left ST', 'Left Ham', 'Left SO', 'Left LG'};

C = linspecer(length(Muscle_ID),'qualitative');

reportPath = ['R:\data_generated\human\uh3_stim\' subjectName '\emgRecruitment_Summary\'];


%% Identify files to load
disp('Please Select the Ripple Data Folder');
[emgFilenames, emgPathname] = uigetfile([path_UH3 subjectName '\data\Open Loop\Trellis\*.nev'],'Pick files','MultiSelect', 'on');
disp(['User selected ', fullfile(emgPathname, emgFilenames)]);

% disp('Please Select the Folder with the Stim Trial Info:');
trialPathname = uigetdir([path_UH3 subjectName '\'] , 'OpenLoop Stim Trial Info');

for f = 1:length(emgFilenames)
    trialFilenames{f} = [emgFilenames{f}(1:end-10) '.mat']; %2014a friendly, no erase function (-4 for just file, -10 for _x# thing)
end


%% Initialize Analysis Parameters

% Filtering Settings:
fs = 30e3;
lowCut = 75; %lowest frequency to pass
highCut = 7500; %highest frequency to pass
Norder = 2;
Wp = [lowCut, highCut]/(.5*fs);
[b,a]=butter(Norder, Wp);

% Smoothing Settings:
Fsmooth = 10; %Lowpass cutoff frequency for smoothing the rectified signal 10 vs 50
NorderS =  2;
dt = 1/fs;
[bS,aS] = butter(NorderS,Fsmooth/(0.5*fs));

% Set Epoch Window for Stim-triggered Averaging
dtPre = 50e-3; % msec before stim onset 
dtPost = 100e-3; % msec after stim onset 
nSampsPre = floor(dtPre*fs);
nSampsPost = floor(dtPost*fs);


%% Perform Analysis/Electrode
for eIdx = 1:length(elecNames)
    setName = ['Elec ' char(string(elecNames(eIdx)))];
    setPath = [reportPath setName '\'];
    mkdir(setPath); %mkdir(reportPath);
    disp(['Evaluating ' runType ' for ' setName]);
    
    switch setName
        case 'Elec 1'
            if strcmp(runType,'LC')
                fIndices = find(contains(trialFilenames,'Set002'));
            elseif strcmp(runType,'RC')
                fIndices = find(contains(trialFilenames,'Set003'));
            else
                disp(['No Files Found for' setName]);
            end
        case 'Elec 2'
            if strcmp(runType,'LC')
                fIndices = find(contains(trialFilenames,'Set004'));
            elseif strcmp(runType,'RC')
                fIndices = find(contains(trialFilenames,'Set005'));
            else
                disp(['No Files Found for' setName]);
            end
        case 'Elec 3'
            if strcmp(runType,'LC')
                fIndices = find(contains(trialFilenames,'Set006'));
            elseif strcmp(runType,'RC')
                fIndices = find(contains(trialFilenames,'Set007')|contains(trialFilenames,'Set008'));
            else
                disp(['No Files Found for' setName]);
            end
        case 'Elec 4'
            if strcmp(runType,'LC')
                fIndices = find(contains(trialFilenames,'Set009'));
            elseif strcmp(runType,'RC')
                fIndices = find(contains(trialFilenames,'Set010'));
            else
                disp(['No Files Found for' setName]);
            end
        case 'Elec 5'
            if strcmp(runType,'LC')
                fIndices = find(contains(trialFilenames,'Set011'));
            elseif strcmp(runType,'RC')
                fIndices = find(contains(trialFilenames,'Set012')|contains(trialFilenames,'Set013')|contains(trialFilenames,'Set014'));
            else
                disp(['No Files Found for' setName]);
            end
        case 'Elec 6'
            if strcmp(runType,'LC')
                fIndices = find(contains(trialFilenames,'Set015'));
            elseif strcmp(runType,'RC')
                fIndices = find(contains(trialFilenames,'Set016')|contains(trialFilenames,'Set017'));
            else
                disp(['No Files Found for' setName]);
            end
        case 'Elec 7'
            if strcmp(runType,'LC2')
                fIndices = find(contains(trialFilenames,'Set018'));
            elseif strcmp(runType,'LC')
                fIndices = find(contains(trialFilenames,'Set019'));
            elseif strcmp(runType,'RC')
                fIndices = find(contains(trialFilenames,'Set020')|contains(trialFilenames,'Set021'));
            else
                disp(['No Files Found for' setName]);
            end
        case 'Elec 8'
            if strcmp(runType,'LC2')
                fIndices = find(contains(trialFilenames,'Set022'));
            elseif strcmp(runType,'LC')
                fIndices = find(contains(trialFilenames,'Set023'));
            elseif strcmp(runType,'RC')
                fIndices = find(contains(trialFilenames,'Set024')|contains(trialFilenames,'Set025')|contains(trialFilenames,'Set026'));
            else
                disp(['No Files Found for' setName]);
            end
        case 'Elec 25'
            if strcmp(runType,'LC')
                fIndices = find(contains(trialFilenames,'Set029'));
            elseif strcmp(runType,'RC')
                fIndices = find(contains(trialFilenames,'Set030'));
            else
                disp(['No Files Found for' setName]);
            end
                
        case 'Elec 26'
            if strcmp(runType,'LC')
                fIndices = find(contains(trialFilenames,'Set031'));
            elseif strcmp(runType,'RC')
                fIndices = find(contains(trialFilenames,'Set032'));
            else
                disp(['No Files Found for' setName]);
            end
                
        case 'Elec 27'
            if strcmp(runType,'LC')
                fIndices = find(contains(trialFilenames,'Set033'));
            elseif strcmp(runType,'RC')
                fIndices = find(contains(trialFilenames,'Set034'));
            else
                disp(['No Files Found for' setName]);
            end
                
        case 'Elec 28'
            if strcmp(runType,'LC')
                fIndices = find(contains(trialFilenames,'Set035'));
            elseif strcmp(runType,'RC')
                fIndices = find(contains(trialFilenames,'Set036'));
            else
                disp(['No Files Found for' setName]);
            end
                
        case 'Elec 29'
            if strcmp(runType,'LC')
                fIndices = find(contains(trialFilenames,'Set037'));
            elseif strcmp(runType,'RC')
                fIndices = find(contains(trialFilenames,'Set038'));
            else
                disp(['No Files Found for' setName]);
            end
                
        case 'Elec 30'
            if strcmp(runType,'LC')
                fIndices = find(contains(trialFilenames,'Set039'));
            elseif strcmp(runType,'RC')
                fIndices = find(contains(trialFilenames,'Set040'));
            else
                disp(['No Files Found for' setName]);
            end
                
        case 'Elec 31'
            if strcmp(runType,'LC')
                fIndices = find(contains(trialFilenames,'Set041'));
            elseif strcmp(runType,'RC')
                fIndices = find(contains(trialFilenames,'Set042'));
            else
                disp(['No Files Found for' setName]);
            end
                
        case 'Elec 32'
            if strcmp(runType,'LC')
                fIndices = find(contains(trialFilenames,'Set043'));
            elseif strcmp(runType,'RC')
                fIndices = find(contains(trialFilenames,'Set044')|contains(trialFilenames,'Set045'));
            else
                disp(['No Files Found for' setName]);
            end
                
        case 'Elec 33'
            if strcmp(runType,'LC')
                fIndices = find(contains(trialFilenames,'Set046'));
            elseif strcmp(runType,'RC')
                fIndices = find(contains(trialFilenames,'Set047')|contains(trialFilenames,'Set048'));
            else
                disp(['No Files Found for' setName]);
            end
                
        case 'Elec 34'
            if strcmp(runType,'LC')
                fIndices = find(contains(trialFilenames,'Set049'));
            elseif strcmp(runType,'RC')
                fIndices = find(contains(trialFilenames,'Set050'));
            else
                disp(['No Files Found for' setName]);
            end
                
        case 'Elec 35'
            if strcmp(runType,'LC')
                fIndices = find(contains(trialFilenames,'Set051'));
            elseif strcmp(runType,'RC')
                fIndices = find(contains(trialFilenames,'Set052'));
            else
                disp(['No Files Found for' setName]);
            end
                
        case 'Elec 36'
            if strcmp(runType,'LC')
                fIndices = find(contains(trialFilenames,'Set053'));
            elseif strcmp(runType,'RC')
                fIndices = find(contains(trialFilenames,'Set054'));
            else
                disp(['No Files Found for' setName]);
            end
        case 'Elec 37'
        case 'Elec 38'
        case 'Elec 39'
        case 'Elec 40'
        otherwise
                disp('TWAT ARE YOU THINKING?!?!?!?!');
    end
    
    if exist('fIndices')
        disp(fIndices);
    else 
        continue;
    end
    %% Parse Files:
    for fnum = 1:length(fIndices)
        disp(['Loading file ' num2str(fnum) ' of ' num2str(length(fIndices))]);
    
        %% Get Basic File Info
        trial(fIndices(fnum)).Name = erase(trialFilenames{fIndices(fnum)}, '.mat'); 
        [ns_status, hFile] = ns_OpenFile([emgPathname emgFilenames{fIndices(fnum)}]); 
        [ns_RESULT, nsFileInfo] = ns_GetFileInfo(hFile);
        ns_CloseFile(hFile);

        trial(fIndices(fnum)).Date = join(string([nsFileInfo.Time_Month nsFileInfo.Time_Day nsFileInfo.Time_Year]),'-');
        trial(fIndices(fnum)).Time = join(string([nsFileInfo.Time_Hour nsFileInfo.Time_Min nsFileInfo.Time_Sec nsFileInfo.Time_MilliSec]),':');

        %% Load Data
        [analogData,timeVec] = read_continuousData([emgPathname emgFilenames{fIndices(fnum)}], 'raw' , chanPort+[1:length(Muscle_ID)*2]); %128 vs 256 | 8 vs. 16
        %% 'Data' -> Bipolar EMG
        for i = 1:2:size(analogData,1)
            bipolar = analogData(i:i+1,:);
            emg(round(i/2),:) = diff(flipud(bipolar));
        end
        %% Filtering (Denoising) & Smoothing
        emgF = filtfilt(b,a,emg')';  
        emgS = filtfilt(bS,aS,emgF')';

        disp('Filtering and Epoching....')

        %% Load Trial Stim Parameters
        trialInfo = load(fullfile(trialPathname,trialFilenames{fIndices(fnum)}));
            % Identify Stim Events Channel
            spinalElec(fIndices(fnum)) = cell2mat(trialInfo.Stim_params(1).SpinalElecs);
            nevStimCh(fIndices(fnum)) = trialInfo.Stim_params(1).NSChannels{1,1}(1);   
            % Set Pulse Width, Stim Duration, & Stim Frequency
            stimAmp(fIndices(fnum)) = cell2mat(trialInfo.Stim_params(1).SpinalElecAmps);
            pulseWidth(fIndices(fnum)) = cell2mat(trialInfo.Stim_params(1).PulseWidth); %provided in ms
            stimDuration(fIndices(fnum)) = cell2mat(trialInfo.Stim_params(1).Duration)/1000; %in s, provided in ms
            stimFrequency(fIndices(fnum)) = cell2mat(trialInfo.Stim_params(1).Frequency); %this is provided in Hz
            numStims(fIndices(fnum)) = stimDuration(fIndices(fnum))*stimFrequency(fIndices(fnum));


            [stimEvts,stchannels] = read_stimEvents([emgPathname emgFilenames{fIndices(fnum)}],nevStimCh(fIndices(fnum)));
            stims = floor(cell2mat(stimEvts)*fs);
       
        stimLength = (pulseWidth(fIndices(fnum))/1000)*fs + 2; %Pulse Width in ms, interstim interval is 2 ticks of clock (2/fs)
        trial(fIndices(fnum)).stimLength = stimLength;

        if length(stims) == 2*numStims(fIndices(fnum))
            stimBegins = stims(1:2:end);
            stimEnds = stims(2:2:end);
        else
            stimBegins = stims;
        end
        %% Build Epoch Time Vector

        trial(fIndices(fnum)).epochTimeVec = linspace(-dtPre, (stimLength/fs) + dtPost, nSampsPre + nSampsPost + stimLength);

        %% Epoch data into Windows 
        for i = 1:size(emg,1)
            asd = emgF(i,:);
            emgStim{i} = cell2mat(arrayfun(@(x) asd((stimBegins(x)-nSampsPre):(stimBegins(x)+(stimLength+nSampsPost-1))), 1:length(stimBegins), 'UniformOutput',false)');
            qwe = emgS(i,:);
            emgStim_Sm{i} = cell2mat(arrayfun(@(x) asd((stimBegins(x)-nSampsPre):(stimBegins(x)+(stimLength+nSampsPost-1))), 1:length(stimBegins), 'UniformOutput',false)');
        end

        trial(fIndices(fnum)).emgStim = emgStim;
        trial(fIndices(fnum)).emgStim_Sm = emgStim_Sm;

        %% Get RMS acros epochs
        trial(fIndices(fnum)).epochRMS = cellfun(@rms,emgStim,'UniformOutput', false);

        %% Extract "Mean" Trace & Response
        stimEnd = nSampsPre + stimLength + artifactBuffer*fs;    
        for m = 1:length(Muscle_ID)
            % Mean Trace
            trial(fIndices(fnum)).meanTrace{m} = mean(emgStim{1,m});
            trial(fIndices(fnum)).meanTrace_Sm{m} = mean(emgStim_Sm{1,m});
            % Response
            trial(fIndices(fnum)).responseRMS(m,:) = mean(trial(fIndices(fnum)).epochRMS{m}(stimEnd:stimEnd+(rmsWindow*fs)));
            trial(fIndices(fnum)).responseP2P(m,:) = mean(trial(fIndices(fnum)).emgStim{m}(stimEnd:stimEnd+(rmsWindow*fs)));
        end

        %% Plot RMS Window for Verification
        disp('Plotting RMS Window for Verification...')

        hF2 = figure('visible', 'off');
            for c = 1:length(Muscle_ID)

             %%Plot /Channel
             plot(trial(fIndices(fnum)).epochTimeVec,trial(fIndices(fnum)).epochRMS{c},'Color',C(c,:));hold on;
                  vline([0 stimLength/fs],'r-'); %vline(stimEnd/fs,'g-');

            end
            ylabel('RMS Amplitude (uV)');
            xlabel('Time (Sec)');
            title([trial(fIndices(fnum)).Name ' - Epoched RMS'],'Interpreter', 'none');
            vline([artifactBuffer artifactBuffer+rmsWindow], '-g');
            legend(Muscle_ID);    
        saveas(hF2,[setPath 'RMS_window-' trial(fIndices(fnum)).Name ' - ' num2str(stimAmp(fIndices(fnum))/1000) 'mA_' ...
                num2str(pulseWidth(fIndices(fnum))) 'ms_' num2str(stimFrequency(fIndices(fnum))) 'Hz_' num2str(stimDuration(fIndices(fnum))) 's.png']);

        %% clean variables for later plotting.

        clearvars -except emgFilenames emgPathname reportPath Muscle_ID fs b a bS aS ...
                    dtPre dtPost fnum spinalElec nevStimCh stimAmp pulseWidth...
                    stimDuration stimFrequency responseRMS trialName ...
                    trialPathname trialFilenames trial subjectName setName ...
                    setDescrpt setPath artifactBuffer rmsWindow ...
                    nSampsPre nSampsPost chanPort C iter eIdx ...
                    elecNames path_UH3 runType fIndices
                
        close all

        disp(['End File Processing for' setName ]);
    end
    
    clear fIndices
end

%% Export Summary Table:
disp('Writing Summary Table\n');

varNames = {'TrialName','SpinalElectrodes', 'NEV_StimChannel', 'StimAmplitudes_mA', ...
    'PulseWidths_ms', 'StimDurations_s', 'StimFrequencies_Hz', 'Record_Date', 'Record_Time'};

TrialNames = [{trial(:).Name}'];
TrialDates = [{trial(:).Date}'];
TrialTimes = [{trial(:).Time}'];

TrialTable = table(TrialNames, spinalElec', nevStimCh', (stimAmp/1000)', pulseWidth',...
    stimDuration', stimFrequency', TrialDates, TrialTimes, 'VariableNames', varNames)

 
writetable(TrialTable,[reportPath '\' strrep(datestr(now),':','_') '_' setDescrpt '_Summary.txt'],'Delimiter','tab');
writetable(TrialTable,[reportPath '\' strrep(datestr(now),':','_') '_' setDescrpt '_Summary.xls']);

disp('End File Grep, Begin Plotting Things....\n');

clear setName

%% Plotting Loop
for eIdx2 = 1:length(elecNames)
    setName = ['Elec ' char(string(elecNames(eIdx2)))];
    switch setName
        case 'Elec 1'
            if strcmp(runType,'LC')
                fIndices = find(contains(trialFilenames,'Set002'));
            elseif strcmp(runType,'RC')
                fIndices = find(contains(trialFilenames,'Set003'));
            else
                disp(['No Files Found for' setName]);
            end
        case 'Elec 2'
            if strcmp(runType,'LC')
                fIndices = find(contains(trialFilenames,'Set004'));
            elseif strcmp(runType,'RC')
                fIndices = find(contains(trialFilenames,'Set005'));
            else
                disp(['No Files Found for' setName]);
            end
        case 'Elec 3'
            if strcmp(runType,'LC')
                fIndices = find(contains(trialFilenames,'Set006'));
            elseif strcmp(runType,'RC')
                fIndices = find(contains(trialFilenames,'Set007')|contains(trialFilenames,'Set008'));
            else
                disp(['No Files Found for' setName]);
            end
        case 'Elec 4'
            if strcmp(runType,'LC')
                fIndices = find(contains(trialFilenames,'Set009'));
            elseif strcmp(runType,'RC')
                fIndices = find(contains(trialFilenames,'Set010'));
            else
                disp(['No Files Found for' setName]);
            end
        case 'Elec 5'
            if strcmp(runType,'LC')
                fIndices = find(contains(trialFilenames,'Set011'));
            elseif strcmp(runType,'RC')
                fIndices = find(contains(trialFilenames,'Set012')|contains(trialFilenames,'Set013')|contains(trialFilenames,'Set014'));
            else
                disp(['No Files Found for' setName]);
            end
        case 'Elec 6'
            if strcmp(runType,'LC')
                fIndices = find(contains(trialFilenames,'Set015'));
            elseif strcmp(runType,'RC')
                fIndices = find(contains(trialFilenames,'Set016')|contains(trialFilenames,'Set017'));
            else
                disp(['No Files Found for' setName]);
            end
        case 'Elec 7'
            if strcmp(runType,'LC2')
                fIndices = find(contains(trialFilenames,'Set018'));
            elseif strcmp(runType,'LC')
                fIndices = find(contains(trialFilenames,'Set019'));
            elseif strcmp(runType,'RC')
                fIndices = find(contains(trialFilenames,'Set020')|contains(trialFilenames,'Set021'));
            else
                disp(['No Files Found for' setName]);
            end
        case 'Elec 8'
            if strcmp(runType,'LC2')
                fIndices = find(contains(trialFilenames,'Set022'));
            elseif strcmp(runType,'LC')
                fIndices = find(contains(trialFilenames,'Set023'));
            elseif strcmp(runType,'RC')
                fIndices = find(contains(trialFilenames,'Set024')|contains(trialFilenames,'Set025')|contains(trialFilenames,'Set026'));
            else
                disp(['No Files Found for' setName]);
            end
        case 'Elec 25'
            if strcmp(runType,'LC')
                fIndices = find(contains(trialFilenames,'Set029'));
            elseif strcmp(runType,'RC')
                fIndices = find(contains(trialFilenames,'Set030'));
            else
                disp(['No Files Found for' setName]);
            end
                
        case 'Elec 26'
            if strcmp(runType,'LC')
                fIndices = find(contains(trialFilenames,'Set031'));
            elseif strcmp(runType,'RC')
                fIndices = find(contains(trialFilenames,'Set032'));
            else
                disp(['No Files Found for' setName]);
            end
                
        case 'Elec 27'
            if strcmp(runType,'LC')
                fIndices = find(contains(trialFilenames,'Set033'));
            elseif strcmp(runType,'RC')
                fIndices = find(contains(trialFilenames,'Set034'));
            else
                disp(['No Files Found for' setName]);
            end
                
        case 'Elec 28'
            if strcmp(runType,'LC')
                fIndices = find(contains(trialFilenames,'Set035'));
            elseif strcmp(runType,'RC')
                fIndices = find(contains(trialFilenames,'Set036'));
            else
                disp(['No Files Found for' setName]);
            end
                
        case 'Elec 29'
            if strcmp(runType,'LC')
                fIndices = find(contains(trialFilenames,'Set037'));
            elseif strcmp(runType,'RC')
                fIndices = find(contains(trialFilenames,'Set038'));
            else
                disp(['No Files Found for' setName]);
            end
                
        case 'Elec 30'
            if strcmp(runType,'LC')
                fIndices = find(contains(trialFilenames,'Set039'));
            elseif strcmp(runType,'RC')
                fIndices = find(contains(trialFilenames,'Set040'));
            else
                disp(['No Files Found for' setName]);
            end
                
        case 'Elec 31'
            if strcmp(runType,'LC')
                fIndices = find(contains(trialFilenames,'Set041'));
            elseif strcmp(runType,'RC')
                fIndices = find(contains(trialFilenames,'Set042'));
            else
                disp(['No Files Found for' setName]);
            end
                
        case 'Elec 32'
            if strcmp(runType,'LC')
                fIndices = find(contains(trialFilenames,'Set043'));
            elseif strcmp(runType,'RC')
                fIndices = find(contains(trialFilenames,'Set044')|contains(trialFilenames,'Set045'));
            else
                disp(['No Files Found for' setName]);
            end
                
        case 'Elec 33'
            if strcmp(runType,'LC')
                fIndices = find(contains(trialFilenames,'Set046'));
            elseif strcmp(runType,'RC')
                fIndices = find(contains(trialFilenames,'Set047')|contains(trialFilenames,'Set048'));
            else
                disp(['No Files Found for' setName]);
            end
                
        case 'Elec 34'
            if strcmp(runType,'LC')
                fIndices = find(contains(trialFilenames,'Set049'));
            elseif strcmp(runType,'RC')
                fIndices = find(contains(trialFilenames,'Set050'));
            else
                disp(['No Files Found for' setName]);
            end
                
        case 'Elec 35'
            if strcmp(runType,'LC')
                fIndices = find(contains(trialFilenames,'Set051'));
            elseif strcmp(runType,'RC')
                fIndices = find(contains(trialFilenames,'Set052'));
            else
                disp(['No Files Found for' setName]);
            end
                
        case 'Elec 36'
            if strcmp(runType,'LC')
                fIndices = find(contains(trialFilenames,'Set053'));
            elseif strcmp(runType,'RC')
                fIndices = find(contains(trialFilenames,'Set054'));
            else
                disp(['No Files Found for' setName]);
            end
        case 'Elec 37'
        case 'Elec 38'
        case 'Elec 39'
        case 'Elec 40'
        otherwise
                disp('TWAT ARE YOU THINKING?!?!?!?!');
    end
    
    if exist('fIndices')
        disp(fIndices);
    else 
        continue;
    end
    %% Plotting Stim Averages
    [sortedStim,I] = sort(stimAmp(fIndices));
    [sortedPulse,I2] = sort(pulseWidth(fIndices));
    [sortedFreq,I3] = sort(stimFrequency(fIndices));

    for c = 1:length(Muscle_ID)
        disp(['Plotting Stim Averages: ' setName ]);

        hF(c) = figure; %('visible', 'off');
        tmpsup = suptitle([Muscle_ID{c} ' - ' setDescrpt ' -  Stim Traces']);
        set(tmpsup, 'Interpreter', 'none') 
        maximize;

        subCount = 1; 
        for fnum2 = 1:iter:length(fIndices)
            hFs(fnum2) = subplot(round(length(fIndices)/iter)+1,1,subCount);
    %         hFs(fnum) = subplot(10,1,subCount);
            %Plot /Channel
            plot(trial(fIndices(I(fnum2))).epochTimeVec,trial(fIndices(I(fnum2))).emgStim{c});hold on;
            plot(trial(fIndices(I(fnum2))).epochTimeVec,trial(fIndices(I(fnum2))).meanTrace{c},'k','LineWidth',1.2);
            vline([0    trial(fIndices(I(fnum2))).stimLength/fs],'r-');
            tmptit = title([trial(fIndices(I(fnum2))).Name ' - ' num2str(stimAmp(fIndices(I(fnum2)))/1000) 'mA, ' ...
                num2str(pulseWidth(fIndices(I(fnum2)))) 'ms, ' num2str(stimFrequency(fIndices(I(fnum2)))) 'Hz, ' num2str(stimDuration(fIndices(I(fnum2)))) 's']); 
            set(hFs(fnum2),'FontSize',7);
            set(tmptit, 'Interpreter', 'none','FontSize',5);
            hFs(fnum2).XLim = [trial(fIndices(I(fnum2))).epochTimeVec(1) trial(fIndices(I(fnum2))).epochTimeVec(end)];

            subCount = subCount + 1;
        end
        linkaxes([hFs],'xy');
        suplabel('Amplitude(uV)','y');
        suplabel('Time(sec)','x');
        set(hF(c),'Position', [1343 49 570 947]);

        disp('Saving to Set Folder')
        savefig(hF(c),[setPath Muscle_ID{c} '_EpochedStimTraces_' setName '_' setDescrpt]);
        saveas(hF(c),[setPath Muscle_ID{c} '_EpochedStimTraces_' setName '_' setDescrpt '.png']);
    end

    close all

    %% Plotting Recruitment Curves
        response = [trial(fIndices).responseRMS]';
        % response = [trial(:).responseP2P]';

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

            subplot(pRow,pCol, pCol*(pI-1)+1)
            for m = 1+mCount*(pI-1):mCount+mCount*(pI-1)
                yyA = smooth(sortedStim/1000,response(I,m),'sgolay');
                plot(sortedStim/1000,yyA,'Color',C(m,:),'LineWidth',3);
        %         scatter(stimAmp/1000, response(I,m)');
                hold on;
            end
                ylabel('(post-Stim) Mean RMS  (uV)')
                xlabel('Stimulation Amplitude (mA)')
                title(['By Amplitude']);

                % See Other Code for PW & Freq plots
                
            if pI == 1
                legend([Muscle_ID(1:mCount)],'Location','northwest', 'Orientation','vertical');
            elseif pI == 2
                legend([Muscle_ID((mCount)+1:length(Muscle_ID))],'Location','northwest', 'Orientation','vertical');
            end
        %     
        end
        tmptit = suptitle([setName ' - ' setDescrpt ' -  Recruitment Curves']);
        set(tmptit, 'Interpreter', 'none');
        % % Pix_SS = get(0,'screensize');
        hF3.Position = [680 49 572 947];

        saveas(hF3,[reportPath '\' strrep(datestr(now),':','_')  '_RecruitmentCurves' '_' setName '_' setDescrpt  '.png']);
        savefig(hF3,[reportPath '\' strrep(datestr(now),':','_') '_RecruitmentCurves_' setName '_' setDescrpt ]);

    
    disp(['End of ' setName ' ' setDescrpt]);       
    
    close all
    
    clear fIndices
end





