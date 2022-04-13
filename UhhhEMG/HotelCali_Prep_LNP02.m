%% Rehash 2 for Each Subject : LSP02b LSP05

%% Subject Info and Load Data Files
subjectName = 'LNP02'; setName = 'RCs';

        Muscle_ID = {'Right VM', 'Right RF', 'Right VL', 'Right BF', 'Right ST', ...
        'Right TA', 'Right SO', 'Right LG','Left VM', 'Left RF','Left VL',...
        'Left BF', 'Left ST', 'Left TA', 'Left MG', 'Left LG'};
    
%         uiopen('C:\Users\dsarma\Box\HotelCali_paperFigs\resources\LSP05_nev_list.csv',1);

C = linspecer(length(Muscle_ID)/2,'qualitative');

%% Set Paths and Experiment Names
path_datatank = ['R:\data_raw\human\uh3_stim\' subjectName '\data\Trellis\'];
path_stimparams = ['R:\data_raw\human\uh3_stim\' subjectName '\data\Open Loop\stim\'];

reportPath = ['D:\FigRescources\UH3\HotelCali\' subjectName '\']; 
setPath = [reportPath setName '\'];
mkdir(setPath); 

%% Initialize Analysis Parameters

% Filtering Settings:
fs = 30e3;
nyq = 0.5*fs;

%Bandpass
lowCut = 30; %lowest frequency to pass
highCut = 800; %highest frequency to pass
Norder = 2;
Wp = [lowCut, highCut]/(nyq);
[b,a]=butter(Norder, Wp);

%Notch Filter
d = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
               'DesignMethod','butter','SampleRate',fs);
           
%Smoothing Settings:
Fsmooth = 10; %Lowpass cutoff frequency for smoothing the rectified signal 10 vs 50
NorderS =  2;
dt = 1/fs;
[bS,aS] = butter(NorderS,Fsmooth/(nyq));

% % externalStimulator = 'no'; %no = NanoStims
multipolar = 'no';
chanPort = 128; % 128 or 256 based on which Grapevine Port (B or C)
artifactBuffer = 0.005; % post stim artifact settle lag in s, so 0.01 is 10ms
rmsWindow = 0.05; % post stim RMS Window in s, so 0.05 is 50ms

%% Identify files to load
% disp('Please Select the Folder with the Stim Trial Info:');

%% Identify files to load
disp('Please Select the Ripple Data Folder');
[tmpFilenames, emgPathname] = uigetfile([path_datatank '*.nf6'],'Pick files','MultiSelect', 'on');
% % disp(['User selected ', fullfile(emgPathname, emgFilenames)]);

%% Make File List
for f = 1:length(filenames)
    t = find(contains(tmpFilenames',filenames(f)));
    emgFilenames{f} = tmpFilenames(t);
    trialFilenames{f} = [filenames{f} '.mat']; %2014a friendly, no erase function (-4 for just file, -10 for _x# thing)
end
trialPathname = path_stimparams;
%% Load Files and Process
for fnum = 1:length(emgFilenames)
    
    disp(['Loading file ' num2str(fnum) ' of ' num2str(length(emgFilenames))]);
    
    %% Get Basic File Info
    trial(fnum).Name = erase(trialFilenames{fnum}, '.mat'); 
    [ns_status, hFile] = ns_OpenFile([char(emgPathname) char(emgFilenames{fnum})]); 
    [ns_RESULT, nsFileInfo] = ns_GetFileInfo(hFile);
    ns_CloseFile(hFile);
    
    trial(fnum).Date = join(string([nsFileInfo.Time_Month nsFileInfo.Time_Day nsFileInfo.Time_Year]),'-');
    trial(fnum).Time = join(string([nsFileInfo.Time_Hour nsFileInfo.Time_Min nsFileInfo.Time_Sec nsFileInfo.Time_MilliSec]),':');
    
    %% Load Data
    [analogData,timeVec] = read_continuousData([char(emgPathname) char(emgFilenames{fnum})], 'raw' , chanPort+[1:length(Muscle_ID)*2]); %128 vs 256 | 8 vs. 16
    %% 'Data' -> Bipolar EMG
    for i = 1:2:size(analogData,1)
        bipolar = analogData(i:i+1,:);
        emg(round(i/2),:) = diff(flipud(bipolar));
    end
    %% Filtering (Denoising) & Smoothing
    emgF = filtfilt(d,emg')';
    emgF = filtfilt(b,a,emgF')';  
    emgS = filtfilt(bS,aS,emgF')';
    
    trial(fnum).emg = emg;
    trial(fnum).emgF = emgF;
    trial(fnum).timeVec = timeVec;
    disp('Filtering and Epoching....');
    
    %% Load Trial Stim Parameters
    trialInfo = load(fullfile(trialPathname,trialFilenames{fnum}));
        % Identify Stim Events Channel
        switch multipolar
            case 'yes' 
                trial(fnum).spinalElec = trialInfo.Stim_params(1).SpinalElecs{1}(1);
                trial(fnum).anodeElec = trialInfo.Stim_params(1).SpinalElecs{1}(2:end);
                        trial(fnum).stimAmp = trialInfo.Stim_params(1).SpinalElecAmps{1}(1);
            otherwise
                trial(fnum).spinalElec = cell2mat(trialInfo.Stim_params(1).SpinalElecs);
                trial(fnum).stimAmp = cell2mat(trialInfo.Stim_params(1).SpinalElecAmps);
                trial(fnum).anodeElec = [];
        end
        trial(fnum).nevStimCh = trialInfo.Stim_params(1).NSChannels{1,1}(1);   
        % Set Pulse Width, Stim Duration, & Stim Frequency
        trial(fnum).pulseWidth = cell2mat(trialInfo.Stim_params(1).PulseWidth); %provided in ms
        trial(fnum).stimDuration = cell2mat(trialInfo.Stim_params(1).Duration)/1000; %in s, provided in ms
        trial(fnum).stimFrequency = cell2mat(trialInfo.Stim_params(1).Frequency); %this is provided in Hz
        trial(fnum).numStims = trial(fnum).stimDuration*trial(fnum).stimFrequency;
        
        [stimEvts,stchannels] = read_stimEvents(char(fullfile(emgPathname, emgFilenames{fnum})),trial(fnum).nevStimCh);
        stims = floor(cell2mat(stimEvts)*fs);
        
        stimLength = (trial(fnum).pulseWidth/1000)*fs + 2; %Pulse Width in ms, interstim interval is 2 ticks of clock (2/fs)
    trial(fnum).stimLength = stimLength;
    
    if length(stims) == 2*trial(fnum).numStims
        stimBegins = stims(1:2:end);
        stimEnds = stims(2:2:end);
    else
        stimBegins = stims;
    end
    
    trial(fnum).stims = stims;
    
    clearvars -except emgFilenames emgPathname reportPath Muscle_ID fs b a bS aS d nyq...
                fnum trialName trialPathname trialFilenames trial subjectName setName ...
                setDescrpt setPath artifactBuffer rmsWindow ...
                externalStimulator nSampsPre nSampsPost chanPort C multipolar
end

%% Get details

clearvars -except reportPath Muscle_ID ...
                  trial subjectName setName ...
                  setDescrpt setPath C fs
                  

for fnum = 1:length(trial)
    freq = trial(fnum).stimFrequency;
    artifactBuffer = 0.005;
    
    switch freq
        case 1
            dtPre = 50e-3; % msec before stim onset 
            dtPost = 100e-3; % msec after stim onset 
            rmsWindow = 0.05; % post stim RMS Window in s, so 0.05 is 50ms          
        case 2
            dtPre = 50e-3; % msec before stim onset 
            dtPost = 100e-3; % msec after stim onset 
            rmsWindow = 0.05; % post stim RMS Window in s, so 0.05 is 50ms
        case 5
            dtPre = 50e-3; % msec before stim onset 
            dtPost = 100e-3; % msec after stim onset 
            rmsWindow = 0.05; % post stim RMS Window in s, so 0.05 is 50ms
        case 10
            dtPre = 50e-3; % msec before stim onset 
            dtPost = 99e-3; % msec after stim onset 
            rmsWindow = 0.05; % post stim RMS Window in s, so 0.05 is 50ms
        case 20
            dtPre = 50e-3; % msec before stim onset 
            dtPost = 49e-3; % msec after stim onset 
            rmsWindow = 0.030; % post stim RMS Window in s, so 0.05 is 50ms
        case 200
            dtPre = 50e-3; % msec before stim onset 
            dtPost = 4e-3; % msec after stim onset 
            rmsWindow = 0.050; % post stim RMS Window in s, so 0.05 is 50ms
        otherwise 
            dtPre = 50e-3; % msec before stim onset 
            dtPost = 4e-3; % msec after stim onset 
            rmsWindow = 0.005; % post stim RMS Window in s, so 0.05 is 50ms;
    end
    
     %% Build Epoch Time Vector
    nSampsPre = floor(dtPre*fs);
    nSampsPost = floor(dtPost*fs);
    
    trial(fnum).epochTimeVec = linspace(-dtPre, (trial(fnum).stimLength/fs) + dtPost, nSampsPre + nSampsPost + trial(fnum).stimLength);
    stimBegins = trial(fnum).stims;
    stimEnd = nSampsPre + trial(fnum).stimLength + artifactBuffer*fs;  
    %% Epoch data into Windows 
    disp(['Epoching' num2str(fnum)]);
    for i = 1:size(trial(fnum).emg,1)
        asd = trial(fnum).emgF(i,:);
        emgStim{i} = cell2mat(arrayfun(@(x) asd((stimBegins(x)-nSampsPre):(stimBegins(x)+(trial(fnum).stimLength+nSampsPost-1))), 1:length(stimBegins), 'UniformOutput',false)');
    end
    
    
    trial(fnum).emgStim = emgStim;
    trial(fnum).baseline = trial(fnum).emgF(:,1:stimBegins(1));
    
    %% Get RMS acros epochs
    disp('Get Response')
    trial(fnum).epochRMS = cellfun(@rms,emgStim,'UniformOutput', false);
    trial(fnum).epochP2P = cellfun(@peak2peak,emgStim,'UniformOutput', false);
    
    
        for m = 1:length(Muscle_ID)
        % Mean Trace
        trial(fnum).meanTrace{m} = mean(emgStim{1,m});
        trial(fnum).meanbase(m,:) = mean(trial(fnum).baseline(m,:));
        % Response
        if freq > 20
            trial(fnum).responseRMS(m) = mean(rms(trial(fnum).emgF(m,stimBegins(end):stimBegins(end)+(rmsWindow*fs))));
            trial(fnum).responseRMS(m) = mean(peak2peak(trial(fnum).emgF(m,stimBegins(end):stimBegins(end)+(rmsWindow*fs))));
        else
            trial(fnum).responseRMS(m) = mean(trial(fnum).epochRMS{m}(stimEnd:stimEnd+(rmsWindow*fs)));
            trial(fnum).responseP2P(m) = mean(trial(fnum).epochP2P{m}(stimEnd:stimEnd+(rmsWindow*fs)));
        end
    end
end
disp('saving')
save([reportPath '\' setName '_' setDescrpt 'v2.mat'], '-v7.3');
