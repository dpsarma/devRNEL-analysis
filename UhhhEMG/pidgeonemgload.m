%% UH3 Data Loading/Vizualization Script

%% 1 - Load Data Files: Raw EMG & Stim Events Loaded Separately
% Still need to add folder indexing and automatic plotting elements.
% Muscle_ID = {'Right GM', 'Right RF', 'Right VL', 'Right BF', 'Right SE', 'Right TA', 'Right SO', 'Right LG', 'Left GM', 'Left RF', 'Left VL', 'Left BF', 'Left SE', 'Left VM', 'Left MG', 'Left LG'};
Muscle_ID = {'Right Knee', 'Right RF', 'Right L-LG', 'Right BF', 'Right MG', 'Right TA', 'Right SO', 'Right LG'};
% Muscle_ID = {'Left VM', 'Left TA'};
% tmp = 'R:\data_raw\test\UH3_Test\6_29_emgTest\stimParams\ACNtest_Ssn011_Set001_Blk001_Trl010.mat';
% trialInfo = load(tmp);
% nevStimChannel = trialInfo.Stim_params.NSChannels{1,1}(1);

% R:\data_raw\human\uh3_stim\testing\EMG
% R:\data_raw\test\UH3_Test\openLoop
% R:\data_raw\human\uh3_stim\testing2\TestTrial
% R:\data_raw\test\UH3_Test\5_1_18_emgTest\trellis %09
% R:\data_raw\human\uh3_stim\LSP01\data\Open Loop\Trellis\LSP01_Ssn003_Set032_Blk001_Trl034
[analogData,timeVec] = read_continuousData('R:\data_raw\test\UH3_Test\6_29_emgTest\trellis\ACNtest_Ssn011_Set001_Blk001_Trl010.nev', 'raw',[129:132]); %'hi-res',[257:276]
% [stimEvts,stchannels] = read_stimEvents('R:\data_raw\test\UH3_Test\6_29_emgTest\trellis\ACNtest_Ssn011_Set001_Blk001_Trl010.nev',nevStimChannel);%nevStimChannel); %vs 34
% tmp = 'R:\data_raw\human\uh3_stim\LSP01\data\Open Loop\stim\LSP01_Ssn003_Set017_Blk001_Trl012.mat';
%% 2 - Filtering Parameters
fs = 30000; %Hi-Res is 2K, Raw is 30K (put in logic for switching)
lowCut = 75; %lowest frequency to pass
highCut = 7500; %highest frequency to pass
% Fsmooth = 5; %lowpass cutoff frequency for smoothing the rectified signal if needed
Norder = 1;
Wp = [lowCut, highCut]/(.5*fs);
[b,a]=butter(Norder, Wp);


%% 3 - Stim Epoching Parameters
dt = 1/fs;
dtPre = 0.95; %seconds before stim onset to use
dtPost = 1; %seconds after stim onset to use
nSampsPre = floor(dtPre*fs);
nSampsPost = floor(dtPost*fs);

stims = stimEvts{1,1};
stimEnd = floor(stims(end)*fs);
stimBegin = floor(stims(1)*fs);
stimLength = stimEnd-stimBegin+1;
tRef = linspace(-dtPre, (stimLength/fs)+dtPost, nSampsPre+nSampsPost+stimLength); %epoch "time"


%% 4- Differencing
c = 1;
emg = zeros(size(analogData, 1)/2,size(analogData, 2));
for i = 1:2:size(analogData, 1)
   emg(c,:) = analogData(i,:) - analogData(i+1,:); 
   c = c + 1;
end
   
%% 5 -  Filtering & Stim Epoching
emgF = filtfilt(b,a,emg');
emgF = emgF';
% [B,A] = butter(Norder, Fsmooth/(.5*fs));
% emgSm = filtfilt(B,A,abs(emgF));

emgStim=emgF(:,stimBegin-nSampsPre:stimEnd+nSampsPost); %emgF
% emgRaw = emg(:,stimBegin-(preStim_window*fs):stimEnd+(emgtime_window*fs));


%% 6 - Plotting

% Stim Plot (16 channels, epoched)
% figure;
% % Muscle_ID = {'Right TFL', 'Right RF', 'Right VL', 'Right BF', 'Right SE', 'Right TA', 'Right SO', 'Right LG', 'Left TFL', 'Left RF', 'Left VL', 'Left BF', 'Left SE', 'Left TA', 'Left SO', 'Left LG'};
% for i = 1:size(emgStim,1)
%     h(i) = subplot(2,size(emgStim,1)/2,i);
%     plot(tRef,emgStim(i,:),'g');
%     title(Muscle_ID(i));
%     vline(0,'r-', 'Stim Train Begin');
%     vline(stimLength/fs,'r-', 'Stim Train End');
%     axis([tRef(1) tRef(end) -inf inf])
%     ylabel('EMG (uV)');
%     xlabel('Time (s)');
% end
% suptitle('Stim Train + EMG Response (Filtered)')
% maximize;

%Filter Comparison Plot (16 Channels full timeseries)
% stims =  stimEvts{1,1};
% emgF = emgF';
figure;
% Muscle_ID = {'Right TFL', 'Right RF', 'Right VL', 'Right BF', 'Right SE', 'Right TA', 'Right SO', 'Right LG', 'Left TFL', 'Left RF', 'Left VL', 'Left BF', 'Left SE', 'Left TA', 'Left SO', 'Left LG'};
for i = 1:size(emg,1)
    h(i) = subplot(2,size(emg,1)/2,i);
    hold on;
    plot(timeVec,emg(i,:));
    plot(timeVec,emgF(i,:),'g');
%     vline([stims]);
%     plot(timeVec,emgSm(i,:),'c');
    title(Muscle_ID(i));
    ylabel('EMG (uV)');
    xlabel('Time (s)');
    legend('raw','filt')%'smooth'
    vline([stims]);
%     ylim([-4000 4000])
%     vline(0,'r-', 'Stim Pulse');
%     vline(stimLength/fs,'r-');
%     axis([tRef(1) tRef(end) -inf inf])
end
% for i = 1:size(emg,1)
%     subplot(4,5,10+i)
%     hold on;
% %     plot(timeVec,emg(i,:));
%     plot(timeVec,emgF(i,:),'r');
% %     plot(timeVec,emgSm(i,:),'r');
%     title(Muscle_ID(i));
%     ylabel('EMG (uV)');
%     xlabel('Time (s)');
%     legend('filt')
% %     vline(0,'r-', 'Stim Pulse');
% %     vline(stimLength/fs,'r-');
% %     axis([tRef(1) tRef(end) -inf inf])
% end
suptitle('Raw EMG vs Filtered(75:7500Hz)')
maximize;


%%Stim Averaging

% stims = stimEvts{1,1};
% stims.freq = trialInfo.Stim_params.Frequency;
% stims.plswdth = trialInfo.Stim_params.PulseWidth; %s or ms???
dt2 = 1/fs;
dtPre2 = 0.020; %20 mseconds before stim onset to use
dtPost2 = 0.050; %50 mseconds after stim onset to use
nSampsPre2 = floor(dtPre2*fs);
nSampsPost2 = floor(dtPost2*fs);

 % Set Pulse Width, Stim Duration, & Stim Frequency
pulseWidth = cell2mat(trialInfo.Stim_params.PulseWidth); %this is provided in ms
stimDuration = cell2mat(trialInfo.Stim_params.Duration)/1000; %this is provided in ms
stimFrequency = cell2mat(trialInfo.Stim_params.Frequency); %this is provided in Hz
numStims = size(stims,2); %stimDuration*stimFrequency;
     
stimBegins2 = floor(stims*fs);
stimLength2 = (2*pulseWidth/1000)*fs + 2; %Pulse Width in ms, interstim interval is 2 ticks of clock (2/fs)
tRef2 = linspace(-dtPre2, (stimLength2/fs)+dtPost2, nSampsPre2+nSampsPost2+stimLength2);
emgTmp = [];

for c = 1:size(emgF,1)
    display('Extracting Epochs...');
    for n = 1:numStims
        startSamp = (stimBegins2(n) - nSampsPre2);
        stopSamp = (startSamp + nSampsPre2 + stimLength2 + nSampsPost2);
        emgTmp(n,:)=emgF(c,startSamp:stopSamp-1);
    end
    emgAvg{c} = emgTmp;
    figure;plot(tRef2,emgTmp);hold on;
    vline([0 stimLength2/fs]);
    title([Muscle_ID{c} ': Stim Epochs']);
%     emgStimRMS(c,:,:) = rms(emgTmp);
end

tmpMat = [];
for c=1:size(emgAvg,2)
    
    tmpMat = emgAvg{1,c};
    emgStimRMS(c,:) = rms(tmpMat);
end

% Stim Plot (16 channels, epoched)
figure;
% Muscle_ID = {'Right TFL', 'Right RF', 'Right VL', 'Right BF', 'Right SE', 'Right TA', 'Right SO', 'Right LG', 'Left TFL', 'Left RF', 'Left VL', 'Left BF', 'Left SE', 'Left TA', 'Left SO', 'Left LG'};
for i = 1:size(emgStimRMS,1)
    h(i) = subplot(2,size(emgStimRMS,1)/2,i);
    plot(tRef2,emgStimRMS(i,:));
    title(Muscle_ID(i));
    vline(0,'r-', 'Stim Pulse');
    vline(stimLength2/fs,'r-');
    axis([tRef2(1) tRef2(end) -inf inf])
    ylabel('EMG (uV)');
    xlabel('Time (s)');
end
suptitle('StimTriggered RMS EMG Responses')
maximize;

%% Mean vs. RMS
tmpMat = [];
for c=1:size(emgAvg,2)
    
    tmpMat = emgAvg{1,c};
    emgStimMean(c,:) = mean(tmpMat);
end
% Stim Plot (16 channels, epoched)
figure;
% Muscle_ID = {'Right TFL', 'Right RF', 'Right VL', 'Right BF', 'Right SE', 'Right TA', 'Right SO', 'Right LG', 'Left TFL', 'Left RF', 'Left VL', 'Left BF', 'Left SE', 'Left TA', 'Left SO', 'Left LG'};
for i = 1:size(emgStim,1)
    h(i) = subplot(2,size(emgStim,1)/2,i);
    plot(tRef2,emgStimRMS(i,:));hold on;
    plot(tRef2,emgStimMean(i,:),'g');
    title(Muscle_ID(i));
    vline(0,'r-', 'Stim Pulse');
    vline(stimLength2/fs,'r-');
    axis([tRef2(1) tRef2(end) -inf inf])
    ylabel('EMG (uV)');
    xlabel('Time (s)');
    legend('RMS', 'Mean');
end
suptitle('StimTriggered RMS (blue) vs. Avg (green) EMG Responses')
maximize;


% %Plot Raw Monopolar Plot
% figure;
% for i = 1:size(analogData,1)
%     subplot(16,2,i)
%     hold on;
%     plot(timeVec,analogData(i,:));
%     ylabel('EMG (uV)');
%     xlabel('Time (s)');
% 
% end
% suptitle('Analog')
% maximize;



