%H-Reflex Data:

% Muscle_ID = {'Right Knee', 'Right RF', 'Right L-LG', 'Right BF', 'Right MG', 'Right TA', 'Right SO', 'Right LG'};
Muscle_ID = {'Right SO'};
% fullfileName =  'R:\data_raw\human\uh3_stim\testing\UH3_Test\hReflex\datafile0002.nev';
% fullfileName =  'R:\data_raw\human\uh3_stim\testing\devTest2\hiAmp opp-ground\datafile0042.nev';
% fullfileName =  'R:\data_raw\human\uh3_stim\testing\devTest2\no blanking\datafile0020.nev';
% fullfileName =  'R:\data_raw\human\uh3_stim\testing\devTest2\with blanking\datafile0040.nev';
% fullfileName =  'R:\data_raw\human\uh3_stim\testing\devTest1\datafile0002.nev';
fullfileName =  'R:\data_raw\human\uh3_stim\testing\MonicaTest\Trellis\MonicaTest_Ssn007_Set001_Blk001_Trl012.nev';
fullfileName =  'R:\data_raw\human\uh3_stim\testing\devTest3\datafile0019.nev';

% trialMat = 'R:\data_raw\human\uh3_stim\testing\MonicaTest\StimMats\StimRampUp-6mA-1ms-1Hz\MonicaTest_Ssn007_Set001_Blk001_Trl012.mat';
% trialInfo = load(trialMat);

%% 1 - Load Data
[analogData,timeVec] = read_continuousData([fullfileName], 'raw',[159:160]);

% % W/Ripple
[stimEvts,stchannels] = read_stimEvents([fullfileName],1);
stims = stimEvts{1}(1:2:end);

% W/External TTL
% digEvts = read_digitalEvents([fullfileName],1);
% stims = digEvts.timeStamp{1}(1:2:end);

%% 2 - Filtering Parameters
fs = 30000; %Hi-Res is 2K, Raw is 30K (put in logic for switching)
lowCut = 75; %lowest frequency to pass
highCut = 7500; %highest frequency to pass
% Fsmooth = 5; %lowpass cutoff frequency for smoothing the rectified signal if needed
Norder = 1;
Wp = [lowCut, highCut]/(.5*fs);
[b,a]=butter(Norder, Wp);

%% 3- Differencing & Filtering
c = 1;
emg = zeros(size(analogData, 1)/2,size(analogData, 2));
for i = 1:2:size(analogData, 1)
   emg(c,:) = analogData(i,:) - analogData(i+1,:); 
   c = c + 1;
end

emgF = filtfilt(b,a,emg');
emgF = emgF';

g=1;
figure;plot(timeVec,emg(g,:));hold on;plot(timeVec,emgF(g,:),'g');
line([stims' stims']',repmat([-1000  1000],length(stims),1)','Color',[1,0,0])
ylabel('EMG (uV)');xlabel('Time (s)');title(Muscle_ID(g));legend('raw','filt','stim');

%% 3 - Stim Epoching
 % Set Epoch Window for Stim-triggered Averaging

dtPre2 = 20/1000; %mseconds before stim onset to use
dtPost2 = 50/1000; %mseconds after stim onset to use

% Set Pulse Width, Stim Duration, & Stim Frequency
pulseWidth = 1; %this is provided in ms
% % stimDuration = cell2mat(trialInfo.Stim_params.Duration)/1000; %this is provided in ms
% % stimFrequency = cell2mat(trialInfo.Stim_params.Frequency); %this is provided in Hz
numStims = length(stims);

% Using TrialInfo from .mat
%         pulseWidth = cell2mat(trialInfo.Stim_params.PulseWidth); %this is provided in ms
%         stimDuration = cell2mat(trialInfo.Stim_params.Duration)/1000; %this is provided in ms
%         stimFrequency = cell2mat(trialInfo.Stim_params.Frequency); %this is provided in Hz
%         numStims = stimDuration*stimFrequency;

emgTmp = [];

 for c = 1:size(emgF,1)
    display(['Extracting Epochs...' cell2mat(Muscle_ID(c))]);
    [emgTmp, tRef, stimPulse] = epoch_stimEMG(emg(c,:),stims,fs,[dtPre2 dtPost2],'pulses', pulseWidth, numStims);
    emgPulseEpochs{c} = emgTmp;

    %Plot /Channel
    figure;plot(tRef,emgTmp);hold on;
    vline([0 stimPulse/fs],{'r-','r-'},{'','1ms stim pulse'});
    title([Muscle_ID{c} ': Stim Epochs']);
 end

