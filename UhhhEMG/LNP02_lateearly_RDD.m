%%%%   Quick and Dirty Late/Early RDD %
%LNP02 - E4

Muscle_ID = {'Right TFL', 'Right RF', 'Right VL', 'Right BF', 'Right ST', ...
        'Right TA', 'Right SO', 'Right LG','Left TFL', 'Left RF','Left VL',...
        'Left BF', 'Left ST', 'Left TA', 'Left MG', 'Left LG'};
% Need to double check whether the muscles are VM or TFL.
%% Initialize Analysis Parameters

% Filtering Settings:
fs = 2e3; %2e3; %hi-res: 2k; hifreq: 7.5k
nyq = 0.5*fs;

%Bandpass
lowCut = 50; %lowest frequency to pass
highCut = 150; %highest frequency to pass
Norder = 2;
Wp = [lowCut, highCut]/(nyq);
[b,a]=butter(Norder, Wp);

%Notch Filter
d = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
               'DesignMethod','butter','SampleRate',fs);

%% Loading Files

for lset = 5:6

switch lset
    case 1
        
        trialname = 'LNP02_Ssn026_Set001_Blk001_Trl009';
        disp(trialname);
        nevfile = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Trellis\NIP_EMG\allStimfiles\LNP02_Ssn026_Set001_Blk001_Trl009_x0481.nev';
        nf3file = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Trellis\NIP_EMG\LNP02_Ssn026_Set001_Blk001_Trl009_x0019.nf3';
        nf6file = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Trellis\NIP_EMG\LNP02_Ssn026_Set001_Blk001_Trl009_x0019.nf6';
        stimfile = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Open Loop\stim\LNP02_Ssn026_Set001_Blk001_Trl009.mat';
        chanbuff = 128;
        date = '10-04';

    case 2
        trialname = 'LNP02_Ssn026_Set014_Blk001_Trl003';
        disp(trialname);
        nevfile = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Trellis\NIP_EMG\allStimfiles\LNP02_Ssn026_Set014_Blk001_Trl003_x0653.nev';
        nf3file = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Trellis\NIP_EMG\LNP02_Ssn026_Set014_Blk001_Trl003_x0191.nf3';
        nf6file = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Trellis\NIP_EMG\LNP02_Ssn026_Set014_Blk001_Trl003_x0191.nf6';
        stimfile = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Open Loop\stim\LNP02_Ssn026_Set014_Blk001_Trl003.mat';
        chanbuff = 128;
        date = '10-04';

    case 3
        trialname = 'LNP02_Ssn107_Set003_Blk001_Trl011';
        disp(trialname);
        nevfile = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Trellis\NIP_EMG\allStimfiles\LNP02_Ssn107_Set003_Blk001_Trl011_x2066.nev';
        nf3file = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Trellis\NIP_EMG\LNP02_Ssn107_Set003_Blk001_Trl011_x0197.nf3';
        nf6file = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Trellis\NIP_EMG\LNP02_Ssn107_Set003_Blk001_Trl011_x0197.nf6';
        stimfile = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Open Loop\stim\LNP02_Ssn107_Set003_Blk001_Trl011.mat';
        chanbuff = 0;
        date = '11-20';
    
    case 4
        trialname = 'LNP02_Ssn107_Set008_Blk001_Trl002';
        disp(trialname);
        nevfile = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Trellis\NIP_EMG\allStimfiles\LNP02_Ssn107_Set008_Blk001_Trl002_x2147.nev';
        nf3file = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Trellis\NIP_EMG\LNP02_Ssn107_Set008_Blk001_Trl002_x0279.nf3';
        nf6file = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Trellis\NIP_EMG\LNP02_Ssn107_Set008_Blk001_Trl002_x0279.nf6';
        stimfile = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Open Loop\stim\LNP02_Ssn107_Set008_Blk001_Trl002.mat';
        chanbuff = 0;
        date = '11-20';   

    case 5
        trialname = 'LNP02_Ssn077_Set001_Blk001_Trl005';
        disp(trialname);
        nevfile = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Trellis\NIP_EMG\allStimfiles\LNP02_Ssn077_Set001_Blk001_Trl005_x1333.nev';
        nf3file = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Trellis\NIP_EMG\LNP02_Ssn077_Set001_Blk001_Trl005_x0067.nf3';
        nf6file = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Trellis\NIP_EMG\LNP02_Ssn077_Set001_Blk001_Trl005_x0067.nf6';
        stimfile = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Open Loop\stim\LNP02_Ssn077_Set001_Blk001_Trl005.mat';
        chanbuff = 128;
        date = '11-05';

    case 6
        trialname = 'LNP02_Ssn063_Set001_Blk001_Trl015';
        disp(trialname);
        nevfile = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Trellis\NIP_EMG\allStimfiles\LNP02_Ssn063_Set001_Blk001_Trl015_x1138.nev';
        nf3file = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Trellis\NIP_EMG\LNP02_Ssn063_Set001_Blk001_Trl015_x0329.nf3';
        nf6file = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Trellis\NIP_EMG\LNP02_Ssn063_Set001_Blk001_Trl015_x0329.nf6';
        stimfile = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Open Loop\stim\LNP02_Ssn063_Set001_Blk001_Trl015.mat';
        chanbuff = 128;
        date = '11-2';

    case 7
        trialname = 'LNP02_Ssn063_Set002_Blk001_Trl014';
        disp(trialname);
        nevfile = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Trellis\NIP_EMG\allStimfiles\LNP02_Ssn063_Set002_Blk001_Trl014_x1152.nev';
        nf3file = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Trellis\NIP_EMG\LNP02_Ssn063_Set002_Blk001_Trl014_x0343.nf3';
        nf6file = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Trellis\NIP_EMG\LNP02_Ssn063_Set002_Blk001_Trl014_x0343.nf6';
        stimfile = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Open Loop\stim\LNP02_Ssn063_Set002_Blk001_Trl014.mat';
        chanbuff = 128;
        date = '11-2';

    case 8
        trialname = 'LNP02_Ssn026_Set005_Blk001_Trl013';
        disp(trialname);
        nevfile = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Trellis\NIP_EMG\allStimfiles\LNP02_Ssn026_Set005_Blk001_Trl013_x0539.nev';
        nf3file = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Trellis\NIP_EMG\LNP02_Ssn026_Set005_Blk001_Trl013_x0077.nf3';
        nf6file = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Trellis\NIP_EMG\LNP02_Ssn026_Set005_Blk001_Trl013_x0077.nf6';
        stimfile = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Open Loop\stim\LNP02_Ssn026_Set005_Blk001_Trl013.mat';
        chanbuff = 128;
        date = '10-4';

    case 9
        trialname = 'LNP02_Ssn026_Set012_Blk001_Trl007';
        disp(trialname);
        nevfile = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Trellis\NIP_EMG\allStimfiles\LNP02_Ssn026_Set012_Blk001_Trl007_x0635.nev';
        nf3file = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Trellis\NIP_EMG\LNP02_Ssn026_Set012_Blk001_Trl007_x0173.nf3';
        nf6file = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Trellis\NIP_EMG\LNP02_Ssn026_Set012_Blk001_Trl007_x0173.nf6';
        stimfile = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Open Loop\stim\LNP02_Ssn026_Set012_Blk001_Trl007.mat';
        chanbuff = 128;
        date = '10-4';

    case 10
        trialname = 'LNP02_Ssn040_Set003_Blk001_Trl002';
        disp(trialname);
        nevfile = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Trellis\NIP_EMG\allStimfiles\LNP02_Ssn040_Set003_Blk001_Trl002_x0426.nev';
        nf3file = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Trellis\NIP_EMG\LNP02_Ssn040_Set003_Blk001_Trl002_x0030.nf3';
        nf6file = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Trellis\NIP_EMG\LNP02_Ssn040_Set003_Blk001_Trl002_x0030.nf6';
        stimfile = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Open Loop\stim\LNP02_Ssn040_Set003_Blk001_Trl002.mat';
        chanbuff = 128;
        date = '10-14';

    case 11
        trialname = 'LNP02_Ssn040_Set005_Blk001_Trl010';
        disp(trialname);
        nevfile = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Trellis\NIP_EMG\allStimfiles\LNP02_Ssn040_Set005_Blk001_Trl010_x0468.nev';
        nf3file = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Trellis\NIP_EMG\LNP02_Ssn040_Set005_Blk001_Trl010_x0072.nf3';
        nf6file = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Trellis\NIP_EMG\LNP02_Ssn040_Set005_Blk001_Trl010_x0072.nf6';
        stimfile = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Open Loop\stim\LNP02_Ssn040_Set005_Blk001_Trl010.mat';
        chanbuff = 128;
        date = '10-14';

end

%% Filtering (Denoising) & Smoothing
[emg,timeVec] = read_continuousData(nf3file, 'hi-res' , chanbuff+[1:length(Muscle_ID)]);
    emgF = filtfilt(d,emg')';
    emgF = filtfilt(b,a,emgF')'; 

    trialInfo = load(stimfile);
    nevStimCh = trialInfo.Stim_params(1).NSChannels{1,1}(1);
    pulseWidth = cell2mat(trialInfo.Stim_params(1).PulseWidth);
    stimDuration = cell2mat(trialInfo.Stim_params(1).Duration)/1000;
    stimFrequency = cell2mat(trialInfo.Stim_params(1).Frequency);
    stimAmp = trialInfo.Stim_params(1).SpinalElecAmps{1}(1);
    [stimEvts,stchannels] = read_stimEvents(nevfile, nevStimCh);
        stims = floor(cell2mat(stimEvts)*fs);

 figure; 
   for m = 9:16
       ax(lset) = nexttile;
       plot(timeVec,emgF(m,:)); hold on; vline(stimEvts{1, 1});
       box off;
       title(Muscle_ID(m));
       xlabel('time(s)');
       ylabel('EMG(mV)')
%        ylim([-10 10])
   end
   sgtitle([date ' - ' num2str(stimAmp/1000) 'mA - ' num2str(stimFrequency) 'Hz'],'Interpreter','none');
end

% % linkaxes([ax(:)],'y'); 