subjectName = 'LSP02b';

path_datatank = 'R:\data_raw\human\uh3_stim\LSP02b\data\Open Loop\Trellis\';
reportPath = ['R:\data_generated\human\uh3_stim\' subjectName '\emgRecruitment_Summary\'];
% path_stimparams = 'R:\data_generated\human\uh3_stim\genSummary\';

Muscle_ID = {'Right VM', 'Right RF', 'Right VL', 'Right BF', 'Right ST', ...
    'Right TA', 'Right SO', 'Right LG', 'Left VM', 'Left RF', 'Left VL',...
    'Left BF', 'Left ST', 'Left Ham', 'Left SO', 'Left LG'};


chanPort = 128;
% C = linspecer(length(Muscle_ID),'qualitative');

% Filtering Settings:
fs = 30e3;
lowCut = 75; %lowest frequency to pass
highCut = 7500; %highest frequency to pass
Norder = 2;
Wp = [lowCut, highCut]/(.5*fs);
[b,a]=butter(Norder, Wp);

%Smoothing Settings:
Fsmooth = 10; %Lowpass cutoff frequency for smoothing the rectified signal 10 vs 50
NorderS =  2;
dt = 1/fs;
[bS,aS] = butter(NorderS,Fsmooth/(0.5*fs)); 

disp('Please Select the Ripple Data Folder');
[emgFilenames, emgPathname] = uigetfile([path_datatank subjectName '\*.nev'],'Pick files','MultiSelect', 'on');
disp(['User selected ', fullfile(emgPathname, emgFilenames)]);

for f = 1:length(emgFilenames)
    trialFilenames{f} = [emgFilenames{f}(1:end-8)]; %2014a friendly, no erase function (-4 for just file, -10 for _x# thing)
end

for fnum = 2:length(emgFilenames)

%% Load Data
    [analogData,timeVec] = read_continuousData([emgPathname emgFilenames{fnum}], 'raw' , chanPort+[1:length(Muscle_ID)*2]); %128 vs 256 | 8 vs. 16
    %% 'Data' -> Bipolar EMG
    for i = 1:2:size(analogData,1)
        bipolar = analogData(i:i+1,:);
        emg(round(i/2),:) = diff(flipud(bipolar));
    end
    %% Filtering (Denoising) & Smoothing
    emgF = filtfilt(b,a,emg')';  
    emgS = filtfilt(bS,aS,emgF')';
    
    disp('Filtering and Epoching....')
    
    hF(fnum) = figure;
    maximize;
    tmptit = suptitle([join([' EMG of Phantom Movements: ' trialFilenames(fnum)])]);
    set(tmptit,'Interpreter', 'none');
    index = reshape(1:length(Muscle_ID), 2, 8).';
    for m = 1:length(Muscle_ID)
        hS(m) = subplot(8,2,index(m))
            
        plot(timeVec,emgF(m,:));
        xlabel('Time (sec)');
        ylabel('(V)');
        title([Muscle_ID(m)]);
    end
    linkaxes([hS],'xy');
    
    saveas(hF(fnum),[reportPath 'EMGPhantom_' cell2mat(trialFilenames(fnum)) '.png']);
    savefig(hF(fnum),[reportPath 'EMGPhantom_' cell2mat(trialFilenames(fnum))]);
    
     
    clearvars -except emgPathname emgFilenames chanPort Muscle_ID ...
        b a bS aS fnum trialFilenames reportPath
end
    