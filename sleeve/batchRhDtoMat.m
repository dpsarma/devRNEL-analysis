%% Convert & Process HD-EMG '.rhd' Intan files into .mat files
%% THis is specific for the BCI-WristExtension Experiment

%%Filtering Parameters
fs = 10000; %Hi-Res is 2K, Raw is 30K (put in logic for switching)
lowCut = 20; %75 lowest frequency to pass
highCut = 450; %7500 highest frequency to pass
% Fsmooth = 5; %lowpass cutoff frequency for smoothing the rectified signal if needed
Norder = 1;
Wp = [lowCut, highCut]/(.5*fs);
[b,a]=butter(Norder, Wp);


%Load data folder
disp('Please Select the Intan Data Folder');
dataPath = uigetdir('R:\data_raw\human\emg_sleeve\');
files = dir(fullfile(dataPath, '*.rhd'));

%% Load Data.  Need to Make Loop for ALL files.
tCnt = 0;
for k=1:length(files)
    disp(['Processing: ' files(k).name]); 
    read_IntanRHD_v2(dataPath, files(k).name);
    t_data = t_amplifier;
%     timeDig = t_dig;
    t_adc = t_board_adc;
    raw = amplifier_data(10:end,:);
    filtSig = filtfilt(b,a,raw')';
    rectSig = abs(filtSig);
    if size(rectSig,1) > size(rectSig,2)
        filtrect_data = rectSig';
    else
        filtrect_data = rectSig;
    end
    trigger = board_dig_in_data(1,:);
    force = board_adc_data(1,:);
    tCnt = 1+tCnt;
    
    figure; n = 1; plot(t_data,raw(n,:)); hold on; plot(t_data,filtSig(n,:)); plot(t_data,filtrect_data(n,:)), plot(t_dig,trigger); plot(t_adc,force*300);

    if strncmp(files(k).name, 'Covert',6)
        savName = extractBefore(files(k).name, 22);%12
    else
        savName = extractBefore(files(k).name, 21);%11
    end
    
    %add '_' for Day 1
    save(char([fullfile([erase(dataPath,'Intan') 'Mat'],savName) '_' num2str(tCnt) '.mat']),...
        't_data','t_dig', 't_adc','filtrect_data','trigger','force','-v7.3');
    
    if tCnt == 12
        tCnt = 0;
    elseif tCnt == 11 && strncmp(files(k).name, 'Free',4)
        tCnt = 0;
    end
    
    
    clear t_data t_dig t_adc filtrect_data trigger force savName
    disp(['NEXT FILE - trial Count - ' num2str(tCnt)]);
end