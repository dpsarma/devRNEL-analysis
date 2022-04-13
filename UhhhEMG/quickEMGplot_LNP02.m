Muscle_ID = {'Right VM', 'Right RF', 'Right VL', 'Right BF', 'Right ST', ...
    'Right TA', 'Right SO', 'Right LG','Left VM', 'Left RF', 'Left VL',...
    'Left BF', 'Left ST', 'Left TA', 'Left MG', 'Left LG'};}

%% 1 - Load Data
[analogData,timeVec] = read_continuousData('R:\data_raw\human\uh3_stim\testing\devTest1\datafile0001.nev', 'raw',[129:144]);
    
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

%% 4 -  Plotting

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