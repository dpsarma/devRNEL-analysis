% Muscle_ID = {'Right GM', 'Right RF', 'Right VL', 'Right BF', 'Right SE', 'Right TA', 'Right SO', 'Right LG', 'Left GM', 'Left RF', 'Left VL', 'Left BF', 'Left SE', 'Left VM', 'Left MG', 'Left LG'};
% Muscle_ID = {'Left Gastroc', 'RTA'};
Muscle_ID = {'Right Knee', 'Right RF', 'Right M-MG', 'Right L', 'Right MG', 'Right TA', 'Right SO', 'Right LG'};

Muscle_ID = {'Right VM', 'Right RF', 'Right VL', 'Right BF', 'Right ST', ...
    'Right TA', 'Right SO', 'Right LG','Left VM', 'Left RF', 'Left VL',...
    'Left BF', 'Left ST', 'Left TA', 'Left MG', 'Left LG'};

%%Filtering Parameters
fs = 30000; %Hi-Res is 2K, Raw is 30K (put in logic for switching)
lowCut = 75; %75 lowest frequency to pass
highCut = 750; %7500 highest frequency to pass
% Fsmooth = 5; %lowpass cutoff frequency for smoothing the rectified signal if needed
Norder = 1;
Wp = [lowCut, highCut]/(.5*fs);
[b,a]=butter(Norder, Wp);


%Load data folder
disp('Please Select the Ripple Data Folder');
[emgFilenames, emgPathname] = uigetfile('*ns5','Pick files','MultiSelect', 'on');

%% Load Data.  Need to Make Loop for ALL files.
for i=1%:length(emgFilenames)
    [analogData,timeVec] = read_continuousData([emgPathname emgFilenames{i}], 'raw',[129:144]); %'raw',[257:276]
%     [stimEvts,stchannels] = read_stimEvents([emgPathname emgFilenames{i}],1);
    digEvts = read_digitalEvents([emgPathname, emgFilenames{i}],1);
    stimOn_Idx = digEvts.timeStamp{1}(1:2:end);

    c = 1;
    emg = zeros(size(analogData, 1)/2,size(analogData, 2));
    for n = 1:2:size(analogData, 1)
        emg(c,:) = analogData(n,:) - analogData(n+1,:); 
        c = c + 1;
    end
    
    emgDiffd = diff(analogData); figure;plot(timeVec,emg);hold on; plot(timeVec,emgDiffd,'c'); legend('Manual Bipolar EMG', 'EMG post diff function');
    hline([stimOn_Idx'/1.0e+06])
    line([stimOn_Idx' stimOn_Idx']',repmat([-1000  1000],length(stimOn_Idx),1)','Color',[1,0,0])
%     vline([stimEvts{1}]);
    
    emgF = filtfilt(b,a,emg');
    emgF = emgF';
    
    % Plotting;
%     figure;
    for h = 1:size(emg,1)
        figure;
%         subplot(8,2,h);
        hold on;
%         plot(timeVec,emg(h,:));
        plot(timeVec,emgF(h,:),'g');
        line([stimOn_Idx' stimOn_Idx']',repmat([-1000  1000],length(stimOn_Idx),1)','Color',[1,0,0])
%         vline([stimEvts{1}]);
%         vline([stims]);
%         plot(timeVec,emgSm(i,:),'c');
        title(Muscle_ID(h));
        ylabel('EMG (uV)');
        xlabel('Time (s)');
        legend('raw','filt','stim')%'smooth'
%         ylim([-4000 4000])
%         vline([stimBegin/fs stimEnd/fs]);
    %     vline(0,'r-', 'Stim Pulse');
    %     vline(stimLength/fs,'r-');
    %     axis([tRef(1) tRef(end) -inf inf])
    end
    suptitle(emgFilenames{i});
%     [titString '| ' Muscle_ID(i)]
    
%     idstring = ['R:\data_generated\human\uh3_stim\LSP01\EMGresponses\raw\' erase(emgFilenames{i}, '.ns5') '-'];
% %     saveallopenfigs([tmp,'\EMG-']);
%     figHandles = get(0,'Children');
%     for iter = numel(figHandles):-1:1
%         saveas(figHandles(iter), [idstring num2str(figHandles(iter).Number) '.png']);
%         saveas(figHandles(iter), [idstring num2str(figHandles(iter).Number) '.fig']);
%         close(figHandles(iter));

figure; 
plot(timeVec, emgF); hold on;
% vline([stimEvts{1}]);
 line([stimOn_Idx' stimOn_Idx']',repmat([-1000  1000],length(stimOn_Idx),1)','Color',[1,0,0])
ylabel('EMG (uV)');
xlabel('Time (s)');
legend('filtered EMG','stim');
suptitle(emgFilenames{i});


end   
%     display('Image Save Complete, '); i
