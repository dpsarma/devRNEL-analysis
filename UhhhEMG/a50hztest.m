Muscle_ID = {'Right VM', 'Right RF', 'Right VL', 'Right BF', 'Right ST', ...
    'Right TA', 'Right SO', 'Right LG','Left VM', 'Left RF', 'Left VL',...
    'Left BF', 'Left ST', 'Left TA', 'Left MG', 'Left LG'};
twin = 0.1; %100 ms or 3 sec.
% Filtering Settings:
fs = 30e3;
lowCut = 75; %lowest frequency to pass
highCut = 750; %highest frequency to pass
Norder = 2;
Wp = [lowCut, highCut]/(.5*fs);
[b,a]=butter(Norder, Wp);    

[analogData,timeVec] = read_continuousData('\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LSP05\data\Trellis\LSP05_Ssn014_Set013_Blk001_Trl003_x0269.nev', 'raw',128+[1:length(Muscle_ID)*2]); %'raw',[257:276]   
stimEvts = read_stimEvents('\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LSP05\data\Trellis\LSP05_Ssn014_Set013_Blk001_Trl003_x0269.nev',29);
stims = floor(cell2mat(stimEvts)*fs);

    c = 1;
    emg = zeros(size(analogData, 1)/2,size(analogData, 2));
    for n = 1:2:size(analogData, 1)
        emg(c,:) = analogData(n,:) - analogData(n+1,:); 
        c = c + 1;
    end

    
% %     emgDiffd = diff(analogData); figure;plot(timeVec,emg);hold on; plot(timeVec,emgDiffd,'c'); legend('Manual Bipolar EMG', 'EMG post diff function');
% % %     hline([stimOn_Idx'/1.0e+06])
% % %     line([stimOn_Idx' stimOn_Idx']',repmat([-1000  1000],length(stimOn_Idx),1)','Color',[1,0,0])
% % %     vline([stimEvts{1}]);
    
    emgF = filtfilt(b,a,emg');
    emgF = emgF';
    t_strt = floor(stims(1)-0.05*fs);
    t_stop = floor(stims(1)+twin*fs); %twin 1.25
    % Plotting;
    figure;% nexttile;
    for h = 1:8%9:size(emg,1)
%         figure; box off
nexttile; box off
        
%         plot(timeVec,emg(h,:));
        plot(timeVec(t_strt:t_stop),emgF(h,t_strt:t_stop));
        hold on;
        vline([timeVec(stims(1:6))]);% timeVec(stims(50))]);
        title(Muscle_ID(h));
        ylabel('EMG (uV)');
        xlabel('Time (s)');
        ylim([-600 600]);
%         vline([stimBegin/fs stimEnd/fs]);
    %     vline(0,'r-', 'Stim Pulse');
    %     vline(stimLength/fs,'r-');
    %     axis([tRef(1) tRef(end) -inf inf])
    end
    sgtitle('50hz');
%     [titString '| ' Muscle_ID(i)]
    
%     idstring = ['R:\data_generated\human\uh3_stim\LSP01\EMGresponses\raw\' erase(emgFilenames{i}, '.ns5') '-'];
% %     saveallopenfigs([tmp,'\EMG-']);
%     figHandles = get(0,'Children');
%     for iter = numel(figHandles):-1:1
%         saveas(figHandles(iter), [idstring num2str(figHandles(iter).Number) '.png']);
%         saveas(figHandles(iter), [idstring num2str(figHandles(iter).Number) '.fig']);
%         close(figHandles(iter));

% % figure; 
% % plot(timeVec(1:2*fs), emgF(1:2*fs)); hold on;
% % % vline([stimEvts{1}]);
% %  line([stimOn_Idx' stimOn_Idx']',repmat([-1000  1000],length(stimOn_Idx),1)','Color',[1,0,0])
% % ylabel('EMG (uV)');
% % xlabel('Time (s)');
% % legend('filtered EMG','stim');
% % sgtitle(emgFilenames{i});

% % linkaxes([fh], 'xy');

