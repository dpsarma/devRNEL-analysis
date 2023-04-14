% % %% Plot and Save figures for EMG Data for Closed Loop Trials LG side
% % 
% % 
% % %% Get Files
% % % Identify files to load
% % disp('Please Select the Ripple Data Folder');
% % [tmpFilenames, emgPathname] = uigetfile(['C:/data/LL_UH3/ClosedLoop/' '*.pkl'],'Pick files','MultiSelect', 'on');
% % nevPath = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Closed Loop\';
% % insolePath = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Closed Loop\Insole_data\';
% % 
% % % Make File List
% % for f = 1:length(tmpFilenames)
% %     emgFilenames{f} = tmpFilenames(f);
% %     trialFilenames{f} = erase(emgFilenames{f}, '.pkl'); 
% % end
% % out=regexp(emgPathname,'\','split'); trialtype = out(end-1);
% % 
% % mLabels = {'Right TFL', 'Right RF', 'Right TA', 'Right SO', 'Right LG', 'Right VL',...
% %     'Left TFL', 'Left RF', 'Left VL', 'Right BF', 'Left BF', 'Left ST', 'Left TA',...
% %     'Right ST', 'Left SO', 'Left LG'};
% % chan_remap = [1 2 6 10 14 3 4 5 7 8 9 11 12 13 15 16]; %To match actual Delsys order
% % 
% % %% Process Files
% % % % figH = figure; maximize;
% % p = []; r = [];
for f = 1%:length(emgFilenames)
    %% Find Files
% %     close all
    clearvars -except f trialFilenames nevPath insolePath mLabels chan_remap...
        emgFilenames emgPathname trialtype stepcycle triallength acti p r
    filename = trialFilenames{f};
    fstruct = dir([nevPath cell2mat(filename) '*.nev']); %stim files
    istruct  = dir([insolePath cell2mat(filename) '*.json']); %insole files
 
    %% Get EMG Data
    disp('loading EMG');
    data = py.pickle.load( py.open([emgPathname cell2mat(emgFilenames{f})], 'rb'));
    py.scipy.io.savemat('C:/data/tmp.mat', mdict=struct('data', data));
    load('C:\data\tmp.mat');
    data.emg = double(data.emg);

    % Filtering Settings:
    fs = 2000;
% %     lowCut = 75; %lowest frequency to pass
% %     highCut = 750; %highest frequency to pass
% %     Norder = 2;
% %     Wp = [lowCut, highCut]/(.5*fs);
% %     [b,a]=butter(Norder, Wp);
    emg = highpass(data.emg', 65, fs)';
    z = emg(chan_remap,:)*1000*1000;
    z = zscore(z')';
%     emg = filtfilt(b,a,data.emg')';
    tend = length(data.time)/fs;
% %     emg_time = linspace(0,tend,length(data.time));
    
    
    %% Get Insole Data
    disp('extract insole data')
    tmpf = py.open([istruct.folder '\' istruct.name]); % Opening JSON file, needs 'rb'?
    insole_data = py.json.load(tmpf);
    py.scipy.io.savemat('C:/data/tmp_insole.mat', mdict=struct('insole_data', insole_data));
    load('C:\data\tmp_insole.mat');
    for i = 1:length(insole_data)
        t_left(i) = double(insole_data{i}.time_L)/1000;
        left(:,i) = insole_data{i}.pressure_L;
        cop_L(i) = insole_data{i}.cop_L;
        t_right(i) = double(insole_data{i}.time_R)/1000;
        right(:,i) = insole_data{i}.pressure_R;
        cop_R(i) = insole_data{i}.cop_R;
    end
    % Normalize to zero and convert to sec
    t_left = (t_left - t_left(1));
    t_right = (t_right - t_right(1));
    % Find phases Left
        mask = logical(cop_L(:).');    %(:).' to force row vector
        stops_L = strfind([false, mask], [0 1]);
        starts_L = strfind([mask, false], [1 0]);
        centers_L = mean([starts_L;stops_L]);
    %Find phases Right
        mask = logical(cop_R(:).');    %(:).' to force row vector
        stops_R = strfind([false, mask], [0 1]);
        starts_R = strfind([mask, false], [1 0]);
        centers_R = mean([starts_R;stops_R]);
        
    %% Get Stims
    if isempty(fstruct)
        disp([filename ' does not have .nev']);
        stimStatus = 2;
    else
        disp([filename ' has nev']);
            nevFilename = [fstruct(1).folder '\' fstruct(1).name];
            trialInfo = load([nevPath cell2mat(filename) '.mat']);

        if isfield(trialInfo.Stim_params, 'NSChannels')
            nevStimCh = trialInfo.Stim_params(1).NSChannels{1,1}(1);
            pulseWidth = trialInfo.Stim_params.PulseWidth{1}(1);
            stimLength = pulseWidth*fs + 2;
            stimFreq = trialInfo.Stim_params.Frequency{1}(1);
            %get stim events
            [stimEvts] = read_stimEvents(nevFilename,nevStimCh);
            stims = floor(cell2mat(stimEvts)*fs);
            stimStatus = 1;
        else
            disp('No Stim Detected');  
            stimStatus = 0;
        end
    end  

     emg_time = linspace(0,size(data.time,2)/fs,length(data.time)); %% EMG Time based on insole time length?

% %     for m = 1:size(emg,1)
% %         for s = 1:length(starts_L)-1
% %             acti(f,m,s) = peak2peak(emg(chan_remap(m), find(emg_time > t_left(starts_L(s)), 1):find(emg_time > t_left(starts_L(s+1)),1)));
% % %             gg{f,m,s} = ['No Stim - ' mLabels{chan_remap(m)}];
% %         end
% %     end

    e = z(16,:); t = emg_time;
    e = abs(e - mean(e));
    idx = findchangepts(e,'Statistic','mean','MaxNumChanges',20);
    figure;
    plot(t,e);
    hold on;
    vline(t(idx), 'r-');
    title({'Extract Epoch Times from Left LG', filename{1}},'interpreter', 'none');
    xlabel('time');
    ylabel('uV');
    
    figure; tiledlayout(4,2); maximize;

    for m = [10 2 12 4 14 6 16 8]  %[3 4 6 8]%[9 1 11 3 12 4 14 6 16 8]
%         mx_cop = max(max(z((m),:)));
        nexttile;
        plot(t, z(m,:))
        hold on;
        plot(t,envelope(z(m,:),100,'rms'),'LineWidth',2);
        hold on;
        vline(t(idx), 'r-');
        title(mLabels(chan_remap(m)));
    end
    sgtitle({filename{1}, 'Gait Trace'},'interpreter', 'none')
    epoch = [];
    %% Extract Phase from LG
    figure; tiledlayout(4,2); maximize;
    for m = [10 2 12 4 14 6 16 8] %[2 4 6 8]
        env = envelope(abs(z(m,:)),100,'rms');
        nexttile;
        
        idx_vec = [3 5 7 9 11];

        switch filename{1}
            case 'LNP02_CL_Ssn073_Set001_Blk001_Trl001'
                  display(filename);
            case 'LNP02_CL_Ssn073_Set001_Blk001_Trl002'
                  display(filename);
            case 'LNP02_CL_Ssn073_Set001_Blk001_Trl003'
                  display(filename);
            case 'LNP02_CL_Ssn073_Set001_Blk001_Trl004'
                  display(filename);
            case 'LNP02_CL_Ssn073_Set001_Blk001_Trl005'
                  display(filename);
            case 'LNP02_CL_Ssn073_Set001_Blk001_Trl006'
                  display(filename);
            case 'LNP02_CL_Ssn093_Set001_Blk001_Trl001'
                  display(filename);
            case 'LNP02_CL_Ssn093_Set001_Blk001_Trl002'
                  display(filename);
            case 'LNP02_CL_Ssn093_Set001_Blk001_Trl003'
                  display(filename);
            case 'LNP02_CL_Ssn093_Set001_Blk001_Trl004'
                  display(filename);
            case 'LNP02_CL_Ssn093_Set001_Blk001_Trl005'
                  display(filename);
            case 'LNP02_CL_Ssn093_Set001_Blk001_Trl006'
                  display(filename);
            case 'LNP02_CL_Ssn093_Set001_Blk001_Trl007'
                  display(filename);
            case 'LNP02_CL_Ssn097_Set001_Blk001_Trl001'
                  display(filename);
            case 'LNP02_CL_Ssn097_Set001_Blk001_Trl002'
                  display(filename);
            case 'LNP02_CL_Ssn097_Set001_Blk001_Trl003'
                  display(filename);
            case 'LNP02_CL_Ssn097_Set001_Blk001_Trl004'
                    display(filename);
            case 'LNP02_CL_Ssn097_Set001_Blk001_Trl005'
                    display(filename);
            case 'LNP02_CL_Ssn097_Set001_Blk001_Trl006'
                    display(filename);
        end

        for i = idx_vec
            imod1 = round((idx(i+2)-idx(i))*0.2);
            imod2 = round((idx(i+2)-idx(i))*0.3);%0 to 130%
            epoch = env(idx(i)-imod1:idx(i+2)+imod2); % +imod; abs(z(m,idx(i):idx(i+2)));%
            e_tim = linspace(-20,130,length(epoch));
            plot(e_tim, epoch); box off; xlim([-20 130]);
            hold on;
            p(m,end+1) = peak2peak(z(m,idx(i):idx(i+2)));
            r(m,end+1) = rms(z(m,idx(i):idx(i+2)));
        end
%         peak(f,m) = max(p(m,:));
%         peakrms(f,m) = max(r(m,:));
%         mean_env{f,m} = mean(epp);
%         hold on;
%         plot(e_tim, [mean_env(f,m)])
        title(mLabels(chan_remap(m)));
        ylabel('EMG (z-scored)')
    end
    sgtitle({cell2mat(trialtype), cell2mat(filename), 'Right Heel Strike:Heel Strike'},'interpreter', 'none')
    


end

    %% Bar Plot of Peaks
% %     pk = mean(peak)
    figure; nexttile; boxplot(r'); xticklabels(mLabels(chan_remap));
    nexttile; boxplot(p'); xticklabels(mLabels(chan_remap));