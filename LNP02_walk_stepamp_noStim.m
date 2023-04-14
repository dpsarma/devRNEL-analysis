%% Plot and Save figures for EMG Data for Closed Loop Trials No Stim


%% Get Files
% Identify files to load
disp('Please Select the Ripple Data Folder');
[tmpFilenames, emgPathname] = uigetfile(['C:/data/LL_UH3/ClosedLoop/' '*.pkl'],'Pick files','MultiSelect', 'on');
nevPath = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Closed Loop\';
insolePath = '\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Closed Loop\Insole_data\';

% Make File List
for f = 1:length(tmpFilenames)
    emgFilenames{f} = tmpFilenames(f);
    trialFilenames{f} = erase(emgFilenames{f}, '.pkl'); 
end
out=regexp(emgPathname,'\','split'); trialtype = out(end-1);

mLabels = {'Right TFL', 'Right RF', 'Right TA', 'Right SO', 'Right LG', 'Right VL',...
    'Left TFL', 'Left RF', 'Left VL', 'Right BF', 'Left BF', 'Left ST', 'Left TA',...
    'Right ST', 'Left SO', 'Left LG'};
chan_remap = [1 2 6 10 14 3 4 5 7 8 9 11 12 13 15 16]; %To match actual Delsys order

p = cell(16,1); r = cell(16,1); pz = cell(16,1); rz = cell(16,1);
mlist = [8 2 10 12 14]; %[10 2 12 4 14 6 16 8];
%% Process Files
% % figH = figure; maximize;
% % close all;
for f = 1:length(emgFilenames)
    %% Find Files
% %     close all
    filename = trialFilenames{f};
    fstruct = dir([nevPath cell2mat(filename) '*.nev']); %stim files
    istruct  = dir([insolePath cell2mat(filename) '*.json']); %insole files

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
        num_steps = size(centers_R,2);

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
    emg = (highpass(data.emg', 65 , fs)')*1000*1000;
%     emg = data.emg;
    z = emg(chan_remap,:);
    z = zscore(z')';
%     emg = filtfilt(b,a,data.emg')';
    tend = length(data.time)/fs;
% %     emg_time = linspace(0,tend,length(data.time));


%% Assign smoothing type

emg_time = linspace(0,size(data.time,2)/fs,length(data.time)); %% EMG Time based on insole time length?  
if contains(filename, '99')
    e = z(8,:); 
    t = emg_time;
    e = abs(e - mean(e));
    e_env = envelope(e, 150,'rms');
    idx = findchangepts(e_env,'Statistic','mean','MinThreshold', 10);%,'MaxNumChanges',num_steps*2);
% % elseif contains(filename, '93')
% %     e = z(8,:); 
% %     t = emg_time;
% %     e = abs(e - mean(e));
% %     e_env = envelope(e, 150,'rms');
% %     idx = findchangepts(e_env,'Statistic','mean','MinThreshold', 10);%,'MaxNumChanges',num_steps*2);
else 
        e = lowpass(emg(8,:)', 1 , fs)';
% %     e = z(8,:); 
    t = emg_time;
    e = abs(e - mean(e));
    e_env = envelope(e, 150,'rms');
    idx = findchangepts(e_env,'Statistic','mean','MinThreshold', 200);%,'MaxNumChanges',num_steps*2);
end



%% Assign Vars
    switch filename{1}
            case 'LNP02_CL_Ssn068_Set001_Blk001_Trl001'
% %                 display(filename); continue;
            case 'LNP02_CL_Ssn068_Set001_Blk001_Trl002'
                 display(filename);
            case 'LNP02_CL_Ssn068_Set001_Blk001_Trl003'
                display(filename);
            case 'LNP02_CL_Ssn068_Set001_Blk001_Trl004'
                display(filename);                                
            case 'LNP02_CL_Ssn068_Set001_Blk001_Trl005'
                display(filename);
            case 'LNP02_CL_Ssn068_Set001_Blk001_Trl006'
                display(filename);
            case 'LNP02_CL_Ssn068_Set001_Blk001_Trl007'
                display(filename);
            case 'LNP02_CL_Ssn068_Set001_Blk001_Trl008'
                display(filename);
            case 'LNP02_CL_Ssn068_Set001_Blk001_Trl009'
                display(filename);
            case 'LNP02_CL_Ssn075_Set001_Blk001_Trl001'
                display(filename);
            case 'LNP02_CL_Ssn075_Set001_Blk001_Trl002'
                display(filename);
            case 'LNP02_CL_Ssn075_Set001_Blk001_Trl003'
                display(filename);
            case 'LNP02_CL_Ssn075_Set001_Blk001_Trl004'
                display(filename);
            case 'LNP02_CL_Ssn075_Set001_Blk001_Trl005'
                display(filename);
            case 'LNP02_CL_Ssn075_Set001_Blk001_Trl006'
                display(filename);
            case 'LNP02_CL_Ssn092_Set001_Blk001_Trl001'
                display(filename);
            case 'LNP02_CL_Ssn092_Set001_Blk001_Trl002'
                display(filename);
            case 'LNP02_CL_Ssn092_Set001_Blk001_Trl003'
                display(filename);
            case 'LNP02_CL_Ssn092_Set001_Blk001_Trl004'
                display(filename);
            case 'LNP02_CL_Ssn092_Set001_Blk001_Trl005'
                display(filename);
            case 'LNP02_CL_Ssn092_Set001_Blk001_Trl006'
                display(filename);
	        case 'LNP02_CL_Ssn099_Set001_Blk001_Trl001'
                display(filename);
                idx_sta = idx([11 20 26 34 39 45 50]);
                idx_sto = idx([20 26 34 39 45 50 57]);
            case 'LNP02_CL_Ssn099_Set001_Blk001_Trl002'
                display(filename);
                idx_sta = idx([6 11 14]);
                idx_sto = idx([11 14 16]);
            case 'LNP02_CL_Ssn099_Set001_Blk001_Trl003'
                display(filename);
                idx_sta = idx([7 11 16 21 23 28]);
                idx_sto = idx([11 16 21 23 28 30]);
            case 'LNP02_CL_Ssn099_Set001_Blk001_Trl004'
                display(filename);
                idx_sta = idx([4 8 11 13]);
                idx_sto = idx([8 11 13 15]);
            case 'LNP02_CL_Ssn099_Set001_Blk001_Trl005'
                display(filename);
                idx_sta = idx([5 10 15 21 26 29]);
                idx_sto = idx([10 15 21 26 29 32]);
            case 'LNP02_CL_Ssn099_Set001_Blk001_Trl006'
                display(filename);
                idx_sta = idx([1 7 14 20 24]);
                idx_sto = idx([7 14 20 24 26]);
    end

% %     figure;
% %     plot(t,emg(8,:)); hold on;
% %     plot(t,e);
% %     hold on;
% %     plot(t,e_env, 'LineWidth',2);
% %     x = string(1:length(idx));
% %     vline(t(idx), 'g-', cellstr(x));
% %     title('Extract Epoch Times from Right LG');
% %     xlabel('time');
% %     ylabel('uV');

% %     idx_vec = [3 5 7 9 11];

        
    
% %     figure; tiledlayout(4,2); maximize;
% % 
% %     for m = mlist %[3 4 6 8]%[9 1 11 3 12 4 14 6 16 8]
% % %         mx_cop = max(max(z((m),:)));
% %         nexttile;
% %         plot(t, z(m,:))
% %         hold on;
% %         plot(t,envelope(z(m,:),100,'rms'),'LineWidth',2);
% %         hold on;
% %         vline(t(idx_sta), 'r-'); hold on;
% % % %         vline(t(idx_sto), 'm-');
% %         title(mLabels(chan_remap(m)));
% %     end


    epoch = []; 
    %% Extract Phase from LG
% %     figure; tiledlayout(4,2); maximize;
    for m = mlist %[2 4 6 8]
        env = envelope(abs(z(m,:)),100,'rms');
% %         nexttile;
% %         % %  
        ppz = []; rrz = []; pp = []; rr = [];
        for i = 1:length(idx_sta)
            imod1 = round((idx_sto(i)-idx_sta(i))*0.2);
            imod2 = round((idx_sto(i)-idx_sta(i))*0.3);%0 to 130%
% %             epoch = env(idx_sta(i)-imod1:idx_sto(i)+imod2); % +imod; abs(z(m,idx(i):idx(i+2)));%
% %             e_tim = linspace(-20,130,length(epoch));
            epoch = env(idx_sta(i):idx_sto(i)); % +imod; abs(z(m,idx(i):idx(i+2)));%
            e_tim = linspace(0,100,length(epoch));
            if length(e_tim) < 2200
               disp('step is too short')
                continue;
            end
            EP(f).m(m).i(i).epoch = epoch;
            EP(f).m(m).i(i).etim = e_tim;
% %             plot(e_tim, epoch); box off; xlim([-20 130]);
% %             hold on;
            ppz(end+1) = peak2peak(z(m,idx_sta(i):idx_sto(i)));
            rrz(end+1) = rms(z(m,idx_sta(i):idx_sto(i)));
            pp(end+1) = peak2peak(emg(m,idx_sta(i):idx_sto(i)));
            rr(end+1) = rms(emg(m,idx_sta(i):idx_sto(i)));
        end
        pz{m} = [pz{m} ppz];
        rz{m} = [rz{m} rrz];
        p{m} = [p{m} pp];
        r{m} = [r{m} rr];
% %         title(mLabels(chan_remap(m)));
% %         ylabel('EMG (z-scored)')
    end
%     sgtitle({cell2mat(trialtype), cell2mat(filename), 'Right Heel Strike:Heel Strike'},'interpreter', 'none');
clearvars -except f trialFilenames nevPath insolePath mLabels chan_remap...
        emgFilenames emgPathname trialtype stepcycle triallength acti p r pz rz mlist EP
end

    %% Bar Plot of Peaks
    allPeak  = cell2mat(rz);
    figure; boxplot(allPeak'); xticklabels(mLabels(chan_remap(sort(mlist))));

    %% All Epochs
    figure; tiledlayout(length(mlist),1); maximize;
    col_map = [0 20 89; 82 130 192; 245 77 0; 250 138 115]/255;

    for m = mlist
        nexttile;
        m
        if m == 10 || m == 2
            i = 1
        elseif m == 12 || m == 4
            i = 3
        elseif m == 14 || m == 6
            i = 4
        elseif m == 8
            i = 2
        else
            disp('muscle not found');
            continue;
        end

        for f = 1:length(EP)
            for count = 1:size(EP(f).m(m).i,2)
            plot(EP(f).m(m).i(count).etim, EP(f).m(m).i(count).epoch, 'Color', col_map(i,:));  
            hold on;
            end
        end
        xlim([0 100]); xticklabels([]); box off;
        title(mLabels(chan_remap(m)));
% %         ylabel({'Activation', '(z-scored EMG)'});
    end
% %     xticklabels([0 20 40 60 80 100]);
    xlabel('% of Gait Cycle')
    sgtitle({cell2mat(trialtype), '%Gait Phase', 'Right Heel Strike:Heel Strike'},'interpreter', 'none');