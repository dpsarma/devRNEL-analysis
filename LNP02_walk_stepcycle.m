%% Plot and Save figures for EMG Data for Closed Loop Trials


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

mLabels = {"Right TFL", "Right RF", "Right TA", "Right SO", "Right LG", "Right VL",...
    "Left TFL", "Left RF", "Left VL", "Right BF", "Left BF", "Left ST", "Left TA",...
    "Right ST", "Left SO", "Left LG"};
chan_remap = [1 2 6 10 14 3 4 5 7 8 9 11 12 13 15 16]; %To match actual Delsys order

%% Process Files
% % figH = figure; maximize;
for f = 1%:length(emgFilenames)
    %% Find Files
% %     close all
    clearvars -except f trialFilenames nevPath insolePath mLabels chan_remap emgFilenames emgPathname trialtype stepcycle triallength
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
    lowCut = 75; %lowest frequency to pass
    highCut = 750; %highest frequency to pass
    Norder = 2;
    Wp = [lowCut, highCut]/(.5*fs);
    [b,a]=butter(Norder, Wp);
    emg = data.emg;
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

    %% Plot EMG w/INSOLEs    
    disp('plotting')
% %     nexttile;
% %     figH = figure; maximize;
% %     tiledlayout(8,2)
% %     i = 0;
    c = [.875 .875 .875];
    
    emg_time = linspace(0,size(data.time,2)/fs,length(data.time)); %% EMG Time based on insole time length?


    stepcycle(f) = length(starts_R) - 1
    triallength(f) = round(emg_time(end))
    
figure; tiledlayout(4,1); maximize;
    for m = [3 4 6 8]%[9 1 11 3 12 4 14 6 16 8]
        z = emg(chan_remap,:);
        mx_cop = max(max(z(chan_remap(m),:)))*1000;
        nexttile;
        plot(emg_time, z(m,:)*1000)
        hold on;
        mx_cop = max(abs(emg(chan_remap(m),:)))*1000;
        if m >=9
            yyaxis right; 
            plot(t_left, cop_L); ylabel('CoP');
            yyaxis left; hold on;
            for g =  2:length(starts_L)
                x = [t_left(starts_L(g-1)) t_left(stops_L(g)) t_left(stops_L(g)) t_left(starts_L(g-1))]; y = [-mx_cop -mx_cop mx_cop mx_cop];
                 h1 = fill(x, y, c,'FaceAlpha',0.3, 'EdgeColor','none');
                if g == length(starts_L)
                    continue;
                else
                     x2 = [t_left(stops_L(g)) t_left(starts_L(g)) t_left(starts_L(g)) t_left(stops_L(g))]; 
                     h2 = fill(x2, y, 'magenta', 'FaceAlpha',0.05,  'EdgeColor','none');
                end
            end
        else
            yyaxis right; 
            plot(t_right, cop_R); ylabel('CoP');
            yyaxis left; hold on;

            for g =  2:length(starts_R)
                x = [t_right(starts_R(g-1)) t_right(stops_R(g)) t_right(stops_R(g)) t_right(starts_R(g-1))]; y = [-mx_cop -mx_cop mx_cop mx_cop];
                h1 = fill(x, y, c,'FaceAlpha',0.3, 'EdgeColor','none');
                if g == length(starts_R)
                    continue;
                else
                    x2 = [t_right(stops_R(g)) t_right(starts_R(g)) t_right(starts_R(g)) t_right(stops_R(g))]; 
                    h2 = fill(x2, y, 'magenta', 'FaceAlpha',0.05,  'EdgeColor','none');
                end
            end
        end
        box off;
    title(mLabels(chan_remap(m)));
    end
    tit = sgtitle({cell2mat(trialtype), cell2mat(filename)},'interpreter', 'none');
end
%% For each step, calculate the P2P amplitude, then simply make a box plot with those. 

idx =  find(stepcycle==0);
stepcycle(idx) = [];
