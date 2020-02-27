%% Modular Bionics Test Analysis
%created by Dev Sarma, 5/1/2018

%Trial Info, adapted from Trial Notes  (C:\Users\dsarma\Box
%Sync\2018_02_28_L7DRG_)

trialNum = {'Trial 1 – ';...
    'Trial 2 – ';...
    'Trial 3 – ';...
    'Trial 6 – ';...
    'Trial 7 – ';...
    'Trial 8 – ';...
    'Trial 9 – ';...
    'Trial 10 – ';...
    'Trial 11 – ';...
    'Trial 12 – ';...
    'Trial 13 – '};
trialInfo = {'toe squeezing: 5 reps/toe (starting @ lateral toe, moving in & finishing w/ central pad) ';...
    'mtp flexion/extension: 5 reps/toe (gap between each)';...
    'stroking shank: 10 reps each (proxim to dist dx anterior side) + (post side, dist-proxim)';...
    'rob poking around at hair follicle receptor behind knee ';...
    '10 reps of ramp and hold flexion/extension at knee';...
    '10 reps of ramp and hold flexion/extension at hip';...
    'repeat trial 1 to test consistency';...
    'repeat trial 2 to check consistency';...
    'ankle flexion/extension ramp and hold 10 times';...
    'knee';...
    'hip'};

%Master Test Condition struct
Condition{1,length(trialInfo)} = [];

%Load data folder
disp('Please Select the Data Folder');
[filenames, pathname] = uigetfile('*rhd','Pick files','MultiSelect', 'on');

if length(filenames) ~= length(trialInfo);
    disp('You have the wrong number of files for this experiment set!!!!');
    return;
end

%Filtering Parameters
hpFilt = designfilt('highpassiir', ...       % Response type
       'StopbandFrequency',400, ...     % Frequency constraints
       'PassbandFrequency',550, ...
       'StopbandAttenuation',55, ...    % Magnitude constraints
       'PassbandRipple',4, ...
       'DesignMethod','cheby1', ...     % Design method
       'MatchExactly','stopband', ...   % Design method options
       'SampleRate',30000)               % Sample rate
   
%Load and Transform Data Files
for i=length(filenames)
    Condition{i}.trial = string(trialNum(i));
    Condition{i}.info = string(trialInfo(i));
    Condition{i}.file = filenames(i);
    %Load Data
    read_IntanRHD_v2(pathname, filenames{i});
    %Filter Data
    dataF = filtfilt(hpFilt, amplifier_data')';
    %Assign Shanks
    Condition{i}.timeVec = t_amplifier;
    Condition{i}.fs = frequency_parameters.amplifier_sample_rate;
    Condition{i}.shankA = flipud(dataF(17:24,:));
    Condition{i}.shankB = flipud(dataF(9:16,:));
    Condition{i}.shankC = dataF(25:32,:);
    Condition{i}.shankD = dataF(1:8,:);
end


%Plotting Correlations & Covariance
% disp('Plotting Correlations....');
% for i=length(filenames)
%     disp(Condition{i}.info);
%     figure;
%     subplot(2,2,1)
%     imagesc(corr(Condition{i}.shankA')); colorbar; title('Shank D - Correlation'); 
%     subplot(2,2,2)
%     imagesc(corr(Condition{i}.shankB')); colorbar; title('Shank D - Correlation');
%     subplot(2,2,3)
%     imagesc(corr(Condition{i}.shankC')); colorbar; title('Shank D - Correlation');
%     subplot(2,2,4)
%     imagesc(corr(Condition{i}.shankD')); colorbar; title('Shank D - Correlation');
%     suptitle({[Condition{i}.trial],[Condition{i}.info]});
%     maximize;
%     
% end
% 
% disp('Plotting Covariance....');
% for i=1:length(filenames)
%     disp(Condition{i}.info);
%     figure;
%     subplot(2,2,1)
%     imagesc(cov(Condition{i}.shankA')); colorbar; title('Shank A - Covariance'); 
%     subplot(2,2,2)
%     imagesc(cov(Condition{i}.shankB')); colorbar; title('Shank B - Covariance');
%     subplot(2,2,3)
%     imagesc(cov(Condition{i}.shankC')); colorbar; title('Shank C - Covariance');
%     subplot(2,2,4)
%     imagesc(cov(Condition{i}.shankD')); colorbar; title('Shank D - Covariance');
%     suptitle({[Condition{i}.trial],[Condition{i}.info]});
%     maximize;
% end

%% Plotting Correlations by Shank
disp('Plotting Correlations by Shank');
disp('...Shank A...')
figure;
for i=1:length(filenames)
    disp(Condition{i}.info);
    subplot(round(length(trialInfo)/3),floor(length(trialInfo)/3),i)
    imagesc(corr(Condition{i}.shankA')); colorbar; title({[Condition{i}.trial],[Condition{i}.info]});
end
suptitle('Shank A - Correlations for All Conditions');
maximize;

disp('...Shank B...')
figure;
for i=1:length(filenames)
    disp(Condition{i}.info);
    subplot(round(length(trialInfo)/3),floor(length(trialInfo)/3),i)
    imagesc(corr(Condition{i}.shankB')); colorbar; title({[Condition{i}.trial],[Condition{i}.info]});
end
suptitle('Shank B - Correlations for All Conditions');
maximize;

disp('...Shank C...')
figure;
for i=1:length(filenames)
    disp(Condition{i}.info);
    subplot(round(length(trialInfo)/3),floor(length(trialInfo)/3),i)
    imagesc(corr(Condition{i}.shankC')); colorbar; title({[Condition{i}.trial],[Condition{i}.info]});
end
suptitle('Shank C - Correlations for All Conditions');
maximize;

disp('...Shank D...')
figure;
for i=1:length(filenames)
    disp(Condition{i}.info);
    subplot(round(length(trialInfo)/3),floor(length(trialInfo)/3),i)
    imagesc(corr(Condition{i}.shankD')); colorbar; title({[Condition{i}.trial],[Condition{i}.info]});
end
suptitle('Shank D - Correlations for All Conditions');
maximize;


%% Plotting Covariances by Shank
disp('Plotting Covariance by Shank');
disp('...Shank A...')
figure;
for i=1:length(filenames)
    disp([Condition{i}.trial Condition{i}.info]);
    subplot(round(length(trialInfo)/3),floor(length(trialInfo)/3),i)
    imagesc(cov(Condition{i}.shankA')); colorbar; title({[Condition{i}.trial],[Condition{i}.info]});
end
suptitle('Shank A - Covariances for All Conditions');
maximize;

disp('...Shank B...')
figure;
for i=1:length(filenames)
    disp([Condition{i}.trial Condition{i}.info]);
    subplot(round(length(trialInfo)/3),floor(length(trialInfo)/3),i)
    imagesc(cov(Condition{i}.shankB')); colorbar; title({[Condition{i}.trial],[Condition{i}.info]});
end
suptitle('Shank B - Covariances for All Conditions');
maximize;

disp('...Shank C...')
figure;
for i=1:length(filenames)
    disp([Condition{i}.trial Condition{i}.info]);
    subplot(round(length(trialInfo)/3),floor(length(trialInfo)/3),i)
    imagesc(cov(Condition{i}.shankC')); colorbar; title({[Condition{i}.trial],[Condition{i}.info]});
end
suptitle('Shank C - Covariances for All Conditions');
maximize;

disp('...Shank D...')
figure;
for i=1:length(filenames)
    disp([Condition{i}.trial Condition{i}.info]);
    subplot(round(length(trialInfo)/3),floor(length(trialInfo)/3),i)
    imagesc(cov(Condition{i}.shankD')); colorbar; title({[Condition{i}.trial],[Condition{i}.info]});
end
suptitle('Shank D - Covariances for All Conditions');
maximize;

