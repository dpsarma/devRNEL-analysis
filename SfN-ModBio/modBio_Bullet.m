%%ModBio Analysis for Bullet Data
%R:\data_raw\primate\NHP_DRG_Stim\Bullet_Data\Grapevine

%Load data folder
% disp('Please Select the Data Folder');
% [filenames, pathname] = uigetfile('*ns5','Pick files','MultiSelect', 'on');
pathname = 'R:\data_raw\primate\NHP_DRG_Stim\Bullet_Data\Grapevine\';

%Using only files that have .ns5 with raw channels for probes.
tmp = load('BulletRawDataFiles.mat');
filenames = tmp.Rnew;

trialInfo = {'41 - Shaker, RF Center, Full';...
'42 - Shaker, 2-06, RF Center, Full';...
'43 - Shaker, 2-06, 1mm distal, Full';...
'44 - Shaker, 2-06, 2mm distal, Full';...
'45 - Shaker, 2-06, 3mm distal, Full';...
'46 - Shaker, 2-06, 1mm proximal, Full';...
'47 - Shaker, 2-06, 2mm proximal, Full';...
'48 - Shaker, 2-06, 3mm proximal, Full';...
'49 - Shaker, 2-31, RF Center, Full';...
'50 - GARBAGE NOISE';...
'51 - Shaker, 2-31, 1mm proximal, Full';...
'52 - CRASHED, Shaker, 2-31, 2mm proximal, Full';...
'53 - Shaker, 2-31, 2mm proximal, Full';...
'54 - Shaker, 2-31, 3mm proximal, SHORT';...
'55 - Shaker, 2-31, 4mm proximal, SHORT';...
'56 - Shaker, 2-31, 5mm proximal, SHORT';...
'57 - Shaker, 2-31, 6mm proximal, SHORT';...
'58 - Shaker, 2-16, RF Center, Full';...
'59 - Shaker, 2-16, 1mm proximal, SHORT';...
'60 - Shaker, 2-16, 2mm proximal, SHORT';...
'61 - Shaker, 2-16, 3mm proximal, SHORT';...
'62 - Shaker, 2-16, 4mm proximal, SHORT';...
'63 - Cancelled';...
'64 - Shaker, 2-16, 5mm proximal, SHORT';...
'65 - Shaker, 2-16, 6mm proximal, SHORT';...
'66 - Shaker, 2-16, 7mm proximal, SHORT';...
'67 - Shaker, 2-16, 8mm proximal, SHORT';...
'68 - Shaker, 2-16, 9mm proximal, SHORT';...
'69 - Shaker, 2-16, 1mm distal, SHORT';...
'70 - Shaker, 2-16, 2mm distal, SHORT';...
'71 - Shaker, 2-16, 3mm distal, SHORT';...
'72 - Shaker, 2-16, 4mm distal, SHORT';...
'73 - DRUM, 2-06, D4d, shallow';...
'74 - DRUM, 2-06, D4d, deep 300u';...
'75 - Cancelled';...
'76 - DRUM, 2-06, D4d, deep 300u';...
'77 - DRUM, 2-06, D4d, deep 300u';...
'78 - DRUM, 2-06, D4d, deep 600u'};

%Master Test dataFiles struct
dataFiles{1,length(filenames)} = [];

%Filtering Parameters
% hpFilt = designfilt('highpassiir', ...       % Response type
%        'StopbandFrequency',400, ...     % Frequency constraints
%        'PassbandFrequency',550, ...
%        'StopbandAttenuation',55, ...    % Magnitude constraints
%        'PassbandRipple',4, ...
%        'DesignMethod','cheby1', ...     % Design method
%        'MatchExactly','stopband', ...   % Design method options
%        'SampleRate',30000);               % Sample rate
%    
   fstart = 31;%31,18
   fstop = 37;%37,22
   probeChans = 32;
   numProbes = 2;
   
   for probeNum = 1:numProbes
       chanNumShift(probeNum) = 32*(probeNum-1);
   end
   
 %Load and Transform Data Files
 fcount = 0;
for fNum=fstart:fstop%length(filenames)

    dataFiles{fNum}.filename = filenames(fNum);
    dataFiles{fNum}.trialInfo = trialInfo(fNum);
    %Load Data
%     tmpName =fullfile(pathname, filenames(i));
    for probeNum =  1:numProbes
        disp(['Loading Probe Number - ' num2str(probeNum)]);
        [dataFiles{fNum}.shankData{probeNum},dataFiles{fNum}.timeVec] = read_continuousData([pathname filenames{fNum}], 'raw' , chanNumShift(probeNum)+[1:probeChans]);
    end
    fcount = fcount+1;
end


%% Plotting Correlations by Shank
disp('Plotting Correlations by Shank');

shankA = [1:8];
shankB = [9:16];
shankC = [17:24];
shankD = [25:32];
shanks = [shankA; shankB; shankC; shankD];
shankNames = ['Shank A'; 'Shank B'; 'Shank C'; 'Shank D'];
probe = [];

for probeNum = 1:numProbes
    disp(['Calculating for Probe Number - ' num2str(probeNum)]);
    for shankNum = 1:size(shanks,1)
        for j=fstart:fstop
            dataFiles{j}.probe(probeNum).shank(shankNum).corr = corr(dataFiles{j}.shankData{probeNum}([shanks(shankNum,:)],:)');
            dataFiles{j}.probe(probeNum).shank(shankNum).cov = cov(dataFiles{j}.shankData{probeNum}([shanks(shankNum,:)],:)');
        end
    end
end

for probeNum = 1:numProbes
    
    disp(['Showing Probe Number - ' num2str(probeNum)]);
    
    for shankNum = 1:size(shanks,1)
        disp(['...' shankNames(shankNum,:) '...']);
        
        figure;
        hold on;
        for j=1:fcount
%             if fcount < j
%                 figIdx = j - abs(j-fcount);
%             else
%                 figIdx = j;
%             end
            fNum = (fstart-1) + j;
            subplot(round(fcount/2),2,j); %figIdx
            imagesc(dataFiles{fNum}.probe(probeNum).shank(shankNum).corr); colorbar; title([dataFiles{fNum}.filename]); %trialInfo
            disp([dataFiles{fNum}.filename]);
        end
        suptitle({['Probe ' num2str(probeNum) ': ' shankNames(shankNum,:) ' - Correlations for subset "Shaker" dataFiles']});
        maximize;
    end
end

for probeNum = 1:numProbes
    
    disp(['Showing Probe Number - ' num2str(probeNum)]);
    
    for shankNum = 1:size(shanks,1)
        disp(['...' shankNames(shankNum,:) '...']);
        
        figure;
        hold on;
        for j=1:fcount
            fNum = (fstart-1) + j;
            subplot(round(fcount/2),2,j); %figIdx
            imagesc(dataFiles{fNum}.probe(probeNum).shank(shankNum).cov); colorbar; title([dataFiles{fNum}.filename]); %trialInfo
            disp([dataFiles{fNum}.filename]);
        end
        suptitle({['Probe ' num2str(probeNum) ': ' shankNames(shankNum,:) ' - Covariances for subset "Shaker" dataFiles']});
        maximize;
    end
end


