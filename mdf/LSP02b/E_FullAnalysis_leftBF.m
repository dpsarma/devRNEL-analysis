%% UH3 - LSP02 anaylsis (5) - Looking at JUST 1 Muscle (left BF)
%----------------------------------------
% tag that is missing from all EMG trials is "Sitting" vs. "Standing"
% RUN THIS IN SECTIONS!!!!!!!
omdfc = mdf.init('auto'); %% Startup MDF


%% Environment Setup
summary =  mdf.load('mdf_type','summary');

datasave = 'R:\users\dsarma\LSP02b\matData\mdftest\';
figuresave = 'R:\users\dsarma\LSP02b\Figures\mdftest_leftBF\';



%% MDF Trial Labels Setup
datestrings = ['2018-10-08'; '2018-10-09'; '2018-10-11'; '2018-10-16'; '2018-10-19'; '2018-10-24'];
dayIDs = ['Day 02'; 'Day 03'; 'Day 04'; 'Day 07'; 'Day 09'; 'Day 12'];
musclabels = summary.emgLabel;
electrodenums = [1:48];
freqs = [1 2 5];%Hz
pulses = [200 250 400 500 600 750 800 1000];%us
for i=1:length(electrodenums)
            if (electrodenums(i)<10)
                electrodestrings{i} = ['Unipolar: 0' num2str(electrodenums(i))];
            else
                electrodestrings{i} = ['Unipolar: ' num2str(electrodenums(i))];
            end 
end

%% Specify Load List
freqstring = freqs(1);
pulsestring = pulses(1);
m = 12;
% Use this to automate for future automated call

%% Load Relevant UUID_List 
load('R:\users\dsarma\LSP02b\matData\mdftest\uuid lists\mdf_uuidlist_1Hz_200us.mat'); %NOTE THIS RELOADS THE MDF TRIAL LABEL (load in to a variable to ID specifc variables)

% broke at 2,35,28
%% Load Data & Parse
for m = 1:length(musclabels)
for f = 4%1:size(datestrings,1) %All Days
    for e = [6 33]%1:length(electrodenums) %All Electrodes
        disp(['Starting: ' char(musclabels(m)) ' ' dayIDs(f,:) ' ' electrodestrings(e)])
        if isempty(uuid_list{e,f})
                disp(['No Trials for: ' dayIDs(f,:) ', ' electrodestrings(e)]);
                continue;
        end
        % Load & Parse Trial EMG
        for tnum = 1:length(uuid_list{e,f}) %All Trials for Given Electrode & Day
            
            %Load Trial & Mask EMG
            file = mdf.load(char(uuid_list{e,f}(tnum)));
            disp(['Processing: ' file.trialname]);
            
            if (f == 2 && e == 35 && tnum == 28) || (f == 4 && e == 33 && tnum == 52)
                disp(['Processing: ' file.trialname 'missing StimT']);
                continue;
            end
            
            zz= file.EMG(m); 
            emgraw =  zz.data.wf; emgtime = zz.data.time;
            
            % Bandpass Data
            fs = 30e3; flow = 70; fhigh = 7500; %alternatively [10,7500]
            emgfiltd = bandpass(emgraw,[flow,fhigh],fs);
            
            %Rectify & Smooth Data
            emgrect = abs(emgfiltd); %emgrectsm = smooth(emgrect);
            
            % Select Stims & Params
            if length(file.stim.stimT) > 15
                stimons = floor(file.stim.stimT(1:2:end)*fs);%time to sample#
                stimoffs = floor(file.stim.stimT(2:2:end)*fs);
            else
                stimons = floor(file.stim.stimT*fs);
            end
            
            % Assign Stim Params
            stimamp = file.stimParams.cathAmp/1000;%mA
            stimfreq = file.stimParams.freq;%Hz
            stimpulse = file.stimParams.width/1000/1000;%us
            
            % Set Epoch Window for Stim-triggered Averaging
            dtPre = 50e-3; % msec before stim onset 
            dtPost = 100e-3; % msec after stim onset 
            nSampsPre = floor(dtPre*fs);
            nSampsPost = floor(dtPost*fs);
            stimlength = floor(2*stimpulse*fs)+2;

            epochtime = linspace(-dtPre, (stimpulse) + dtPost, nSampsPre + nSampsPost + stimlength );

            % Epoch Data into Segments
            disp(['Epoching: ' file.trialname])
            emgepochs = cell2mat(arrayfun(@(x) emgfiltd((stimons(x)-nSampsPre):(stimons(x)+(stimlength+nSampsPost-1))), 1:length(stimons), 'UniformOutput',false)');
            emgepochs_rect = cell2mat(arrayfun(@(x) emgrect((stimons(x)-nSampsPre):(stimons(x)+(stimlength+nSampsPost-1))), 1:length(stimons), 'UniformOutput',false)');
            
            % Find Mean 
            meantrace = mean(emgepochs);
            meantrace_rect = mean(emgepochs_rect);
            
            % Calculate Values for Recruitment Curve and Onsets
            rmswindow = floor(50e-3 * fs); %in ms
            blanksamps = floor(5e-3 * fs); %in ms post stim
            windsamps = nSampsPre + stimlength + blanksamps;

            epochRMS = rms(emgepochs_rect(:,windsamps:windsamps+rmswindow)')'; 
            epochP2P = peak2peak(emgepochs(:,windsamps:windsamps+rmswindow)')';  

            meanRMS = mean(epochRMS); errRMS = std(epochRMS)/sqrt(length(epochRMS));
            meanP2P = mean(epochP2P); errP2P = std(epochP2P)/sqrt(length(epochP2P));
            
            % Calculate Onset Values
            [pks,locs,widths,proms] = findpeaks(emgepochs_rect(windsamps:windsamps+rmswindow), epochtime(windsamps:windsamps+rmswindow), 'Annotate','extents', 'WidthReference','halfheight'); 
            [mpks,mlocs, mwidths, mproms] = findpeaks(meantrace_rect(windsamps:windsamps+rmswindow), epochtime(windsamps:windsamps+rmswindow), 'Annotate','extents', 'WidthReference','halfheight'); 
            
            %% Hash Values w/o object structure
            hash{e,f,tnum}.epochs = emgepochs;
            hash{e,f,tnum}.epochsrect = emgepochs_rect;
            hash{e,f,tnum}.epochtime = epochtime;
            hash{e,f,tnum}.meantrace = meantrace;
            hash{e,f,tnum}.meantrace_rect = meantrace_rect;
            hash{e,f,tnum}.epochRMS = epochRMS;
            hash{e,f,tnum}.epochP2P = epochP2P;
            hash{e,f,tnum}.meanRMS = meanRMS;
            hash{e,f,tnum}.meanP2P = meanP2P;
            hash{e,f,tnum}.errRMS = errRMS;
            hash{e,f,tnum}.errP2P = errP2P;
            hash{e,f,tnum}.peaks = [pks,locs,widths,proms];
            hash{e,f,tnum}.mpeaks = [mpks,mlocs, mwidths, mproms];

            hash{e,f,tnum}.stimAmp = stimamp;
            hash{e,f,tnum}.stimfreq = stimfreq;
            hash{e,f,tnum}.stimpulse = stimpulse;
        end
    end
end
save(char([datasave char(musclabels(m)) '_AllDays_AllElectrodes_' num2str(freqstring) 'Hz_' num2str(pulsestring) 'us.mat']),'hash','-v7.3');
clear hash
end
%% Plotting Stage