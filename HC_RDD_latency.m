%% Load a trial struct
clear all
close all
clc

subj = 'LSP05';
%% Set Basic params
switch subj
    case 'LSP02b'
        load('C:\figs\UH3\HotelCali\LSP02b\LSP02b_PAD_Elec1_1Hz.mat');
        Muscle_ID = {'Right VM', 'Right RF', 'Right VL', 'Right BF', 'Right ST', ...
                    'Right TA', 'Right SO', 'Right LG', 'Left VM', 'Left RF', 'Left VL',...
                    'Left BF', 'Left ST', 'Left Ham', 'Left SO', 'Left LG'};
        fs = 30e3;

    case 'LSP05'
        load('C:\figs\UH3\HotelCali\LSP05\LSP05_PAD_Elec9_20Hz.mat')
        Muscle_ID = {'Right VM', 'Right RF', 'Right VL', 'Right BF', 'Right ST', ...
                    'Right TA', 'Right SO', 'Right LG','Left VM', 'Left RF','Left VL',...
                    'Left BF', 'Left ST', 'Left TA', 'Left MG', 'Left LG'};
        fs = 30e3;

    case 'LNP02'
        load('C:\figs\UH3\HotelCali\LNP02\LNP02_PAD_Elec4_20Hz.mat');
        Muscle_ID = {'Right VM', 'Right RF', 'Right VL', 'Right BF', 'Right ST', ...
                    'Right TA', 'Right SO', 'Right LG','Left VM', 'Left RF','Left VL',...
                    'Left BF', 'Left ST', 'Left TA', 'Left MG', 'Left LG'};
        fs = 7.5e3; % or 2e3;
end
freqs = unique([trial.stimFrequency]);
pwidths = unique([trial.pulseWidth]);
elecs = unique([trial.spinalElec]);
msubset = 1:16;

%% Set Windows
nyq = 0.5*fs;
artifactBuffer = 5e-3;
check_p = pwidths(1);
check_f = freqs(1);
switch check_f
        case 1
            dtPre = 50e-3; % msec before stim onset 
            dtPost = 100e-3; % msec after stim onset 
            rmsWindow = 0.05; % post stim RMS Window in s, so 0.05 is 50ms     
            if subj == 'LNP02'
                njit = 25e-3; 
            end
        case 2
            dtPre = 50e-3; % msec before stim onset 
            dtPost = 100e-3; % msec after stim onset 
            rmsWindow = 0.05; % post stim RMS Window in s, so 0.05 is 50ms
            if subj == 'LNP02'
                rmsWindow = 115e-3; %100ms
                dtPost = 150e-3;
                njit = 5e-3; 
            end
        case 5
            dtPre = 50e-3; % msec before stim onset 
            dtPost = 100e-3; % msec after stim onset 
            rmsWindow = 0.05; % post stim RMS Window in s, so 0.05 is 50ms
            if subj == 'LNP02'
                rmsWindow = 115e-3; %100ms
                dtPost = 150e-3;
                njit = 0e-3; 
            end
        case 10
            dtPre = 50e-3; % msec before stim onset 
            dtPost = 95e-3; % msec after stim onset 
            rmsWindow = 0.05; % post stim RMS Window in s, so 0.05 is 50ms
            if subj == 'LNP02'
                rmsWindow = 75e-3;
                njit = 0e-3; 
            end
        case 20
            dtPre = 50e-3; % msec before stim onset 
            dtPost = 45e-3; % msec after stim onset 
            rmsWindow = 0.040; % post stim RMS Window in s, so 0.05 is 50ms
            if subj == 'LNP02'
                njit = 0e-3; 
            end

        case 200
            dtPre = 50e-3; % msec before stim onset 
            dtPost = 4e-3; % msec after stim onset 
            rmsWindow = 0.050; % post stim RMS Window in s, so 0.05 is 50ms
        otherwise 
            dtPre = 10e-3; % msec before stim onset 
            dtPost = 20e-3; % msec after stim onset 
            rmsWindow = 0.005; % post stim RMS Window in s, so 0.05 is 50ms;
end
nSampsPre = floor(dtPre*fs);
nSampsPost = floor(dtPost*fs);
rmsWindow = floor(rmsWindow*fs);

if subj == 'LNP02'
    tmptrial = trial; clear trial;
%     trial.epochTimeVec = tmptrial.epochTimeVec;
    for fnum = 1:length(tmptrial)
        trial(fnum).spinalElec = tmptrial(fnum).spinalElec;
        trial(fnum).stimAmp = tmptrial(fnum).stimAmp;
        trial(fnum).stims = tmptrial(fnum).stims;
        trial(fnum).emgF = tmptrial(fnum).emgF;
        freq = tmptrial(fnum).stimFrequency;
        stimBegins = tmptrial(fnum).stims;
        stimBegins=stimBegins+floor(njit*fs);%Adjusting for nomad jitter
        stimLength = (tmptrial(fnum).pulseWidth)*fs + 2; %Pulse Width in ms, interstim interval is 2 ticks of clock (2/fs)
        trial(fnum).stimLength = stimLength;
        stimEnd = nSampsPre + stimLength + floor(artifactBuffer*fs); 
        
        trial(fnum).epochTimeVec = linspace(-dtPre, (trial(fnum).stimLength/fs) + dtPost, nSampsPre + nSampsPost + trial(fnum).stimLength);
        %% Epoch data into Windows 
        disp(['Epoching' num2str(fnum)]);
        for i = 1:size(tmptrial(fnum).emg,1)
            asd = tmptrial(fnum).emgF(i,:);
            emgStim{i} = cell2mat(arrayfun(@(x) asd((stimBegins(x)-nSampsPre):(stimBegins(x)+(stimLength+nSampsPost-1))), 1:length(stimBegins), 'UniformOutput',false)');
        end  
        trial(fnum).emgStim = emgStim;
        trial(fnum).baseline = tmptrial(fnum).emgF(:,1:stimBegins(1));
        
        %% Get RMS across epochs
        disp('Get Response')
        trial(fnum).epochRMS = cellfun(@rms,emgStim,'UniformOutput', false);
        trial(fnum).epochP2P = cellfun(@peak2peak,emgStim,'UniformOutput', false);
        for m = 1:length(Muscle_ID)
            % Mean Trace        
            trial(fnum).meanTrace{m} = mean(emgStim{1,m});
            trial(fnum).meanbase(m,:) = mean(trial(fnum).baseline(m,:));
            % Response
            if freq > 20
                trial(fnum).responseRMS(m) = mean(rms(tmptrial(fnum).emgF(m,stimBegins(end):stimBegins(end)+(rmsWindow))));
                trial(fnum).responseP2P(m) = mean(peak2peak(tmptrial(fnum).emgF(m,stimBegins(end):stimBegins(end)+(rmsWindow))));
            else
                trial(fnum).responseRMS(m) = mean(trial(fnum).epochRMS{m}(stimEnd:stimEnd+(rmsWindow)));
                trial(fnum).responseP2P(m) = mean(trial(fnum).epochP2P{m}(stimEnd:stimEnd+(rmsWindow)));
            end
        end
    end
end
clear tmptrial;

disp(['Working On: ' num2str(check_p) 'ms ' num2str(check_f) 'Hz']);

for i = 1:length(msubset)
    m = msubset(i);

    for t = 1:length(trial)
        disp(['Processing ' num2str(trial(t).stimAmp) 'mA - ' Muscle_ID(m)]);
        
        stimBegins = trial(t).stims;
        if subj == 'LNP02'
            stimBegins=stimBegins+floor(njit*fs); %Adjusting for nomad jitter
        end
        stimEnd = nSampsPre + trial(t).stimLength + floor(artifactBuffer*fs);
        meta(t).elec = elecs;
        meta(t).freq = check_f;
        meta(t).pulse = check_p; 
        meta(t).Amp = trial(t).stimAmp;
        meta(t).muscle(m).Muscle_ID = Muscle_ID{m};
        meta(t).muscle(m).meanTrace = trial(t).meanTrace{1, m};
        meta(t).muscle(m).epochTimeVec = trial(t).epochTimeVec;

        
        figure; nexttile; %maximize;
        plot(trial(t).epochTimeVec, trial(t).meanTrace{1, m},'LineWidth',2);
        hold on;
        plot(trial(t).epochTimeVec, abs(trial(t).meanTrace{1, m}),'LineWidth',2);
        vline([0.02 0.025]);
        title(strcat('e',num2str(trial(t).spinalElec),': ', num2str(trial(t).stimAmp), 'mA -', Muscle_ID(m)));

        pause(0.01);
        prompt = 'Is there a response? Y/N [Y]: ';
        x = input(prompt,'s');
        if contains(x,'n','IgnoreCase',true)
            meta(t).muscle(m).exist = 0;
            meta(t).muscle(m).p2pResponse = NaN;
            meta(t).muscle(m).peakResponse = NaN;
            meta(t).muscle(m).peak_latency = NaN; 
            meta(t).muscle(m).onset = NaN;
            meta(t).muscle(m).offset =  NaN;
            meta(t).muscle(m).width = NaN;
            close;
        else
            meta(t).muscle(m).exist = 1;
            meta(t).muscle(m).p2pResponse = peak2peak(trial(t).meanTrace{m}(stimEnd:end));
            [meta(t).muscle(m).peakResponse, mi] = max(abs(trial(t).meanTrace{m}(stimEnd:end)));
            meta(t).muscle(m).peak_latency = mi/fs; %trial(t).epochTimeVec((stimEnd-1) + mi);
            
            tmp = abs(trial(t).meanTrace{m}(stimEnd:(stimEnd-1) + mi));
            idx = find(tmp <= 0.25*meta(t).muscle(m).p2pResponse,1, 'last');
%             idx = find(abs(trial(t).meanTrace{m}(stimEnd:end)) == 0.25*(meta(t).muscle(m).p2pResponse));
            meta(t).muscle(m).onset = idx/fs;
            tmp = abs(trial(t).meanTrace{m}((stimEnd-1) + mi:end));
            idx = find(tmp <= 0.25*meta(t).muscle(m).p2pResponse,1, 'first');
            meta(t).muscle(m).offset = (mi+idx)/fs;
            meta(t).muscle(m).width = meta(t).muscle(m).offset - meta(t).muscle(m).onset;                    
            
            disp('Plotting Trace of First Two Stims Trace');
            disp(num2str([stimBegins(1:2)]))
            for s = 1:2
                timbegin = trial(t).stims(s)+floor(artifactBuffer*fs);
                meta(t).muscle(m).RDDp2p(s) = peak2peak(trial(t).emgF(m,timbegin:timbegin+rmsWindow));
                [meta(t).muscle(m).RDDpeakR(s),latency] = max(abs(trial(t).emgF(m,timbegin:timbegin+rmsWindow)));
                meta(t).muscle(m).RDDpkLat(s) = (latency/fs)*1000; %artifactBuffer
                
%                 onset = find(abs(trial(t).emgF(m,timbegin:timbegin+latency)) == 0.5*meta(t).muscle(m).RDDp2p(s));
                tmp = abs(trial(t).emgF(m,(stimEnd-1) + mi:end));
                onset = find(tmp <= 0.25*meta(t).muscle(m).RDDp2p(s),1, 'last');
                meta(t).muscle(m).RDDpkOns(s) = (onset/fs)*1000; %+ artifactBuffer
            end  

            close all;
            meta(t).muscle(m).first2trace = trial(t).emgF(m,stimBegins(1)-2*nSampsPre:stimBegins(2)+2*nSampsPost);
            meta(t).muscle(m).first2timVec = linspace(-2*nSampsPre/fs, (stimBegins(2)-stimBegins(1)+2*nSampsPost)/fs, 2*nSampsPre+stimBegins(2)-stimBegins(1)+2*nSampsPost+1);
            figure; plot(meta(t).muscle(m).first2timVec, meta(t).muscle(m).first2trace); hold on; vline([0 (stimBegins(2)-stimBegins(1))/fs]);
        end

        disp([num2str(t) ' of ' num2str(length(trial))]);
    end
end

reportPath = ['C:\figs\UH3\HotelCali\' subj '\'];
save([reportPath '\' subj '_metaB_e' num2str(elecs) '_' num2str(check_f) 'Hz.mat'], 'meta', '-v7.3');
