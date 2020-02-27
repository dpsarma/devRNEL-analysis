%% UH3 - LSP02 anaylsis (1) - MDF Data Load
%----------------------------------------
% tag that is missing from all EMG trials is "Sitting" vs. "Standing"

%% Example Processing for Single Electrodes at a time by Data Trial:
% Selecting Day 02, Electrode 9, Muscle Left ST
e = 9; f = 1; m = 13; tnum = 10;

for m = m%1:lengthmusclabels
    for f = f; %1:size(datestrings,1)
        for e = e; %1:length(electrodenums)

            if isempty(uuid_list{e,f})
                disp(['No Trials for: ' dayIDs(f,:) ', ' electrodestrings(e)]);
                continue;
            end

            for tnum = 1:length(uuid_list{e,f})

                % Select Data wf
                file = mdf.load(char(uuid_list{e,f}(tnum)));
                disp(['Processing: ' file.trialname])

                zz = file.EMG(m);
                emgraw =  zz.data.wf;
                emgtime = zz.data.time;

                % Bandpass Data
                fs = 30e3; flow = 10; fhigh = 7500; %was previously [75,7500]
                emgfiltd = bandpass(emgraw,[flow,fhigh],fs);

                %Rectify & Smooth Data
                emgrect = abs(emgfiltd);
                %emgrectsm = smooth(emgrect);

                % Select Stims & Params
                if length(file.stim.stimT) > 15
                stimons = floor(file.stim.stimT(1:2:end)*fs);%time to sample#
                stimoffs = floor(file.stim.stimT(2:2:end)*fs);
                else
                stimons = floor(file.stim.stimT*fs);
                end
                stimamp = file.stimParams.cathAmp/1000;%mA
                stimfreq = file.stimParams.freq;%Hz
                stimpulse = file.stimParams.width/1000/1000;%in us

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

                %%% THESE EPOCHS Could BE A CHILD OBJECT

                % Find Mean 
                meantrace = mean(emgepochs);
                meantrace_rect = mean(emgepochs_rect);
                %%% To Save

                % % % Verify by plotting
                % % figure;plot(epochtime,emgepochs); hold on; plot(epochtime,meantrace,'k');
                % % figure;plot(epochtime,emgepochs_rect); hold on; plot(epochtime,meantrace_rect,'k');
                % % figure;
                % % for i = 1:size(emgepochs,1)
                % %     subplot(15,1,i);
                % %     plot(epochtime,emgepochs(i,:)); hold on
                % %     plot(epochtime,meantrace);
                % % end

                % Calculate Values for Recruitment Curve and Onsets
                rmswindow = floor(50e-3 * fs); %in ms
                blanksamps = floor(5e-3 * fs); %in ms post stim
                windsamps = nSampsPre + stimlength + blanksamps;

                epochRMS = rms(emgepochs(:,windsamps:windsamps+rmswindow)')'; %%%% SAVE This???
                epochP2P = peak2peak(emgepochs(:,windsamps:windsamps+rmswindow)')';  %%%% SAVE This???

                meanRMS = mean(epochRMS); errRMS = std(epochRMS)/sqrt(length(epochRMS));
                meanP2P = mean(epochP2P); errP2P = std(epochP2P)/sqrt(length(epochP2P));

                %% Hash Values w/o object structure
                hash{e,f,tnum}.epochs = emgepochs;
                hash{e,f,tnum}.epochsrect = emgepochs_rect;
                hash{e,f,tnum}.epochRMS = epochRMS;
                hash{e,f,tnum}.epochP2P = epochP2P;
                hash{e,f,tnum}.meanRMS = meanRMS;
                hash{e,f,tnum}.meanP2P = meanP2P;
                hash{e,f,tnum}.errRMS = errRMS;
                hash{e,f,tnum}.errP2P = errP2P;

                hash{e,f,tnum}.stimAmp = stimamp;
                hash{e,f,tnum}.stimfreq = stimfreq;
                hash{e,f,tnum}.stimpulse = stimpulse;

            end
        save(char([datasave char(musclabels(13)) '_' dayIDs(f,:) '_E-' num2str(e) '.mat']),'hash','-v7.3');    
        end
%         save(char([datasave char(musclabels(13)) '_' dayIDs(f,:) '_allelectrodes.mat']),'hash','-v7.3');
    end
end

