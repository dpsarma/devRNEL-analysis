artifactBuffer2 = 0.005; % post stim artifact settle lag in s, so 0.01 is 10ms
rmsWindow2 = 0.100; % post stim RMS Window in s, so 0.05 is 50ms
% % Set Epoch Window for Stim-triggered Averaging
% dtPre = 10e-3; % msec before stim onset 
% dtPost = 10e-3; % msec after stim onset 
% nSampsPre = floor(dtPre*fs);
% nSampsPost = floor(dtPost*fs);
iter = 2; %iteration count to reduce number of stim-trig average plots

for fnum = 1:length(emgFilenames)
    
    emg = trial(fnum).emg;
    emgF = trial(fnum).emgF;
    timeVec = trial(fnum).timeVec; 
    
    stimLength = trial(fnum).stimLength;
     stims = trial(fnum).stims;
     stimBegins = stims(1:2:end);
     stimEnds = stims(2:2:end);
    
    %% Epoch data into Windows 
    for i = 1:size(emg,1)
        asd = trial(fnum).emgF(i,:);
        emgTrResp{i} = asd(stims(end)+(artifactBuffer2*fs):(stims(end)+(artifactBuffer2*fs)+(rmsWindow2*fs)));
    end
       
    %% Get RMS acros epochs
    trial(fnum).epochTrRMS = cellfun(@rms,emgTrResp,'UniformOutput', false);
    trial(fnum).epochTrP2P = cellfun(@peak2peak,emgTrResp,'UniformOutput', false);
    %% Extract "Mean" Trace & Response
    stimEnd = nSampsPre + stimLength + artifactBuffer*fs;    
    for m = 1:length(Muscle_ID)
        % Response
        trial(fnum).responseTrRMS(m,:) = max(trial(fnum).epochTrRMS{m});
        trial(fnum).responseTrP2P(m,:) = max(trial(fnum).epochTrP2P{m});
    end
    
    clearvars -except emgFilenames emgPathname reportPath Muscle_ID fs b a bS aS ...
                dtPre dtPost fnum spinalElec nevStimCh stimAmp pulseWidth...
                stimDuration stimFrequency responseRMS trialName ...
                trialPathname trialFilenames trial subjectName setName ...
                setDescrpt setPath artifactBuffer rmsWindow ...
                externalStimulator nSampsPre nSampsPost chanPort C iter artifactBuffer2 rmsWindow2
end



% % %% Plotting Stim Averages
[sortedStim,I] = sort(stimAmp);
[sortedPulse,I2] = sort(pulseWidth);
[sortedFreq,I3] = sort(stimFrequency);

lastSamp = fs*3;

for c = 1:length(Muscle_ID)
    disp('Plotting Stim Averages');
    
    hF(c) = figure; %('visible', 'off');
    tmpsup = suptitle([Muscle_ID{c} ' - ' setDescrpt ' -  Stim Traces']);
    set(tmpsup, 'Interpreter', 'none') 
    maximize;
    
    subCount = 1; 
    for fnum = 1:iter:length(emgFilenames)
        hFs(fnum) = subplot(round(length(emgFilenames)/iter)+1,1,subCount);
%         hFs(fnum) = subplot(10,1,subCount);
        %Plot /Channel
        plot(trial(I(fnum)).timeVec(1:lastSamp),trial(I(fnum)).emgF(c,1:lastSamp));hold on;
%         plot(trial(I(fnum)).epochTimeVec,trial(I(fnum)).meanTrace{c},'k','LineWidth',1.2);
        vline([trial(I(fnum)).stims(1)/fs    trial(I(fnum)).stims(end)/fs],'r-');
        tmptit = title([trial(I(fnum)).Name ' - ' num2str(stimAmp(I(fnum))/1000) 'mA, ' ...
            num2str(pulseWidth(I(fnum))) 'us, ' num2str(stimFrequency(I(fnum))) 'Hz, ' num2str(stimDuration(I(fnum))) 's']); 
        set(hFs(fnum),'FontSize',7);
        set(tmptit, 'Interpreter', 'none','FontSize',5);
        hFs(fnum).XLim = [0 trial(I(fnum)).timeVec(lastSamp)];
        
        subCount = subCount + 1;
    end
    linkaxes([hFs],'xy');
    suplabel('Amplitude(uV)','y');
    suplabel('Time(sec)','x');
%     set(hF(c),'Position', [1343 49 570 947]);
    
    disp('Saving to Set Folder')
    savefig(hF(c),[setPath Muscle_ID{c} '_filtEMGTraces_' setName '_' setDescrpt]);
    saveas(hF(c),[setPath Muscle_ID{c} '_filtEMGTraces_' setName '_' setDescrpt '.png']);
end

close all


%% Plotting Recruitment Curves

% % %% Plotting Stim Averages
[sortedStim,I] = sort(stimAmp);
[sortedPulse,I2] = sort(pulseWidth);
[sortedFreq,I3] = sort(stimFrequency);

% response = [trial(:).responseTrRMS]';
response = [trial(:).responseTrP2P]';

switch length(Muscle_ID)
    case 16
        pRow = 2; pCol = 1;
        mCount = length(Muscle_ID)/2;
        pIdx = 2;
    case 8
        pRow = 1; pCol = 3;
        mCount = length(Muscle_ID);
        pIdx = 1;
    otherwise
        disp('How Would you like to Plot the Recruitment Curves?');
        return;
end

hF3 = figure; %maximize;
for pI = 1:pIdx
    
    hF3s(pI) = subplot(pRow,pCol, pCol*(pI-1)+1);
    for m = 1+mCount*(pI-1):mCount+mCount*(pI-1)
        yyA = smooth(sortedStim/1000,response(I,m),'sgolay');
        plot(sortedStim/1000,yyA,'Color',C(m,:),'LineWidth',3);
%         scatter(stimAmp/1000, response(I,m)');
        hold on;
    end
        ylabel('(post-Stim) P2P  (uV)')
        xlabel('Stimulation Amplitude (mA)')
        title(['By Amplitude']);

%     subplot(pRow,pCol, pCol*(pI-1)+2)
%     for m = 1+mCount*(pI-1):mCount+mCount*(pI-1)
% %         scatter(pulseWidth,response(:,m)')
%         yyP = smooth(sortedPulse,response(I2,m),'sgolay');
%         plot(sortedPulse,yyP,'Color',C(m,:),'LineWidth',3);
%         hold on;
%     end
%         ylabel('(post-Stim) Mean RMS  (uV)')
%         xlabel('Stimulation PulseWidth (us)')
%         title('By PulseWidth');
%     %     
%     subplot(pRow,pCol, pCol*(pI-1)+3)
%     for m = 1+mCount*(pI-1):mCount+mCount*(pI-1)
% %         scatter(stimFrequency,response(:,m)')
%         yyF = smooth(sortedFreq,response(I3,m),'sgolay');
%         plot(sortedFreq,yyF,'Color',C(m,:),'LineWidth',3);
%         hold on;
%     end
%         ylabel('(post-Stim) Mean RMS  (uV)')
%         xlabel('Stimulation Frequency (Hz)')
%         title('By Frequency');
        
    if pI == 1
        legend([Muscle_ID(1:mCount)],'Location','northwest', 'Orientation','vertical');
    elseif pI == 2
        legend([Muscle_ID((mCount)+1:length(Muscle_ID))],'Location','northwest', 'Orientation','vertical');
    end
%     
end
tmptit = suptitle([setName ' - ' setDescrpt ' -  Recruitment Curves (P2P)']);
set(tmptit, 'Interpreter', 'none');
% % Pix_SS = get(0,'screensize');
hF3.Position = [680 49 572 947];
% linkaxes([hF3s],'xy');
    
saveas(hF3,[reportPath '\' strrep(datestr(now),':','_')  '_TrainRecruitmentCurvesP2P_2_' setName '_' setDescrpt  '.png']);
savefig(hF3,[reportPath '\' strrep(datestr(now),':','_') '_TrainRecruitmentCurvesP2P_2_' setName '_' setDescrpt ]);


%% RMS RC
response = [trial(:).responseTrRMS]';
% response = [trial(:).responseTrP2P]';

switch length(Muscle_ID)
    case 16
        pRow = 2; pCol = 1;
        mCount = length(Muscle_ID)/2;
        pIdx = 2;
    case 8
        pRow = 1; pCol = 3;
        mCount = length(Muscle_ID);
        pIdx = 1;
    otherwise
        disp('How Would you like to Plot the Recruitment Curves?');
        return;
end

hF3 = figure; %maximize;
for pI = 1:pIdx
    
    hF3s(pI) = subplot(pRow,pCol, pCol*(pI-1)+1);
    for m = 1+mCount*(pI-1):mCount+mCount*(pI-1)
        yyA = smooth(sortedStim/1000,response(I,m),'sgolay');
        plot(sortedStim/1000,yyA,'Color',C(m,:),'LineWidth',3);
%         scatter(stimAmp/1000, response(I,m)');
        hold on;
    end
        ylabel('(post-Stim) RMS  (uV)')
        xlabel('Stimulation Amplitude (mA)')
        title(['By Amplitude']);

%     subplot(pRow,pCol, pCol*(pI-1)+2)
%     for m = 1+mCount*(pI-1):mCount+mCount*(pI-1)
% %         scatter(pulseWidth,response(:,m)')
%         yyP = smooth(sortedPulse,response(I2,m),'sgolay');
%         plot(sortedPulse,yyP,'Color',C(m,:),'LineWidth',3);
%         hold on;
%     end
%         ylabel('(post-Stim) Mean RMS  (uV)')
%         xlabel('Stimulation PulseWidth (us)')
%         title('By PulseWidth');
%     %     
%     subplot(pRow,pCol, pCol*(pI-1)+3)
%     for m = 1+mCount*(pI-1):mCount+mCount*(pI-1)
% %         scatter(stimFrequency,response(:,m)')
%         yyF = smooth(sortedFreq,response(I3,m),'sgolay');
%         plot(sortedFreq,yyF,'Color',C(m,:),'LineWidth',3);
%         hold on;
%     end
%         ylabel('(post-Stim) Mean RMS  (uV)')
%         xlabel('Stimulation Frequency (Hz)')
%         title('By Frequency');
        
    if pI == 1
        legend([Muscle_ID(1:mCount)],'Location','northwest', 'Orientation','vertical');
    elseif pI == 2
        legend([Muscle_ID((mCount)+1:length(Muscle_ID))],'Location','northwest', 'Orientation','vertical');
    end
%     
end
tmptit = suptitle([setName ' - ' setDescrpt ' -  Recruitment Curves (RMS)']);
set(tmptit, 'Interpreter', 'none');
% % Pix_SS = get(0,'screensize');
hF3.Position = [680 49 572 947];
% linkaxes([hF3s],'xy');
    
saveas(hF3,[reportPath '\' strrep(datestr(now),':','_')  '_TrainRecruitmentCurvesRMS_2_' setName '_' setDescrpt  '.png']);
savefig(hF3,[reportPath '\' strrep(datestr(now),':','_') '_TrainRecruitmentCurvesRMS_2_' setName '_' setDescrpt ]);

