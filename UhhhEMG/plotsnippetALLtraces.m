% % %% Plotting Stim Averages
[sortedStim,I] = sort(stimAmp);
[sortedPulse,I2] = sort(pulseWidth);
[sortedFreq,I3] = sort(stimFrequency);

uniqFreqs = unique(stimFrequency);

for count =length
    freqIndices = find(stimFrequency == uniqFreqs(count))
end

% uniqPW = unique(pulseWidth);
% for count =length
%     pwIndices = find(pulseWidth == uniqPW(count))
% end



for c = 1:length(Muscle_ID)
    disp('Plotting Stim Averages');
    
    hF(c) = figure; %('visible', 'off');
    tmpsup = suptitle([Muscle_ID{c} ' - ' setDescrpt ' -  Stim Traces']);
    set(tmpsup, 'Interpreter', 'none') 
    maximize;
    
    subCount = 1; 
    for fnum = 1:iter:length(emgFilenames)
        hFs(fnum) = subplot(6,5,subCount);
%         hFs(fnum) = subplot(10,1,subCount);
        %Plot /Channel
        plot(trial(I(fnum)).epochTimeVec,trial(I(fnum)).emgStim{c});hold on;
        plot(trial(I(fnum)).epochTimeVec,trial(I(fnum)).meanTrace{c},'k','LineWidth',1.2);
        vline([0    trial(I(fnum)).stimLength/fs],'r-');
        tmptit = title([trial(I(fnum)).Name ' - ' num2str(stimAmp(I(fnum))/1000) 'mA, ' ...
            num2str(pulseWidth(I(fnum))) 'us, ' num2str(stimFrequency(I(fnum))) 'Hz, ' num2str(stimDuration(I(fnum))) 's']); 
        set(hFs(fnum),'FontSize',7);
        set(tmptit, 'Interpreter', 'none','FontSize',5);
        hFs(fnum).XLim = [trial(I(fnum)).epochTimeVec(1) trial(I(fnum)).epochTimeVec(end)];
        
        subCount = subCount + 1;
    end
    linkaxes([hFs],'xy');
    suplabel('Amplitude(uV)','y');
    suplabel('Time(sec)','x');
%     set(hF(c),'Position', [1343 49 570 947]);
    
    disp('Saving to Set Folder')
    savefig(hF(c),[setPath Muscle_ID{c} '_EpochedStimTraces_' setName '_' setDescrpt]);
    saveas(hF(c),[setPath Muscle_ID{c} '_EpochedStimTraces_' setName '_' setDescrpt '.png']);
end

close all