 
for hIdx = 1:5
    stimAmps = [hash(hIdx).standing.stimAmp];
    [sortedStim,I] = sort(stimAmps);
    epochTimeVec = hash(hIdx).standing.epochTimeVec;
    iter = 1;

    disp(['Plotting Case: ' num2str(hIdx)]);

    for c = hash(hIdx).Indices.muscles
        disp('Plotting Stim Averages');

        hF(c) = figure; %('visible', 'off');
        tmpsup = suptitle([Muscle_ID{c} ' - Day ' num2str(hash(hIdx).day) ' - Elec ' num2str(hash(hIdx).elecNum) ' - Standing -  Stim Traces']);
        set(tmpsup, 'Interpreter', 'none') 
        maximize;

        subCount = 1; 
        for fnum = 1:iter:length(hash(hIdx).Indices.standing)
            if length(hash(hIdx).Indices.standing) > 20
                hFs(fnum) = subplot(9,5,subCount);
            else
                hFs(fnum) = subplot(8,2,subCount);
            end
    %         hFs(fnum) = subplot(10,1,subCount);
            %Plot /Channel
            plot(hash(hIdx).standing(I(fnum)).epochTimeVec,hash(hIdx).standing(I(fnum)).emgStim{c});hold on;
            plot(hash(hIdx).standing(I(fnum)).epochTimeVec,hash(hIdx).standing(I(fnum)).meanTrace{c},'k','LineWidth',1.2);
            vline([0    hash(hIdx).standing(I(fnum)).stimLength/fs],'r-');
            tmptit = title([num2str(stimAmps(I(fnum))/1000) 'mA, ' num2str(pulseWidth(I(fnum))) 'us, 1Hz, 10s']); 
            set(hFs(fnum),'FontSize',7);
            set(tmptit, 'Interpreter', 'none','FontSize',5);
            hFs(fnum).XLim = [hash(hIdx).standing(I(fnum)).epochTimeVec(hIdx) hash(hIdx).standing(I(fnum)).epochTimeVec(end)];

            subCount = subCount + 1;
        end
        linkaxes([hFs],'xy');
        suplabel('Amplitude(uV)','y');
        suplabel('Time(sec)','x');
    %     set(hF(c),'Position', [1343 49 570 947]);

        disp('Saving to Set Folder')
        savefig(hF(c),[reportPath Muscle_ID{c} '_EpochedStimTraces_Day' num2str(hash(hIdx).day) '_Elec' num2str(hash(hIdx).elecNum) '_Standing']);
        saveas(hF(c),[reportPath Muscle_ID{c} '_EpochedStimTraces_Day' num2str(hash(hIdx).day) '_Elec' num2str(hash(hIdx).elecNum) '_Standing.png']);
    end
    pause (0.1);
    close all
end

for hIdx = 1:5
    stimAmps = [hash(hIdx).sitting.stimAmp];
    [sortedStim,I] = sort(stimAmps);
    epochTimeVec = hash(hIdx).sitting.epochTimeVec;
    iter = 1;

    disp(['Plotting Case: ' num2str(hIdx)]);

    for c = hash(hIdx).Indices.muscles
        disp('Plotting Stim Averages');

        hF(c) = figure; %('visible', 'off');
        tmpsup = suptitle([Muscle_ID{c} ' - Day ' num2str(hash(hIdx).day) ' - Elec ' num2str(hash(hIdx).elecNum) ' - Sitting -  Stim Traces']);
        set(tmpsup, 'Interpreter', 'none') 
        maximize;

        subCount = 1; 
        for fnum = 1:iter:length(hash(hIdx).Indices.sitting)
            if length(hash(hIdx).Indices.sitting) > 20
                hFs(fnum) = subplot(9,5,subCount);
            else
                hFs(fnum) = subplot(8,2,subCount);
            end
    %         hFs(fnum) = subplot(10,1,subCount);
            %Plot /Channel
            plot(hash(hIdx).sitting(I(fnum)).epochTimeVec,hash(hIdx).sitting(I(fnum)).emgStim{c});hold on;
            plot(hash(hIdx).sitting(I(fnum)).epochTimeVec,hash(hIdx).sitting(I(fnum)).meanTrace{c},'k','LineWidth',1.2);
            vline([0    hash(hIdx).sitting(I(fnum)).stimLength/fs],'r-');
            tmptit = title([num2str(stimAmps(I(fnum))/1000) 'mA, ' num2str(pulseWidth(I(fnum))) 'us, 1Hz, 10s']); 
            set(hFs(fnum),'FontSize',7);
            set(tmptit, 'Interpreter', 'none','FontSize',5);
            hFs(fnum).XLim = [hash(hIdx).sitting(I(fnum)).epochTimeVec(hIdx) hash(hIdx).sitting(I(fnum)).epochTimeVec(end)];

            subCount = subCount + 1;
        end
        linkaxes([hFs],'xy');
        suplabel('Amplitude(uV)','y');
        suplabel('Time(sec)','x');
    %     set(hF(c),'Position', [1343 49 570 947]);

        disp('Saving to Set Folder')
        savefig(hF(c),[reportPath Muscle_ID{c} '_EpochedStimTraces_Day' num2str(hash(hIdx).day) '_Elec' num2str(hash(hIdx).elecNum) '_Sitting']);
        saveas(hF(c),[reportPath Muscle_ID{c} '_EpochedStimTraces_Day' num2str(hash(hIdx).day) '_Elec' num2str(hash(hIdx).elecNum) '_Sitting.png']);
    end
    pause (0.1);
    close all
end


