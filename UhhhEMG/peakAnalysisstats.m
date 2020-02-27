load('sitting_peaks.mat')
%% Characterize Muscle Response Info from peaks.mat

Muscle_ID = {'Right VM', 'Right RF', 'Right VL', 'Right BF', 'Right ST', ...
    'Right TA', 'Right SO', 'Right LG', 'Left VM', 'Left RF', 'Left VL',...
    'Left BF', 'Left ST', 'Left Ham', 'Left SO', 'Left LG'};

path_Save = 'R:\users\dsarma\';
path_UH3 = 'R:\data_raw\human\uh3_stim\';
subjectName = 'LSP02b';
reportPath = ['R:\users\dsarma\' subjectName '\Figures\'];


for m = 1:length(Muscle_ID)
    mResponse.ID(m,:) = Muscle_ID(m);
end


%Max Peak is first item in pks command.

%m is muscle Number
%e is electrode set
%f is trial
for m = 1:length(Muscle_ID)
    disp(['Processing: ' Muscle_ID(m)]);
    for e = 1:length(peakInfo)
        if e == 9
            continue;
        end
        for f = 1:length(peakInfo(e).muscle(m).sitting)
            [mPk{e,f}, I_Pk] = max(peakInfo(e).muscle(m).sitting(f).pks); %max peak value
            tPk{e,f} = peakInfo(e).muscle(m).sitting(f).locs(I_Pk); %onset of max peak
            wPk{e,f} = peakInfo(e).muscle(m).sitting(f).widths(I_Pk);
            width{e,f} = sum([peakInfo(e).muscle(m).sitting(f).widths]);
        end
        stAmps{e,:} = [peakInfo(e).muscle(m).sitting.stimAmp]; %stimAmp of eStim
        pkMns{e,:} = [peakInfo(e).muscle(m).sitting.mean]; %mean/RMS activation
        pkVars{e,:} = [peakInfo(e).muscle(m).sitting.var]; %var of peaks
        pkStd{e,:} = [peakInfo(e).muscle(m).sitting.std]; %std of peaks
    end
% 
% for m = 1%:length(Muscle_ID)
%     for e = 1:length(peakInfo)
%           stAmps{e,:} = [peakInfo(e).muscle(m).sitting.stimAmp]; %stimAmp of eStim
%     end
% end

        allStimAmps = 100:50:6000;
        emptyVessel = zeros(length(allStimAmps),size(stAmps,1))';
        allPeaks = emptyVessel;
        allOnsets = emptyVessel;
        allPeakWidths = emptyVessel;
        allPeakDurations = emptyVessel;
        allMeans = emptyVessel;
        allMeanVars = emptyVessel;
        allMeanStds = emptyVessel;


        for n = 1:length(allStimAmps)
            for f = 1:size(stAmps,1)
                Idx = find([stAmps{f,:}] == allStimAmps(n));
                if isempty(Idx) || isempty(tPk{f,Idx(1)})
                   allPeaks(f,n) = NaN;
                   allOnsets(f,n) = NaN;
                   allPeakWidths(f,n) = NaN;
                   allPeakDurations(f,n) = NaN;
                   allMeans(f,n) = NaN;
                   allMeanVars(f,n) = NaN;
                   allMeanStds(f,n) = NaN;
% %                    f
               else
                   allPeaks(f,n) = mPk{f,Idx(1)};
                   allOnsets(f,n) = tPk{f,Idx(1)};
                   allPeakWidths(f,n) = wPk{f,Idx(1)};
                   allPeakDurations(f,n) = width{f,Idx(1)};
                   allMeans(f,n) = pkMns{f}(Idx(1));
                   allMeanVars(f,n) = pkVars{f}(Idx(1));
                   allMeanStds(f,n) = pkStd{f}(Idx(1));
% %                    f
                end
            end
        end
            
        
        meanMaxPeaks = nanmean(allPeaks); stdErrMaxPeaks = nanstd(allPeaks)/size(allPeaks,1);
        meanOnsets = nanmean(allOnsets); stdErrOnsets = nanstd(allOnsets)/size(allOnsets,1);
        meanPkWidths = nanmean(allPeakWidths); stdErrPkWidths = nanstd(allPeakWidths)/size(allPeakWidths,1);
        meanDuration = nanmean(allPeakDurations); stdErrDuration = nanstd(allPeakDurations)/size(allPeakDurations,1);
        allStimAmps = allStimAmps/1000;
        %% Plotting  Peak, Onset, Width, Duration
        hFig(m) =  figure;maximize;

        % Plot Max Responses
        hs(m,1) = subplot (411);
        b(1) = bar(allStimAmps,meanMaxPeaks');
        hold on;
        er(1) = errorbar(allStimAmps,meanMaxPeaks,stdErrMaxPeaks,'k', 'LineWidth', 2,'Capsize', 10);    
        % er(1).Color = [0 0 0];                            
        er(1).LineStyle = 'none';
        
        %%getting indices where y is valid (not NaN)
        idxValid = ~isnan(meanMaxPeaks);
        %%fitting
%         poly = polyfit(allStimAmps(idxValid),meanMaxPeaks(idxValid),3);
        [poly, gof] = fit(allStimAmps(idxValid)',meanMaxPeaks(idxValid)','poly3','Normalize','on','Robust','Bisquare');
        plot(poly,'r'); 
        legend('Mean Response','stderr',['fit: ' num2str(gof.rsquare)]);
        
        hold off;
        title('Peak Response by Amplitude');
        ylabel('Response Voltage (mV)');
        xlabel('Stimulation Amplitudes (mA)');

        % Plot Onset Times
        hs(m,2) = subplot (412);
        b(2) = bar(allStimAmps,meanOnsets');

        hold on;
        er(2) = errorbar(allStimAmps,meanOnsets,stdErrOnsets,'k', 'LineWidth', 2,'Capsize', 10);    
        % er(2).Color = [0 0 0];                            
        er(2).LineStyle = 'none';
        
        %%getting indices where y is valid (not NaN)
        idxValid = ~isnan(meanOnsets);
        %%fitting
        [poly, gof] = fit(allStimAmps(idxValid)',meanOnsets(idxValid)','poly3','Normalize','on','Robust','Bisquare');
        plot(poly,'r'); 
        legend('Mean Response','stderr',['fit: ' num2str(gof.rsquare)]);
        
        hold off;
        title('Response Onset by Amplitude');
        ylabel('Response Onset (sec)');
        xlabel('Stimulation Amplitudes (mA)');
        
        % Plot Width of Peaks
        hs(m,3) = subplot (413);
        b(3) = bar(allStimAmps,meanPkWidths');

        hold on;
        er(3) = errorbar(allStimAmps,meanPkWidths,stdErrPkWidths,'k', 'LineWidth', 2,'Capsize', 10);    
        % er(3).Color = [0 0 0];                            
        er(3).LineStyle = 'none';
        
        %%getting indices where y is valid (not NaN)
        idxValid = ~isnan(meanPkWidths);
        %%fitting
        %%fitting
        [poly, gof] = fit(allStimAmps(idxValid)',meanPkWidths(idxValid)','poly3','Normalize','on','Robust','Bisquare');
        plot(poly,'r'); 
        legend('Mean Response','stderr',['fit: ' num2str(gof.rsquare)]);
        
        hold off;
        title('Width of Peak by Amplitude');
        ylabel('Peak Duration (sec)');
        xlabel('Stimulation Amplitudes (mA)');

        % Plot Response Duration
        hs(m,4) = subplot (414);
        b(4) = bar(allStimAmps,meanDuration');

        hold on;
        er(4) = errorbar(allStimAmps,meanDuration,stdErrDuration,'k', 'LineWidth', 2,'Capsize', 10);    
        % er(4).Color = [0 0 0];                            
        er(4).LineStyle = 'none';
        
        %%getting indices where y is valid (not NaN)
        idxValid = ~isnan(meanDuration);
        %%fitting
        [poly, gof] = fit(allStimAmps(idxValid)',meanDuration(idxValid)','poly3','Normalize','on','Robust','Bisquare');
        plot(poly,'r'); 
        legend('Mean Response','stderr',['fit: ' num2str(gof.rsquare)]);
        
        hold off;
        tmptit = suptitle([Muscle_ID(m) 'Seated Stim, Response Summary']);
        set(tmptit,'Fontsize', 20)
        title('Response Duration by Amplitude');
        ylabel('Response Duration (sec)');
        xlabel('Stimulation Amplitudes (mA)');


        
        linkaxes([hs],'x');
        hs(m,1).XLim = [0 1.1*max(allStimAmps)];
        
    clearvars -except peakInfo Muscle_ID mResponse m hFig hs reportPath
end

a = num2cell(0:500:6000);
a = cellfun(@(x) x/1000,a,'un',0);
for i = 1:4
% %     linkaxes([hs(:,i)],'y');
% %     limsy= hs(m,i).YLim; hs(m,i).YLim = [0 limsy(2)];
    for c = 1:length(Muscle_ID)
        limsy= hs(c,i).YLim; hs(c,i).YLim = [0 inf];
        hs(c,i).FontSize = 16;
        hs(c,i).XTick = cellfun(@(x) x,a);
    end
    if i==1
        continue
    else
        linkaxes([hs(:,i)],'y');
        limsy= hs(c,i).YLim; hs(c,i).YLim = [0 limsy(2)];
    end
end

for  m = 1:length(Muscle_ID)
        disp(['Saving Muscle: ' Muscle_ID(m)]);
        saveas(hFig(m),[reportPath Muscle_ID{m} '_Seated_ResponseSummary_unlinked.png']);
        savefig(hFig(m),[reportPath Muscle_ID{m} '_Seated_ResponseSummary_unlinked']);
end
close all



