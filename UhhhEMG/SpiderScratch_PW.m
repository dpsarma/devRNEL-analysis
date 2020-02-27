%% Making Spider Plots for PW Mod
% Things we need, MuscleID, Response Value (P2P or RMS), baseline Values, and Stim Amps, and electrode Idx


 %Create a vector that is electrode by muscle: Rbinary = [M1:M16] where each Rbinary element vector is an electrode and each value % is the threshold at which the muscle is activated 
 % First identify for every electrode, the responses across level and then sort.  Can make a separate radar plot that is for each stimulation level with the relative amplitude of the response. 

% Sort Trials By Electrodes 1st

elecs = unique(spinalElec);
RBinary = zeros(length(elecs),length(Muscle_ID));

for ii = 1:length(elecs)
    eIdx = find(spinalElec == elecs(ii));
    
    % Sort Stim Amps by Electrode
    [sortedPulse,I] = sort(pulseWidth(eIdx));
    for mm = 1:length(Muscle_ID)
            for ll =  1:length(eIdx)
                XS(ll,:) = max(trial(eIdx(I(ll))).baseline(mm));
                YS(ll,:) = max(trial(eIdx(I(ll))).responseRMS(mm));
            end
			mBase = mean(XS);
            stdBase = std(XS);
            YS = YS - mBase;
            for ss = 1:length(XS)
                if YS(ss) > 2*stdBase
                	Rstim(ss) = 1;
               	else
               		Rstim(ss) = 0;
               	end
            end
              tmp =  find(Rstim);
              if isempty(tmp) 
                  RBinary(ii,mm) = nan;
              else
                  RBinary(ii,mm) = sortedPulse(tmp(1));
              end
             ResponsetoStim{:,mm} = Rstim;
             clear XS YS Rstim
    end
   
   if ~isempty(anodeElec)
       leglabs{ii} = ['E' num2str(elecs(ii)) '::' num2str(anodeElec(eIdx(1)))]; 
   else
       leglabs{ii} = ['E' num2str(elecs(ii))];
   end
end
%Plot Individual Electrodes
spidermin(1:8) = 100; spdrmax(1:8) = 1000; spdriter = 5; spdrfill = 'off';
for ii = 1:length(elecs)
   hS0 = figure; maximize;
    subplot(121); spider_plot([RBinary(ii,9:16);],...
        'AxesLabels', Muscle_ID(9:16),'AxesInterval', 2,...
        'FillOption', spdrfill,'FillTransparency', 0.1, 'AxesLimits',...
        [spidermin; spdrmax],...
        'AxesInterval', spdriter, 'Direction', 'clockwise');
    title('Residual Limb (left)');
    % Legend properties
    legend(leglabs{ii}, 'Location', 'southoutside', 'Orientation', 'horizontal');
    subplot(122); spider_plot([RBinary(ii,1:8);],...
        'AxesLabels', Muscle_ID(1:8),'AxesInterval', 2,...
        'FillOption', spdrfill,'FillTransparency', 0.1, 'AxesLimits',...
        [spidermin; spdrmax],...
        'AxesInterval', spdriter);
    title('Intact Limb (right)');
    legend(leglabs{ii}, 'Location', 'southoutside', 'Orientation', 'horizontal');

    tmptit = sgtitle(['Thresholds for ' num2str(elecs(ii)) ' - ' setDescrpt]);
    set(tmptit, 'Interpreter', 'none');

    saveas(hS0,[reportPath '\' strrep(datestr(now),':','_')  '_PWThresholdPlot' '_' num2str(elecs(ii)) '.png']);
    savefig(hS0,[reportPath '\' strrep(datestr(now),':','_')  '_PWThresholdPlot' '_' num2str(elecs(ii))]);
    close all

end

% % %% Plot All Electrodes Together
% % hS1 = figure; maximize;
% % spdrfill = 'on';
% %     subplot(121); spider_plot([RBinary(:,9:16);],...
% %         'AxesLabels', Muscle_ID(9:16),'AxesInterval', 2,...
% %         'FillOption', spdrfill,'FillTransparency', 0.1, 'AxesLimits',...
% %         [spidermin; spdrmax],...
% %         'AxesInterval', spdriter, 'Direction', 'clockwise');
% %     title('Residual Limb (left)');
% %     % Legend properties
% %     legend(leglabs, 'Location', 'southoutside', 'Orientation', 'horizontal');
% %     subplot(122); spider_plot([RBinary(:,1:8);],...
% %         'AxesLabels', Muscle_ID(1:8),'AxesInterval', 2,...
% %         'FillOption', spdrfill,'FillTransparency', 0.1, 'AxesLimits',...
% %         [spidermin; spdrmax],...
% %         'AxesInterval', spdriter);
% %     title('Intact Limb (right)');
% %     legend(leglabs, 'Location', 'southoutside', 'Orientation', 'horizontal');
% % 
% %     tmptit = sgtitle(['Thresholds for ' setDescrpt]);
% %     set(tmptit, 'Interpreter', 'none');
% % 
% % saveas(hS1,[reportPath '\' strrep(datestr(now),':','_')  '_ThresholdPlot' '_' setName '.png']);
% % savefig(hS1,[reportPath '\' strrep(datestr(now),':','_')  '_ThresholdPlot' '_' setName]);
% % 
% % close all


%% Plot Responses across all trials in a Set
switch multipolar
    case 'yes'
        hS2 = figure; maximize;
        for mm = 1:length(Muscle_ID)
            for ll =  1:length(trial)
                X(ll,:) = trial(ll).baseline(mm);
                Y(ll,:) = trial(ll).responseRMS(mm);
                Z(ll,:) = trial(ll).responseP2P(mm);
                label{ll,:} =  [num2str(spinalElec(ll)) '::' num2str(anodeElec(ll)) ' - ' num2str(stimAmp(ll)) 'mA'];
            end
            subplot(2,8,mm);
            plot(X); hold on; plot(Y); plot(Z);
            mBase = mean(X);
            stdBase = std(X);
            hline([mBase (mBase+2*stdBase) (mBase-2*stdBase)]);
            xlabel('Trial'); ylabel('EMG in uV)');
            title(Muscle_ID(mm));
            xticks(1:length(trial));
            xtickangle(90);
            xticklabels(label);
            legend('Baseline', 'RMS','P2P');
        end
    otherwise
        hS2 = figure; maximize;
        for mm = 1:length(Muscle_ID)
            for ll =  1:length(trial)
                X(ll,:) = trial(ll).baseline(mm);
                Y(ll,:) = trial(ll).responseRMS(mm);
                Z(ll,:) = trial(ll).responseP2P(mm);
                label{ll,:} =  [num2str(spinalElec(ll)) ' - ' num2str(stimAmp(ll)) 'mA'];
            end
            subplot(2,8,mm);
            plot(X); hold on; plot(Y); plot(Z);
            mBase = mean(X);
            stdBase = std(X);
            hline([mBase (mBase+2*stdBase) (mBase-2*stdBase)]);
            xlabel('Trial'); ylabel('EMG in uV)');
            title(Muscle_ID(mm));
            xticks(1:length(trial));
            xtickangle(90);
            xticklabels(label);
            legend('Baseline', 'RMS','P2P');
        end
end

saveas(hS2,[reportPath '\' strrep(datestr(now),':','_')  '_AllTrialsResponse' '_' setName '.png']);
savefig(hS2,[reportPath '\' strrep(datestr(now),':','_')  '_AllTrialsResponse' '_' setName]);

close all;




% % %% Plotting Stim Averages
elecs = unique(spinalElec);
for ii = 1:length(elecs)
    eIdx = find(spinalElec == elecs(ii));
%     [sortedStim,I] = sort(stimAmp(eIdx));
    [sortedPulse,I2] = sort(pulseWidth(eIdx));
%     [sortedFreq,I3] = sort(stimFrequency(eIdx));
    
    %% Plotting Recruitment Curves
    response = [trial(:).responseRMS]';
% %     response = [trial(:).responseP2P]';

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
            
            for ll =  1:length(trial)
                X(ll,:) = trial(ll).baseline(m);
            end
            mBase = mean(X);
            yyP = smooth(sortedPulse,response(I2,m)-mBase,'sgolay');
            plot(sortedPulse,yyP,'Color',C(m,:),'LineWidth',3);
           
            hold on;
        end
            ylabel('(post-Stim) basesub-Mean RMS  (uV)'); %'(post-Stim) Mean RMS  (uV)'
            xlabel('Stimulation Amplitude (mA)');
            title(['E' num2str(elecs(ii)) ' - By Amplitude']);


        if pI == 1
            legend([Muscle_ID(1:mCount)],'Location','northwest', 'Orientation','vertical');
        elseif pI == 2
            legend([Muscle_ID((mCount)+1:length(Muscle_ID))],'Location','northwest', 'Orientation','vertical');
        end
    %     
    end
    if ~isempty(anodeElec)
        tmptit = sgtitle([setName ' - ' setDescrpt ' - e' num2str(elecs(ii)) '::' num2str(anodeElec(eIdx(1))) ' -  PW Recruitment Curves (RMS)']);
    else
        tmptit = sgtitle([setName ' - ' setDescrpt ' - e' num2str(elecs(ii)) ' -  PW Recruitment Curves (RMS)']);
    end
    set(tmptit, 'Interpreter', 'none');
    % % Pix_SS = get(0,'screensize');
    hF3.Position = [680 49 572 947];
%     linkaxes([hF3s],'xy');

    if ~isempty(anodeElec)
        saveas(hF3,[reportPath '\' strrep(datestr(now),':','_')  '_bsrPWRecruitmentCurvesRMS' '_' setName '_ E' num2str(elecs(ii)) '-' num2str(anodeElec(eIdx(1))) '.png']);
        savefig(hF3,[reportPath '\' strrep(datestr(now),':','_') '_bsrPWRecruitmentCurvesRMS_' setName '_ E' num2str(elecs(ii)) '-' num2str(anodeElec(eIdx(1))) ]);
    else
        saveas(hF3,[reportPath '\' strrep(datestr(now),':','_')  '_bsrPWRecruitmentCurvesRMS' '_' setName '_ E' num2str(elecs(ii)) '.png']);
        savefig(hF3,[reportPath '\' strrep(datestr(now),':','_') '_bsrPWRecruitmentCurvesRMS_' setName '_ E' num2str(elecs(ii)) ]);
    end
end
close all

 


  

clear XS
clear Rstim
clear YS
clear XS
% clear RBinary
% clear ResponsetoStim