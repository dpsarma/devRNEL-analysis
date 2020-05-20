%%%% UH3 EMG - Recruitment Curves - Stability - Testing Version B %%
% %     %Proto code for looping
% % % % for p = 1:length(pwidths)
% % % %     check_p = pwidths(p);
% % % %     disp(['Pulse is ' num2str(check_p)]);
% % % %     for f = 1:length(freqs)
% % % %         check_f = freqs(f);
% % % %         disp(['Freq is ' num2str(check_p)]);
% % % %         eIdx = find(([trial.stimFrequency] == check_f) & ([trial.pulseWidth] == check_p));
% % % %         for t = 1:length(eIdx)
% % % %             n = eIdx(t);
% % % %             for i = 1:length(msubset)
% % % %                 m = msubset(i);
% % % %             end
% % % %         end    
% % % %     end
% % % % end
for r = 1:4

    switch r
        case 1
            % elecs 16: Day 6
            load('D:\FigRescources\UH3\LSP05\rehash\Stability_E16_day6.mat');
            multipolar = 'no';
            anodeElec = [];
            disp('elecs 16: Day 6');
        case 2
            % elecs 16: Day 8
            load('D:\FigRescources\UH3\LSP05\rehash\Stability_E16_day8.mat');
            multipolar = 'no';
            anodeElec = [];
            disp('elecs 16: Day 8');
        case 3
            % elecs 16: Day 12
            load('D:\FigRescources\UH3\LSP05\rehash\Stability_E16_day12.mat');
            multipolar = 'no';
            anodeElec = [];
            disp('elecs 16: Day 12');
        case 4
            % elecs 16: Day 15
            load('D:\FigRescources\UH3\LSP05\rehash\Stability_E16_day15.mat');
            disp('elecs 16: Day 15');
        otherwise
            % multipolar: day9_set3
            disp('Pick a fileset');
    end

% Basic Recruitment Curve Stuff for Frequency and Set on all Elecs
    freqs = unique([trial.stimFrequency]);
    pwidths = unique([trial.pulseWidth]);
    elecs = unique([trial.spinalElec]);



% Swapped RF & VL so 3<->2 & 11<->10
msubset = [9:16]; %or 1:8 for right/intact
msubset([2 3]) = msubset([3 2]);

% Make RC vectors for each subset.
check_p = 0.5;
check_f = 1;
artifactBuffer = 0.005;
disp(['Processing: ' num2str(check_p) 'ms ' num2str(check_f) 'Hz']);

for e = 1:length(elecs)
    response(r).elec = elecs(e);
    response(r).freq = check_f;
    response(r).pulse = check_p;
    response(r).setDescrpt = setDescrpt;
    
    eIdx = find(([trial.stimFrequency] == check_f) & ([trial.pulseWidth] == check_p) & ([trial.spinalElec] == elecs(e)));
    [sortedStim,I] = sort([trial.stimAmp]);

    response(r).StimAmps = sortedStim;
    
    for i = 1:length(msubset)
    m = msubset(i);
    disp(['Processing ' Muscle_ID(m)]);
    response(r).muscle(m).Muscle_ID = Muscle_ID(m);
    
        for t = 1:length(eIdx)
            n = eIdx(I(t));
            stimBegins = trial(n).stims;
            stimEnd = nSampsPre + trial(n).stimLength + artifactBuffer*fs;  
            

            response(r).muscle(m).responseRMS(t) = trial(n).responseRMS(m);
            response(r).muscle(m).p2pResponse(t) = peak2peak(trial(n).meanTrace{m}(stimEnd:end));
            [response(r).muscle(m).peakResponse(t), mi] = max(abs(trial(n).meanTrace{m}(stimEnd:end)));
            response(r).muscle(m).onset(t) = trial(n).epochTimeVec(mi+stimEnd-1);
            response(r).muscle(m).meanbase(t) = mean(trial(n).baseline(m,:));
            response(r).muscle(m).rmsbase(t) = rms(trial(n).baseline(m,:));

            disp([num2str(t) ' of ' num2str(length(eIdx))]);
        end

    end

    % Get Threshold
    for i = 1:length(msubset)
        m = msubset(i);
        response(r).muscle(m).threshIdx = findchangepts(cumsum(response(r).muscle(m).p2pResponse), 'Statistic','rms');
        response(r).muscle(m).threshold = (sortedStim(response(r).muscle(m).threshIdx));
%         
    end
end

end


%Save Response Variable for later


%% Plot RC for Peak2Peak, RMS, and Max Peak (Abs)
setPath = 'D:\FigRescources\UH3\LSP05\rehash\Stability\';
% Using the saved response variable created above
for r = 1:length(response)
    hF(r) = figure; maximize;
    for i = 1:length(msubset)
        m = msubset(i);
        h(i) = subplot(2,length(msubset)/2,i);
        plot(response(r).StimAmps/1000, response(r).muscle(m).peakResponse);
        hold on; plot(response(r).StimAmps/1000, response(r).muscle(m).p2pResponse);
        hold on; plot(response(r).StimAmps/1000, response(r).muscle(m).responseRMS);
        title(response(r).muscle(m).Muscle_ID);
        ylabel('EMG (uV)');
        xlabel('StimAmp (mA)');

        yyaxis right
        plot(response(r).StimAmps/1000, response(r).muscle(m).onset*1000,'-.');
        ylabel('Onset (ms)');

        legend('Abs Peak', 'RMS', 'P2P', 'Onset','Location','north');
    end
    % linkaxes([h(:)],'xy');
    tit = sgtitle([response(r).setDescrpt ' - RCs plotted 3 ways (' num2str(check_p) 'ms ' num2str(check_f) 'Hz)']);
    text = tit.String
    
    saveas(hF(r),[setPath tit.String '.png'])
end


% % %% Plot Baselines
% % hF(2) = figure; maximize;
% % for i = 1:length(msubset)
% %     m = msubset(i);
% %     h2(i) = subplot(2,length(msubset)/2,i);
% %     plot(response(r).StimAmps/1000, response(r).muscle(m).meanbase);
% %     hold on; plot(response(r).StimAmps/1000,  response(r).muscle(m).rmsbase);
% %     hold on; plot(response(r).StimAmps/1000, response(r).muscle(m).p2pResponse);
% %     title(response(r).muscle(m).Muscle_ID);
% %     ylabel('EMG (uV)');
% %     xlabel('StimAmp (mA)');
% %     
% %     legend('mean base', 'RMS base', 'P2P');
% % end
% % % linkaxes([h(:)],'xy');
% % tit{2} = sgtitle(['Baselines+P2P (' num2str(check_p) 'ms ' num2str(check_f) 'Hz)']);
% % text{2} = tit{2}.String


%% Plot CumSum vs P2P
for r = 1:length(response)
    hF(r) = figure; maximize;
    for i = 1:length(msubset)
        m = msubset(i);
        h2(i) = subplot(2,length(msubset)/2,i);
        plot(response(r).StimAmps/1000, response(r).muscle(m).p2pResponse);
        hold on; plot(response(r).StimAmps/1000, cumsum(response(r).muscle(m).p2pResponse));
        title(response(r).muscle(m).Muscle_ID);
        ylabel('EMG (uV)');
        xlabel('StimAmp (mA)');

        legend('P2P','cumsum');
    end
    % linkaxes([h(:)],'xy');
    tit = sgtitle([response(r).setDescrpt ' - PS2 vs Cumsum (' num2str(check_p) 'ms ' num2str(check_f) 'Hz)']);
    text = tit.String
    
    saveas(hF(r),[setPath tit.String '.png'])
end

% % %% Plot Change Points
% % for r = 1:length(response)
% %     for i = 1:length(msubset)
% %         m = msubset(i);
% %         hf4(i) = figure; maximize;
% %         findchangepts(cumsum(response(r).muscle(m).p2pResponse), 'Statistic','rms');
% % %         htit{i} = title(['E' num2str(elecs(r)) ' - ' char(response(r).muscle(m).Muscle_ID)]);
% % %         ylabel('EMG (uV)');
% % %         xlabel('StimAmp (mA)');
% %         
% % %         saveas(hf4(f),[setPath 'E' num2str(elecs(r)) '-Inflection.png'])
% %     end
% %     
% % end
% % 
% % % % close all;

%% Plot Comparison RC


msubset = [10 12 14 16];
elecs = unique([response(:).elec]);%[1 8 16 18 23 29];
C = linspecer(length(response));

for i = 1:length(msubset)
    m = msubset(i);
    h(i) = subplot(1,length(msubset),i);
    for r = 1:length(response)
% %         e = find([response.elec] == elecs(r));
        tmp(m,r) = max(response(r).muscle(m).p2pResponse);
    end
    maxR(m) = max(tmp(m,:));
end

hF = figure; maximize;
for i = 1:length(msubset)
    m = msubset(i);
    h(i) = subplot(1,length(msubset),i);
    for r = 1:length(response)
        e = r;%find([response.elec] == elecs(r));
        plot(response(e).StimAmps/1000, response(e).muscle(m).p2pResponse,'Color',C(r,:), 'LineWidth',2);
        hold on;
    end
    title(response(e).muscle(m).Muscle_ID);
    ylabel('EMG (uV)');
    xlabel('StimAmp (mA)');
        
    tmpleg = legend(response(:).setDescrpt,'Interpreter', 'none', 'location', 'west');
end

% linkaxes([h(:)],'xy');
tit = sgtitle(['Recruitment Curves Across Days']);
saveas(hF,['D:\FigRescources\UH3\LSP05\rehash\Stability\E' num2str(elecs(1)) '_P2PacrossDays.png'])
saveas(hF,['D:\FigRescources\UH3\LSP05\rehash\Stability\E' num2str(elecs(1)) '_P2PacrossDays.svg'])

close all;

%% Plot Comparison RC Stacked


msubset = [10 12 16 14];
% % etmp= [response.elec];
% % elecs = etmp(1:4:end);
elecs = unique([response(:).elec]);%[1 8 16 18 23 29];
C = linspecer(length(response));

for i = 1:length(msubset)
    m = msubset(i);
    for r = 1:length(elecs)
        e = r;%find([response.elec] == elecs(r));
        tmp(m,r) = max(response(e).muscle(m).p2pResponse);
    end
    maxR(m) = max(tmp(m,:));
end

hF = figure;% maximize;
for i = 1:length(msubset)
    m = msubset(i);
    
    for r = 1:length(response)
        h(i,r) = subplot(length(response),length(msubset),length(msubset)*(r-1)+i);
        e = r;%find([response.elec] == elecs(r));
        p = area(response(e).StimAmps/1000, response(e).muscle(m).p2pResponse);
        p.FaceColor = C(r,:);
        hold on;
        if i==1
            ylabel({response(e).setDescrpt,'EMG (uV)'},'Interpreter','none');
        end
        if r==length(elecs)
            xlabel('StimAmp (mA)');
        end
        if r==1
            title(response(e).muscle(m).Muscle_ID);
        end
%         legend(['E' num2str(elecs(r))],'Location','northwest');
        box off
% %         h(i,r).Color = 'none';
    end
    
     
    
end
for i = 1:length(msubset)
linkaxes([h(i,:)],'xy');
end
% % linkaxes([h(:,1:end-2)],'xy');
tit = sgtitle(['Recruitment Curves across days']);
set(0,'defaultAxesFontSize',12)
saveas(hF,['D:\FigRescources\UH3\LSP05\rehash\Stability\Stacked_RCs_E' num2str(elecs(1)) '_P2PacrossDays.png'])
saveas(hF,['D:\FigRescources\UH3\LSP05\rehash\Stability\Stacked_RCs_E' num2str(elecs(1)) '_P2PacrossDays.svg'])

close all;



%% Plot Threshold Comparisons
musclesL = {'VM', 'VL', 'RF', 'BF', 'ST', 'TA', 'MG', 'LG',};
div = 1000;

msubset = [9 10 11 12 13 14 15 16];
clear threshold
for r = 1:length(response)
    threshold(r,:) = [response(r).muscle(msubset).threshold]/div;
end

hS1 = figure; maximize;
spdrfill = 'on';
spdrmax(1:8) =((6000/div)); 
spdrmin(1:8) = ((300/div)); 
spdriter = 6; 

C = linspecer(length(response));

    spider_plot(threshold,...
        'AxesLabels', musclesL(msubset-8),'AxesInterval', 2,...
        'FillOption', spdrfill,'FillTransparency', 0.1, 'AxesLimits',...
        [spdrmin(1:length(msubset)); spdrmax(1:length(msubset))],...
        'AxesInterval', spdriter, 'Direction', 'clockwise',...
        'Color', C);
%     title('Residual Limb (left)');
    % Legend properties
    legend({response(:).setDescrpt}, 'Location', 'westoutside', 'Orientation', 'vertical','Interpreter','none');
     

    tmptit = sgtitle(['Thresholds for ' num2str(elecs(1)) '-' num2str(elecs(end)) ' (1Hz, 500us)']);
    set(tmptit, 'Interpreter', 'none');

% tit2 = sgtitle('Threshold by Muscle for Electrodes 9 & 16');
saveas(hS1,['D:\FigRescources\UH3\LSP05\rehash\Stability\Threshold_e' num2str(elecs(1)) '.png']);
saveas(hS1,['D:\FigRescources\UH3\LSP05\rehash\Stability\Threshold_e' num2str(elecs(1)) '.svg']);

%% Activation at Threshold vs Activation Ratio
sensoryThresh = [2000];

for r = 1:length(response)
    idx = [response(r).muscle(msubset).threshIdx];
    for i = 1:length(idx)
        m = msubset(i);
        p2p_EMGthresh(r,i) = response(r).muscle(m).p2pResponse(idx(i));
        p2p_Sensthresh(r,i) = mean(response(r).muscle(m).p2pResponse(response(r).StimAmps == sensoryThresh));
        EMGThreshRatio(r,i) = p2p_EMGthresh(r,i)/max(response(r).muscle(m).p2pResponse);
        SenseThreshRatio(r,i) = p2p_Sensthresh(r,i)/max(response(r).muscle(m).p2pResponse);
        ActiRatio(r,i) = max(response(r).muscle(m).p2pResponse)/min(response(r).muscle(m).p2pResponse);
    end
end

%% Plotting Spiders - at EMG vs Sensory Threshold
hS2= figure; maximize;
spdrfill = 'on';
spdrmax(1:8) = 3200; %max(max(p2p_thresh)); 
spdrmin(1:8) = 20; %min(min(p2p_thresh)); 
spdriter = 3; 

C = linspecer(length(response));

hS3= subplot(121);

    spider_plot(p2p_EMGthresh,...
        'AxesLabels', musclesL(msubset-8),'AxesInterval', 2,...
        'FillOption', spdrfill,'FillTransparency', 0.1, 'AxesLimits',...
        [spdrmin(1:length(msubset)); spdrmax(1:length(msubset))],...
        'AxesInterval', spdriter, 'Direction', 'clockwise',...
        'Color', C);
%     title('Residual Limb (left)');
    % Legend properties
    legend({response(:).setDescrpt}, 'Location', 'westoutside', 'Orientation', 'vertical','Interpreter','none');
     

    tmptit = title('P2P at PRM Threshold');
    set(tmptit, 'Interpreter', 'none');
% saveas(hS3,'D:\FigRescources\UH3\LSP05\rehash\PAD\P2P_EMGThresh_e9and16.png')

hS4= subplot(122);

    spider_plot(p2p_Sensthresh,...
        'AxesLabels', musclesL(msubset-8),'AxesInterval', 2,...
        'FillOption', spdrfill,'FillTransparency', 0.1, 'AxesLimits',...
        [spdrmin(1:length(msubset)); spdrmax(1:length(msubset))],...
        'AxesInterval', spdriter, 'Direction', 'clockwise',...
        'Color', C);
%     title('Residual Limb (left)');
    % Legend properties
    legend({response(:).setDescrpt}, 'Location', 'westoutside', 'Orientation', 'vertical','Interpreter','none');
     

    tmptit = title('P2P at SENSORY Threshold');
    set(tmptit, 'Interpreter', 'none');

    sgtitle('Activation Profiles at Threshold (1Hz,500us)');
saveas(hS2,['D:\FigRescources\UH3\LSP05\rehash\Stability\P2P_atThreshes_e' num2str(elecs(1)) '.png'])
saveas(hS2,['D:\FigRescources\UH3\LSP05\rehash\Stability\P2P_atThreshes_e' num2str(elecs(1)) '.svg'])
 


%% Plotting Spiders - at EMG vs Sensory Threshold - Activation Ratio
hS5= figure; maximize;
spdrfill = 'on';
spdrmax(1:8) = 1; %max(max(p2p_thresh)); 1200
spdrmin(1:8) = 0; %min(min(p2p_thresh)); 20
spdriter = 5; 

C = linspecer(length(response));

hS6= subplot(121);

    spider_plot(EMGThreshRatio,...
        'AxesLabels', musclesL(msubset-8),'AxesInterval', 2,...
        'FillOption', spdrfill,'FillTransparency', 0.1, 'AxesLimits',...
        [spdrmin(1:length(msubset)); spdrmax(1:length(msubset))],...
        'AxesInterval', spdriter, 'Direction', 'clockwise',...
        'Color', C);
%     title('Residual Limb (left)');
    % Legend properties
    legend({response(:).setDescrpt}, 'Location', 'westoutside', 'Orientation', 'vertical','Interpreter','none');
     

    tmptit = title('P2P at PRM Threshold');
    set(tmptit, 'Interpreter', 'none');
% saveas(hS3,'D:\FigRescources\UH3\LSP05\rehash\PAD\P2P_EMGThresh_e9and16.png')

hS7= subplot(122);

    spider_plot(SenseThreshRatio,...
        'AxesLabels', musclesL(msubset-8),'AxesInterval', 2,...
        'FillOption', spdrfill,'FillTransparency', 0.1, 'AxesLimits',...
        [spdrmin(1:length(msubset)); spdrmax(1:length(msubset))],...
        'AxesInterval', spdriter, 'Direction', 'clockwise',...
        'Color', C);
%     title('Residual Limb (left)');
    % Legend properties
    legend({response(:).setDescrpt}, 'Location', 'westoutside', 'Orientation', 'vertical','Interpreter','none');
     

    tmptit = title('P2P at SENSORY Threshold');
    set(tmptit, 'Interpreter', 'none');

    sgtitle('% Activation at Threshold (1Hz,500us)');
saveas(hS5,['D:\FigRescources\UH3\LSP05\rehash\Stability\ActiRatio_atThreshes_e' num2str(elecs(1)) '.png'])
saveas(hS5,['D:\FigRescources\UH3\LSP05\rehash\Stability\ActiRatio_atThreshes_e' num2str(elecs(1)) '.svg'])

%% Plotting Spiders - Activation Ratio
hS5= figure; maximize;
spdrfill = 'on';
spdrmax(1:8) = 1500; %max(max(p2p_thresh)); 1200
spdrmin(1:8) = 20; %min(min(p2p_thresh)); 20
spdriter = 5; 

C = linspecer(length(response));

    spider_plot(ActiRatio,...
        'AxesLabels', musclesL(msubset-8),'AxesInterval', 2,...
        'FillOption', spdrfill,'FillTransparency', 0.1, 'AxesLimits',...
        [spdrmin(1:length(msubset)); spdrmax(1:length(msubset))],...
        'AxesInterval', spdriter, 'Direction', 'clockwise',...
        'Color', C);

    % Legend properties
    legend({response(:).setDescrpt}, 'Location', 'westoutside', 'Orientation', 'vertical','Interpreter','none');
     

    tmptit = title('Selectivity by Activation Ratio');
    set(tmptit, 'Interpreter', 'none');
saveas(hS5,['D:\FigRescources\UH3\LSP05\rehash\Stability\ActiRatio_e' num2str(elecs(1)) '.png'])
saveas(hS5,['D:\FigRescources\UH3\LSP05\rehash\Stability\ActiRatio_e' num2str(elecs(1)) '.svg'])



%% Plot Onset and AUC

%% PCA for AUC
for iD = 1:length(response)
    for iM = 1 : length(msubset)
        xinfo = response(iD).muscle(msubset(iM)).p2pResponse;
        ampinfo = [ 0 response(iD).StimAmps];
        aucData(iD,iM) = sum(xinfo.*diff(ampinfo));
        onsets(iD,iM) = mean(response(iD).muscle(msubset(iM)).onset(response(iD).muscle(msubset(iM)).threshIdx:end));
        tmpOns = response(iD).muscle(msubset(iM)).onset(response(iD).muscle(msubset(iM)).threshIdx:end);
        stderror(iD,iM) = std( tmpOns ) / sqrt( length(tmpOns) );
    end
end
hs6 = figure;maximize;
tiledlayout(1,3);
nexttile
bar([1:8],aucData);ylabel('Cum EMG (uV*mA)');
legend({response(:).setDescrpt}, 'Location', 'northwest', 'Orientation', 'vertical','Interpreter','none');
title('Gross Recruitment');
xticklabels(musclesL);

nexttile
bar([1:8],onsets);ylabel('Response Time (ms)')
xticklabels(musclesL);
legend({response(:).setDescrpt}, 'Location', 'northwest', 'Orientation', 'vertical','Interpreter','none');
title('Response Onset');

nexttile
boxplot(onsets)

saveas(hS6,['D:\FigRescources\UH3\LSP05\rehash\Stability\AUC_Onset' num2str(elecs(1)) '.png'])
saveas(hS6,['D:\FigRescources\UH3\LSP05\rehash\Stability\AUC_Onset' num2str(elecs(1)) '.svg'])


