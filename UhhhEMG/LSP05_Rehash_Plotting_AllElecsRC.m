% % %%%% UH3 EMG - Recruitment Curves - Quick Load Trial Stats - Testing Version B %%
% % 
for lcall = 4%:4

    switch lcall
        case 1
            % elecs 1-8: day9_set2
            load('D:\DATA\UH3 testing\LSP05\data_gen\Day12-Set_1_E1-8RC_1-6mA_500us_1Hz');
            multipolar = 'no';
            anodeElec = [];
            disp('elecs 1-8: day9_set2');
        case 2
            % elecs 9-16: Day8-Set_1_RC_1-6mA_500us_1Hz
            load('D:\DATA\UH3 testing\LSP05\data_gen\Day8-Set_1_RC_1-6mA_500us_1Hz.mat');
            multipolar = 'no';
            anodeElec = [];
            disp('elecs 9-16: Day8-Set_1_RC_1-6mA_500us_1Hz');
        case 3
            % elecs 18-24: Day8-Set_2_RC_1-6mA_500us_1Hz
            load('D:\DATA\UH3 testing\LSP05\data_gen\Day8-Set_2_RC_1-6mA_500us_1Hz.mat');
            multipolar = 'no';
            anodeElec = [];
            disp('elecs 18-24: Day8-Set_2_RC_1-6mA_500us_1Hz');
        case 4
            % elecs 25-32: day9_set1
            load('D:\DATA\UH3 testing\LSP05\data_gen\day9-set1.mat');
            disp('elecs 25-32: day9_set1');
        otherwise
            % multipolar: day9_set3
            load('D:\DATA\UH3 testing\LSP05\data_gen\day9-set3.mat');
            disp('multipolar: day9_set3');
    end
end
%% Basic Recruitment Curve Stuff for Frequency and Set on all Elecs
    freqs = unique([stimFrequency]);
    pwidths = unique([pulseWidth]);
    elecs = unique([spinalElec]);
    
% %     freqs = unique([trial.stimFrequency]);
% %     pwidths = unique([trial.pulseWidth]);
% %     elecs = unique([trial.spinalElec]);

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

% Swapped RF & VL so 3<->2 & 11<->10
msubset = [1:16]; %or 1:8 for right/intact 1:16 for 
% msubset([2 3]) = msubset([3 2]);

%% Make RC vectors for each subset.
check_p = 0.5;
check_f = 1;
artifactBuffer = 0.005;
disp(['Processing: ' num2str(check_p) 'ms ' num2str(check_f) 'Hz']);

for r = 1:length(elecs)
    response(r).elec = elecs(r);
    response(r).freq = check_f;
    response(r).pulse = check_p; 
    
    eIdx = find(([stimFrequency] == check_f) & ([pulseWidth] == check_p) & ([spinalElec] == elecs(r)));
% %     eIdx = find([trial.spinalElec] == elecs(r));
    [sortedStim,I] = sort([stimAmp(eIdx)]);
% %     [sortedStim,I] = sort([trial(eIdx).stimAmp]);

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
% %     for i = 1:length(msubset)
% %         m = msubset(i);
% %         response(r).muscle(m).threshIdx = findchangepts(cumsum(response(r).muscle(m).p2pResponse), 'Statistic','rms');
% %         response(r).muscle(m).threshold = (sortedStim(response(r).muscle(m).threshIdx));
% %     end

end

%% Get Thresholds
for r = 1:length(elecs)
    for m = 1:length(msubset)
        test_diff = response(r).muscle(m).p2pResponse(end) <= (response(r).muscle(m).p2pResponse(1));
        test_mean = response(r).muscle(m).p2pResponse(end) <= mean(response(r).muscle(m).p2pResponse);
        test_std = 1 > ((response(r).muscle(m).p2pResponse(end) - mean(response(r).muscle(m).p2pResponse))/(2*std(response(r).muscle(m).p2pResponse))) ;

        tmp = test_diff || test_mean; %|| test_std (i);
        response(r).muscle(m).exist = ~tmp;
        if tmp
            response(r).muscle(m).threshIdx = nan;
            response(r).muscle(m).threshold = nan;
        else
            
            for i = 2:length(response(r).muscle(m).p2pResponse)
                x = mean(response(r).muscle(m).p2pResponse(1:i));
                tt(i) = response(r).muscle(m).p2pResponse(i) > (x + std(response(r).muscle(m).p2pResponse(1:i)));
            end
            tt(1) = 0;
            for  i = 3:length(response(r).muscle(m).p2pResponse)
                if tt(i-1) == 0 && tt(i-2) == 1
                    tt(i-2) = 0;
                end
            end
            for  i = 3:length(response(r).muscle(m).p2pResponse)
                if tt(i-1) == 0 && tt(i-2) == 1
                    tt(i-2) = 0;
                end
            end
            idx = find(tt);
            if isempty(idx)
                response(r).muscle(m).threshIdx = nan;
                response(r).muscle(m).threshold = nan;
            else
                response(r).muscle(m).threshIdx = idx(1);
                response(r).muscle(m).threshold = response(r).StimAmps(idx(1));
            end
            clear tt
        end
    end
end
    
% Plot CUm Sum vs P2P then adjust





%Save Response Variable for later


%% Plot RC for Peak2Peak, RMS, and Max Peak (Abs)
setPath = 'D:\FigRescources\UH3\LSP05\rehash\VL_BF_LG\';
% Using the saved response variable created above
for r = 1:length(elecs)
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
    tit = sgtitle(['E' num2str(elecs(r)) ' - RCs plotted 3 ways (' num2str(check_p) 'ms ' num2str(check_f) 'Hz)']);
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
% % 

%% Plot CumSum vs P2P
for r = 1:length(elecs)
    hF(r) = figure; maximize;
    for i = 1:length(msubset)
        m = msubset(i);
        h2(i) = subplot(2,length(msubset)/2,i);
        plot(response(r).StimAmps/1000, response(r).muscle(m).p2pResponse);
        hold on; plot(response(r).StimAmps/1000, cumsum(response(r).muscle(m).p2pResponse));
        plot(response(r).StimAmps/1000,smooth(cumsum(response(r).muscle(m).p2pResponse)));
        title(response(r).muscle(m).Muscle_ID);
        ylabel('EMG (uV)');
        xlabel('StimAmp (mA)');

        legend('P2P','cumsum', 'smooth');
    end
    % linkaxes([h(:)],'xy');
    tit = sgtitle(['E' num2str(elecs(r)) ' - PS2 vs Cumsum (' num2str(check_p) 'ms ' num2str(check_f) 'Hz)']);
    text = tit.String
    
    saveas(hF(r),[setPath tit.String '.png'])
end

% % %% Plot Change Points
% % for r = 1:length(elecs)
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
msubset = [11 12 16];
% % elecs = [1 8 16 18 23 29];
C = linspecer(length(elecs));

for i = 1:length(msubset)
    m = msubset(i);
    h(i) = subplot(1,length(msubset),i);
    for r = 1:length(elecs)
        e = find([response.elec] == elecs(r));
        tmp(e) = max(response(e).muscle(m).p2pResponse);
    end
    maxR(m) = max(tmp);
end

hF = figure; maximize;
for i = 1:length(msubset)
    m = msubset(i);
    h(i) = subplot(1,length(msubset),i);
    for r = 1:length(elecs)
        e = find([response.elec] == elecs(r));
        plot(response(e).StimAmps/1000, response(e).muscle(m).p2pResponse/maxR(m),'Color',C(r,:), 'LineWidth',2);
        hold on;
    end
    title(response(e).muscle(m).Muscle_ID);
    ylabel('EMG (uV)');
    xlabel('StimAmp (mA)');
        
    legend(arrayfun(@num2str,elecs,'UniformOutput',false),'Location','northwest');
end
% linkaxes([h(:)],'xy');
tit = sgtitle(['Recruitment Curves by Muscle for Electrodes ' num2str(elecs(1)) '-' num2str(elecs(end))]);
% % saveas(hF,['D:\FigRescources\UH3\LSP05\rehash\PAD\P2P_RCs_e' num2str(elecs(1)) '-' num2str(elecs(end)) '.png'])
% % savefig(hF,['D:\FigRescources\UH3\LSP05\rehash\PAD\P2P_RCs_e' num2str(elecs(1)) '-' num2str(elecs(end))])
% % 
% % close all;

%% Plot Comparison RC Stacked
% % clear all;
% % 
% % load('D:\FigRescources\UH3\LSP05\rehash\responses_500us_1Hz_allElecsRC.mat')
% % response = response2; clear response2;

labelmuscles = { 'VM', 'VL', 'RF', 'BF', 'ST', 'Ham', 'MG', 'LG'}; %Ham vs TA
msubset = [9 11 10 12 13 14 15 16]; %[11 12 16 14];
% % etmp= [response.elec];
% % elecs = etmp(1:4:end);
% % elecs = [response(15:end).elec];%[response.elec]; %[1 4 7 9 11 18 13 20 16 23 26 29]; %
C = linspecer(length(elecs));

for i = 1:length(msubset)
    m = msubset(i);
    for r = 1:length(elecs)
        e = find([response.elec] == elecs(r));
        tmp(e) = max(response(e).muscle(m).p2pResponse);
    end
    musc_map(m,:) = tmp;
    maxR(m) = max(tmp);
end

hF = figure;% maximize;
for i = 1:length(msubset)
    m = msubset(i);
    
    for r = 1:length(elecs)
        h(i,r) = subplot(length(elecs),length(msubset),length(msubset)*(r-1)+i);
        e = find([response.elec] == elecs(r));
        p = area(response(e).StimAmps/1000, response(e).muscle(m).p2pResponse); %/maxR(m)
        p.FaceColor = C(r,:);
        hold on;
        if i==1
            ylabel('EMG (uV)');
        end
        if r==length(elecs)
            xlabel('StimAmp (mA)');
        end
        if r==1
            title(response(e).muscle(m).Muscle_ID);
        end
%         legend(['E' num2str(elecs(r))],'Location','northwest');
        box off
        grid off
        set(gca,'Visible','off')
        axis off;
        h(i,r).Color = 'none';
    end
    
     
    
end
% % for i = 1:length(msubset)
% % linkaxes([h(i,:)],'xy');
% % end
linkaxes([h],'xy');
tit = sgtitle(['Recruitment Curves by Muscle for Electrodes ' num2str(elecs(1)) '-' num2str(elecs(end))]);
set(0,'defaultAxesFontSize',12)
% % saveas(hF,['D:\FigRescources\UH3\LSP05\rehash\PAD\Stacked_RCs_e' num2str(elecs(1)) '-' num2str(elecs(end)) '.png'])
% % savefig(hF,['D:\FigRescources\UH3\LSP05\rehash\PAD\Stacked_RCs_e' num2str(elecs(1)) '-' num2str(elecs(end))])
% % 
% % close all;



%% Plot Threshold Comparisons
musclesL = {'VM', 'VL', 'RF', 'BF', 'ST', 'TA', 'MG', 'LG',};
div = 1000;


for r = 1:length(response)
    threshold(r,:) = [response(r).muscle(msubset).threshold]/div;
end

hS1 = figure; maximize;
spdrfill = 'off';
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
    legend(arrayfun(@num2str,elecs,'UniformOutput',false), 'Location', 'westoutside', 'Orientation', 'vertical');
     

    tmptit = sgtitle(['Thresholds for ' num2str(elecs(1)) '-' num2str(elecs(end)) ' (1Hz, 500us)']);
    set(tmptit, 'Interpreter', 'none');

% tit2 = sgtitle('Threshold by Muscle for Electrodes 9 & 16');
saveas(hS1,['D:\FigRescources\UH3\LSP05\rehash\PAD\Threshold_e' num2str(elecs(1)) '-' num2str(elecs(end)) '.png']);
savefig(hS1,['D:\FigRescources\UH3\LSP05\rehash\PAD\Threshold_e' num2str(elecs(1)) '-' num2str(elecs(end))]);

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
spdrmax(1:8) = 1200; %max(max(p2p_thresh)); 
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
    legend(arrayfun(@num2str,elecs,'UniformOutput',false), 'Location', 'westoutside', 'Orientation', 'vertical');
     

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
    legend(arrayfun(@num2str,elecs,'UniformOutput',false), 'Location', 'westoutside', 'Orientation', 'vertical');
     

    tmptit = title('P2P at SENSORY Threshold');
    set(tmptit, 'Interpreter', 'none');

    sgtitle('Activation Profiles at Threshold (1Hz,500us)');
saveas(hS2,['D:\FigRescources\UH3\LSP05\rehash\PAD\P2P_atThreshes_e' num2str(elecs(1)) '-' num2str(elecs(end)) '.png'])
savefig(hS2,['D:\FigRescources\UH3\LSP05\rehash\PAD\P2P_atThreshes_e' num2str(elecs(1)) '-' num2str(elecs(end))])
 


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
    legend(arrayfun(@num2str,elecs,'UniformOutput',false), 'Location', 'westoutside', 'Orientation', 'vertical');
     

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
    legend(arrayfun(@num2str,elecs,'UniformOutput',false), 'Location', 'westoutside', 'Orientation', 'vertical');
     

    tmptit = title('P2P at SENSORY Threshold');
    set(tmptit, 'Interpreter', 'none');

    sgtitle('% Activation at Threshold (1Hz,500us)');
saveas(hS5,['D:\FigRescources\UH3\LSP05\rehash\PAD\ActiRatio_atThreshes_e' num2str(elecs(1)) '-' num2str(elecs(end)) 'png'])
savefig(hS5,['D:\FigRescources\UH3\LSP05\rehash\PAD\ActiRatio_atThreshes_e' num2str(elecs(1)) '-' num2str(elecs(end))])

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
    legend(arrayfun(@num2str,elecs,'UniformOutput',false), 'Location', 'westoutside', 'Orientation', 'vertical');
     

    tmptit = title('Selectivity by Activation Ratio');
    set(tmptit, 'Interpreter', 'none');
saveas(hS5,['D:\FigRescources\UH3\LSP05\rehash\PAD\ActiRatio_e' num2str(elecs(1)) '-' num2str(elecs(end)) 'png'])
savefig(hS5,['D:\FigRescources\UH3\LSP05\rehash\PAD\ActiRatio_e' num2str(elecs(1)) '-' num2str(elecs(end))])


%% Plot Onsets and Thresholds
musclesL = {'VM', 'RF', 'VL', 'BF', 'ST', 'TA', 'MG', 'LG',};
msubset = [10 12 16 14];
elecs = [4 7 9 11 18 13 20 16 23 26];
C = linspecer(length(elecs),'qualitative');

for r = 1:length(elecs)
    e = find([response2.elec] == elecs(r));
    for i = 1:length(msubset)
        m = msubset(i);  
        tmpOns = [response2(e).muscle(m).onset(response2(e).muscle(m).threshIdx:end)];
        onset(r,i) = nanmean(tmpOns);
        onserror(r,i) = nanstd( tmpOns ) / sqrt( length(tmpOns) );
        
        threshes(r,i) = response2(e).muscle(m).threshold/1000;
    end
end

%% Plot
hS6 = figure;%maximize;
% % tiledlayout(1,2);
% % nexttile

superbar(onset','E', onserror','BarFaceColor',permute(C, [3 1 2]));
ylabel('Response Onset (ms)');
% tmp = [' ', mlab, ' '];
% xticklabels(tmp);
xticks([1:4]);
xticklabels(musclesL(msubset-8));
legend(arrayfun(@num2str,elecs,'UniformOutput',false), 'Location', 'eastoutside', 'Orientation', 'vertical');
hold on;
% % 
% % nexttile
hS7 = figure;
superbar(threshes','BarFaceColor',permute(C, [3 1 2]));
ylabel('Threshold Amplitude (mA)');
% tmp = [' ', mlab, ' '];
% xticklabels(tmp);
xticks([1:4]);ylim([0,6])
xticklabels(musclesL(msubset-8));
legend(arrayfun(@num2str,elecs,'UniformOutput',false), 'Location', 'eastoutside', 'Orientation', 'vertical');
hold on;

saveas(hS6,['D:\FigRescources\UH3\LSP05\rehash\Onset_barplot.png']);
saveas(hS6,['D:\FigRescources\UH3\LSP05\rehash\Onset_barplot.svg']);

saveas(hS7,['D:\FigRescources\UH3\LSP05\rehash\Threshes_barplot.png']);
saveas(hS7,['D:\FigRescources\UH3\LSP05\rehash\Threshes_barplot.svg']);
