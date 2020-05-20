%%%% UH3 EMG - Recruitment Curves - Quick Load Trial Stats - Testing Version B %%
%% Basic Recruitment Curve Stuff for Frequency and Set on Elec 9 and 16
    freqs = unique([trial.stimFrequency]);
    pwidths = unique([trial.pulseWidth]);
    elecs = unique([trial.spinalElec]);

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
msubset = [9:16]; %or 1:8 for right/intact
msubset([2 3]) = msubset([3 2]);

%% Make RC vectors for each subset.
r = 2; %1 for e9, 2 for e16
check_p = 0.5;
check_f = 1;
disp(['Processing: ' num2str(check_p) 'ms ' num2str(check_f) 'Hz']);

response(r).elec = trial(1).spinalElec;
response(r).freq = check_f;
response(r).pulse = check_p;

eIdx = find(([trial.stimFrequency] == check_f) & ([trial.pulseWidth] == check_p));
[sortedStim,I] = sort([trial(eIdx).stimAmp]);

response(r).StimAmps = sortedStim;
%sensoryThresh = 2000

for i = 1:length(msubset)
    m = msubset(i);
    disp(['Processing ' Muscle_ID(m)]);
    response(r).muscle(m).Muscle_ID = Muscle_ID(m);
    
    for t = 1:length(eIdx)
        n = eIdx(I(t));
        
        response(r).muscle(m).responseRMS(t) = trial(n).responseRMS(m);
        response(r).muscle(m).p2pResponse(t) = peak2peak(trial(n).meanTrace{m}(stimEnd:end));
        [response(r).muscle(m).peakResponse(t), mi] = max(abs(trial(n).meanTrace{m}(stimEnd:end)));
        response(r).muscle(m).onset(t) = trial(n).epochTimeVec(mi+stimEnd-1);
        response(r).muscle(m).meanbase(t) = trial(n).meanbase(m,:);
        response(r).muscle(m).rmsbase(t) = rms(trial(n).baseline(m,:));
        
        disp([num2str(t) ' of ' num2str(length(eIdx))]);
    end

end

%% Get Threshold
for i = 1:length(msubset)
    m = msubset(i);
    hf4(i) = figure; maximize;
    response(r).muscle(m).threshIdx = findchangepts(cumsum(response(r).muscle(m).p2pResponse), 'Statistic','rms');
    response(r).muscle(m).threshold = (sortedStim(response(r).muscle(m).threshIdx));
end

%Save Response Variable for later


%% Plot RC for Peak2Peak, RMS, and Max Peak (Abs)

% Using the saved response variable created above

% % hF(1) = figure; maximize;
% % for i = 1:length(msubset)
% %     m = msubset(i);
% %     h(i) = subplot(2,length(msubset)/2,i);
% %     plot(response(r).StimAmps/1000, response(r).muscle(m).peakResponse);
% %     hold on; plot(response(r).StimAmps/1000, response(r).muscle(m).p2pResponse);
% %     hold on; plot(response(r).StimAmps/1000, response(r).muscle(m).responseRMS);
% %     title(response(r).muscle(m).Muscle_ID);
% %     ylabel('EMG (uV)');
% %     xlabel('StimAmp (mA)');
% %     
% %     yyaxis right
% %     plot(response(r).StimAmps/1000, response(r).muscle(m).onset*1000,'-.');
% %     ylabel('Onset (ms)');
% %     
% %     legend('Abs Peak', 'RMS', 'P2P', 'Onset');
% % end
% % % linkaxes([h(:)],'xy');
% % tit{1} = sgtitle(['RCs plotted 3 ways (' num2str(check_p) 'ms ' num2str(check_f) 'Hz)']);
% % text{1} = tit{1}.String
% % 
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
% % %% Plot CumSum vs P2P
% % hF(3) = figure; maximize;
% % for i = 1:length(msubset)
% %     m = msubset(i);
% %     h2(i) = subplot(2,length(msubset)/2,i);
% %     plot(response(r).StimAmps/1000, response(r).muscle(m).p2pResponse);
% %     hold on; plot(response(r).StimAmps/1000, cumsum(response(r).muscle(m).p2pResponse));
% %     title(response(r).muscle(m).Muscle_ID);
% %     ylabel('EMG (uV)');
% %     xlabel('StimAmp (mA)');
% %     
% %     legend('P2P','cumsum');
% % end
% % % linkaxes([h(:)],'xy');
% % tit{3} = sgtitle(['PS2 vs Cumsum (' num2str(check_p) 'ms ' num2str(check_f) 'Hz)']);
% % text{3} = tit{3}.String
% % 
% % for f = 1:3
% %     saveas(hF(f),[setPath setDescrpt '-' tit{f}.String '.png'])
% % end
% % %% Plot Change Points
% % 
% % for i = 1:length(msubset)
% %     m = msubset(i);
% %     hf4(i) = figure; maximize;
% %     findchangepts(cumsum(response(r).muscle(m).p2pResponse), 'Statistic','rms');
% %     htit{i} = title(char(response(r).muscle(m).Muscle_ID));
% %     ylabel('EMG (uV)');
% %     xlabel('StimAmp (mA)');
% %     
% %     legend('P2P','cumsum');
% % end
% % for f = 1:length(hf4)
% %     saveas(hf4(f),[setPath setDescrpt '-' htit{f}.String '-Inflection.png'])
% % end
% % 
% % close all;

%% Plot Comparison RC
% % C = linspecer(length(response));
% % 
% % msubset = [9:16];
% % hF = figure; maximize;
% % for i = 1:length(msubset)
% %     m = msubset(i);
% %     h(i) = subplot(2,length(msubset)/2,i);
% %     for r = 1:length(response)
% %         plot(response(r).StimAmps/1000, response(r).muscle(m).p2pResponse,'Color',C(r,:));
% %         hold on;
% %     end
% %     title(response(r).muscle(m).Muscle_ID);
% %     ylabel('EMG (uV)');
% %     xlabel('StimAmp (mA)');
% %         
% %     legend('e9', 'e16','Location','northwest');
% % end
% % % linkaxes([h(:)],'xy');
% % tit = sgtitle('Recruitment Curves by Muscle for Electrodes 9 & 16');
% % saveas(hF,'D:\FigRescources\UH3\LSP05\rehash\PAD\P2P_RCs_e9and16.png')
% % savefig(hF,'D:\FigRescources\UH3\LSP05\rehash\PAD\P2P_RCs_e9and16')

%% Plot Threshold Comparisons
% % musclesL = {'VM', 'VL', 'RF', 'BF', 'ST', 'TA', 'MG', 'LG',};
% % div = 1000;
% % 
% % 
% % for r = 1:length(response)
% %     threshold(r,:) = [response(r).muscle(msubset).threshold]/div;
% % end
% % 
% % hS1 = figure; maximize;
% % spdrfill = 'off';
% % spdrmax(1:8) =((6000/div)); 
% % spdrmin(1:8) = ((300/div)); 
% % spdriter = 6; 
% % 
% % C = linspecer(length(response));
% % 
% %     spider_plot(threshold,...
% %         'AxesLabels', musclesL,'AxesInterval', 2,...
% %         'FillOption', spdrfill,'FillTransparency', 0.1, 'AxesLimits',...
% %         [spdrmin; spdrmax],...
% %         'AxesInterval', spdriter, 'Direction', 'clockwise',...
% %         'Color', C);
% % %     title('Residual Limb (left)');
% %     % Legend properties
% %     legend(['E09'; 'E16';], 'Location', 'westoutside', 'Orientation', 'vertical');
% %      
% % 
% %     tmptit = sgtitle('Thresholds for E9 and E16: 1Hz, 500us');
% %     set(tmptit, 'Interpreter', 'none');

% tit2 = sgtitle('Threshold by Muscle for Electrodes 9 & 16');
% % saveas(hS1,'D:\FigRescources\UH3\LSP05\rehash\PAD\Threshold_e9and16.png')
% % savefig(hS1,'D:\FigRescources\UH3\LSP05\rehash\PAD\Threshold_e9and16')

%% Activation at Threshold vs Activation Ratio
sensoryThresh = [3000;2000];

for r = 1:length(response)
    idx = find(response(r).StimAmps==min([response(r).muscle(msubset).threshold])); %[response(r).muscle(msubset).threshIdx];
    for i = 1:length(msubset)
        m = msubset(i);
        p2p_EMGthresh(r,i) = response(r).muscle(m).p2pResponse(idx(1));
        p2p_Sensthresh(r,i) = mean(response(r).muscle(m).p2pResponse(response(r).StimAmps == sensoryThresh(r)));
        EMGThreshRatio(r,i) = p2p_EMGthresh(r,i)/max(response(r).muscle(m).p2pResponse);
        SenseThreshRatio(r,i) = p2p_Sensthresh(r,i)/max(response(r).muscle(m).p2pResponse);
        ActiRatio(r,i) = max(response(r).muscle(m).p2pResponse)/min(response(r).muscle(m).p2pResponse);
    end
end
display('done processing');
% % %% Plotting Spiders - at EMG vs Sensory Threshold
% % hS2= figure; %maximize;
% % spdrfill = 'on';
% % spdrmax(1:8) = 1200; %max(max(p2p_thresh)); 
% % spdrmin(1:8) = 20; %min(min(p2p_thresh)); 
% % spdriter = 2; 
% % 
% % C = linspecer(length(response));
% % 
% % % hS3= subplot(121);
% % 
% %     spider_plot(p2p_EMGthresh,...
% %         'AxesLabels', musclesL,'AxesInterval', 2,...
% %         'FillOption', spdrfill,'FillTransparency', 0.1, 'AxesLimits',...
% %         [spdrmin; spdrmax],...
% %         'AxesInterval', spdriter, 'Direction', 'clockwise',...
% %         'Color', C);
% % %     title('Residual Limb (left)');
% %     % Legend properties
% %     legend(['E09'; 'E16';], 'Location', 'westoutside', 'Orientation', 'vertical');
% %      
% % 
% %     tmptit = title('P2P at PRM Threshold');
% %     set(tmptit, 'Interpreter', 'none');
% % % saveas(hS3,'D:\FigRescources\UH3\LSP05\rehash\PAD\P2P_EMGThresh_e9and16.png')
% % 
% % % hS4= subplot(122);
% % 
% %     spider_plot(p2p_Sensthresh,...
% %         'AxesLabels', musclesL,'AxesInterval', 2,...
% %         'FillOption', spdrfill,'FillTransparency', 0.1, 'AxesLimits',...
% %         [spdrmin; spdrmax],...
% %         'AxesInterval', spdriter, 'Direction', 'clockwise',...
% %         'Color', C);
% % %     title('Residual Limb (left)');
% %     % Legend properties
% %     legend(['E09, 3mA'; 'E16, 2mA';], 'Location', 'southwestoutside', 'Orientation', 'vertical');
% %      
% % 
% %     tmptit = title('P2P at SENSORY Threshold');
% %     set(tmptit, 'Interpreter', 'none');
% % 
% % % %     sgtitle('Activation Profiles at Threshold (1Hz,500us)');
% % saveas(hS2,'D:\FigRescources\UH3\LSP05\rehash\PAD\P2P_atThreshes_e9and16.png')
% % savefig(hS2,'D:\FigRescources\UH3\LSP05\rehash\PAD\P2P_atThreshes_e9and16')
% %  


% % %% Plotting Spiders - at EMG vs Sensory Threshold - Activation Ratio
% % hS5= figure; maximize;
% % spdrfill = 'on';
% % spdrmax(1:8) = 1; %max(max(p2p_thresh)); 1200
% % spdrmin(1:8) = 0; %min(min(p2p_thresh)); 20
% % spdriter = 5; 
% % 
% % C = linspecer(length(response));
% % 
% % hS6= subplot(121);
% % 
% %     spider_plot(EMGThreshRatio,...
% %         'AxesLabels', musclesL,'AxesInterval', 2,...
% %         'FillOption', spdrfill,'FillTransparency', 0.1, 'AxesLimits',...
% %         [spdrmin; spdrmax],...
% %         'AxesInterval', spdriter, 'Direction', 'clockwise',...
% %         'Color', C);
% % %     title('Residual Limb (left)');
% %     % Legend properties
% %     legend(['E09'; 'E16';], 'Location', 'westoutside', 'Orientation', 'vertical');
% %      
% % 
% %     tmptit = title('P2P at PRM Threshold');
% %     set(tmptit, 'Interpreter', 'none');
% % % saveas(hS3,'D:\FigRescources\UH3\LSP05\rehash\PAD\P2P_EMGThresh_e9and16.png')
% % 
% % hS7= subplot(122);
% % 
% %     spider_plot(SenseThreshRatio,...
% %         'AxesLabels', musclesL,'AxesInterval', 2,...
% %         'FillOption', spdrfill,'FillTransparency', 0.1, 'AxesLimits',...
% %         [spdrmin; spdrmax],...
% %         'AxesInterval', spdriter, 'Direction', 'clockwise',...
% %         'Color', C);
% % %     title('Residual Limb (left)');
% %     % Legend properties
% %     legend(['E09, 3mA'; 'E16, 2mA';], 'Location', 'westoutside', 'Orientation', 'vertical');
% %      
% % 
% %     tmptit = title('P2P at SENSORY Threshold');
% %     set(tmptit, 'Interpreter', 'none');
% % 
% %     sgtitle('% Activation at Threshold (1Hz,500us)');
% % saveas(hS5,'D:\FigRescources\UH3\LSP05\rehash\PAD\ActiRatio_atThreshes_e9and16.png')
% % savefig(hS5,'D:\FigRescources\UH3\LSP05\rehash\PAD\ActiRatio_atThreshes_e9and16')

% % %% Plotting Spiders - Activation Ratio
% % hS5= figure; maximize;
% % spdrfill = 'on';
% % spdrmax(1:8) = 1500; %max(max(p2p_thresh)); 1200
% % spdrmin(1:8) = 20; %min(min(p2p_thresh)); 20
% % spdriter = 5; 
% % 
% % C = linspecer(length(response));
% % 
% %     spider_plot(ActiRatio,...
% %         'AxesLabels', musclesL,'AxesInterval', 2,...
% %         'FillOption', spdrfill,'FillTransparency', 0.1, 'AxesLimits',...
% %         [spdrmin; spdrmax],...
% %         'AxesInterval', spdriter, 'Direction', 'clockwise',...
% %         'Color', C);
% % 
% %     % Legend properties
% %     legend(['E09'; 'E16';], 'Location', 'westoutside', 'Orientation', 'vertical');
% %      
% % 
% %     tmptit = title('Selectivity by Activation Ratio');
% %     set(tmptit, 'Interpreter', 'none');
% % saveas(hS5,'D:\FigRescources\UH3\LSP05\rehash\PAD\ActiRatio_e9and16.png')
% % savefig(hS5,'D:\FigRescources\UH3\LSP05\rehash\PAD\ActiRatio_e9and16')

