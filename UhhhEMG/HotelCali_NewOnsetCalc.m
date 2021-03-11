%% Get Data
for lcall = 1%:4

    switch lcall
        case 1
            % elecs 1-8: day9_set2
            load('D:\DATA\UH3 testing\LSP05\data_gen\Day12-Set_1_E1-8RC_1-6mA_500us_1Hz');
            multipolar = 'no';
            anodeElec = [];
            disp('elecs 1-8: day9_set2');
            freqs = unique([stimFrequency]);
            pwidths = unique([pulseWidth]);
            elecs = unique([spinalElec]);
        case 2
            % elecs 9-16: Day8-Set_1_RC_1-6mA_500us_1Hz
            load('D:\DATA\UH3 testing\LSP05\data_gen\Day8-Set_1_RC_1-6mA_500us_1Hz.mat');
            multipolar = 'no';
            anodeElec = [];
            disp('elecs 9-16: Day8-Set_1_RC_1-6mA_500us_1Hz');
            freqs = unique([stimFrequency]);
            pwidths = unique([pulseWidth]);
            elecs = unique([spinalElec]);
        case 3
            % elecs 18-24: Day8-Set_2_RC_1-6mA_500us_1Hz
            load('D:\DATA\UH3 testing\LSP05\data_gen\Day8-Set_2_RC_1-6mA_500us_1Hz.mat');
            multipolar = 'no';
            anodeElec = [];
            disp('elecs 18-24: Day8-Set_2_RC_1-6mA_500us_1Hz');
            freqs = unique([stimFrequency]);
            pwidths = unique([pulseWidth]);
            elecs = unique([spinalElec]);
        case 4
            % elecs 25-32: day9_set1
            load('D:\DATA\UH3 testing\LSP05\data_gen\day9-set1.mat');
            disp('elecs 25-32: day9_set1');
            freqs = unique([stimFrequency]);
            pwidths = unique([pulseWidth]);
            elecs = unique([spinalElec]);
        otherwise
            % multipolar: day9_set3
% %             load('D:\DATA\UH3 testing\LSP05\data_gen\day9-set3.mat');
% %             disp('multipolar: day9_set3');
            load('D:\FigRescources\UH3\HotelCali\LSP02b\RCs_E1-12_33-35v2.mat')
            disp('LSP02bL: e1-12, 33-35');
            
            freqs = unique([trial.stimFrequency]);
            pwidths = unique([trial.pulseWidth]);
            elecs = unique([trial.spinalElec]);
    end
end

%%

    
msubset = [1:16];
    %% Make RC vectors for each subset.
check_p = 0.5;
check_f = 1;
artifactBuffer = 0.005;
dtPre = 50e-3; % msec before stim onset 
dtPost = 100e-3; % msec after stim onset 
rmsWindow = 0.05; % post stim RMS Window in s, so 0.05 is 50ms 

            
disp(['Processing: ' num2str(check_p) 'ms ' num2str(check_f) 'Hz']);

for r = 12%:length(elecs)
    meta(r).elec = elecs(r);
    meta(r).freq = check_f;
    meta(r).pulse = check_p; 
    
% %     eIdx = find(([stimFrequency] == check_f) & ([pulseWidth] == check_p) & ([spinalElec] == elecs(r)));
    eIdx = find([trial.spinalElec] == elecs(r));
% %     [sortedStim,I] = sort([stimAmp(eIdx)]);
    [sortedStim,I] = sort([trial(eIdx).stimAmp]);

    meta(r).StimAmps = sortedStim;
    
    for i = [2 10 3 11 4 12 5 13 9 16]%:length(msubset)
    m = msubset(i);
    disp(['Processing ' Muscle_ID(m)]);
    meta(r).muscle(m).Muscle_ID = Muscle_ID(m);
    
        for t = 1:length(eIdx)
            n = eIdx(I(t));
            stimBegins = trial(n).stims;
            stimEnd = nSampsPre + trial(n).stimLength + artifactBuffer*fs;  
            
% %             for s = 1:trial(n).numStims
% %                 close all;
% %                 figure;maximize;
% %                 plot(trial(n).epochTimeVec, trial(n).emgStim{1,m}(s,:));
% %                 hold on;
% %                 plot(trial(n).epochTimeVec, trial(n).meanTrace{1, m})
% %                 title([Muscle_ID(m) ' - trial ' num2str(t) ' stim ' num2str(s)]);
% %                 [x,y] = ginput(2);
% %                 
% %                 if isempty(x)
% %                     meta(r).muscle(m).onset(t,s) = NaN;
% %                     meta(r).muscle(m).peaktime(t,s) = NaN;
% %                     meta(r).muscle(m).peak(t,s) = NaN;   
% %                 else
% %                     meta(r).muscle(m).onset(t,s) = x(1);
% %                     meta(r).muscle(m).peaktime(t,s) = x(2);
% %                     meta(r).muscle(m).peak(t,s) = x(2);              
% %                 end
% %                                
% %             end
% %             meta(r).muscle(m).onset_mean(t) = nanmean([meta(r).muscle(m).onset(t,:)]);
% %             meta(r).muscle(m).onset_var(t) = nanvar([meta(r).muscle(m).onset(t,:)]);
% %             meta(r).muscle(m).onset_std(t) = nanstd([meta(r).muscle(m).onset(t,:)]);
% %             meta(r).muscle(m).peaktime_mean(t) = nanmean([meta(r).muscle(m).peaktime(t,:)]);
% %             meta(r).muscle(m).peaktime_var(t) = nanvar([meta(r).muscle(m).peaktime(t,:)]);
% %             meta(r).muscle(m).peaktime_std(t) = nanstd([meta(r).muscle(m).peaktime(t,:)]);
% %             meta(r).muscle(m).peak_mean(t) = nanmean([meta(r).muscle(m).peak(t,:)]);
% %             meta(r).muscle(m).peak_var(t) = nanvar([meta(r).muscle(m).peak(t,:)]);
% %             meta(r).muscle(m).peak_std(t) = nanstd([meta(r).muscle(m).peak(t,:)]);
            close all;
            figure;maximize;
            plot(trial(n).epochTimeVec, trial(n).meanTrace{1, m},'LineWidth',2);
            hold on;
            plot(trial(n).epochTimeVec, abs(trial(n).meanTrace{1, m}),'LineWidth',2);
            
                title(strcat('e',num2str(meta(r).elec),' - ', Muscle_ID(m),' - trial ', num2str(t)));
                [x,y] = ginput(2);
                if isempty(x)
                    meta(r).muscle(m).onset(t) = NaN;
                    meta(r).muscle(m).peak_onset(t) = NaN;
                    meta(r).muscle(m).peak(t) = NaN;   
                    meta(r).muscle(m).p2pResponse(t) = NaN;
                    meta(r).muscle(m).peakResponse(t) = NaN;
                    meta(r).muscle(m).peaktime(t) = NaN;
                    meta(r).muscle(m).responseRMS(t) = NaN;
                else
                    meta(r).muscle(m).onset(t) = x(1);
                    meta(r).muscle(m).peak_onset(t) = x(2);
                    meta(r).muscle(m).peak(t) = x(2);  
                    meta(r).muscle(m).p2pResponse(t) = peak2peak(trial(n).meanTrace{m}(stimEnd:end));
                    [meta(r).muscle(m).peakResponse(t), mi] = max(abs(trial(n).meanTrace{m}(stimEnd:end)));
                    meta(r).muscle(m).peaktime(t) = trial(n).epochTimeVec(mi+stimEnd-1);
                    meta(r).muscle(m).responseRMS(t) = trial(n).responseRMS(m);
                end

            meta(r).mean(m).base(t) = mean(trial(n).baseline(m,:));
            meta(r).mean(m).rmsbase(t) = rms(trial(n).baseline(m,:));

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


%% Extract Thresholds
for r = 1:length(elecs)
    for i = 1:length(msubset)
        if isempty(find(~isnan(meta(r).muscle(i).onset)))
            meta(r).muscle(i).exist = 0;
            meta(r).muscle(i).threshold = NaN;
            meta(r).muscle(i).threshIdx = NaN;
        else
            meta(r).muscle(i).exist = 1;
        tmp = find(~isnan(meta(r).muscle(i).onset));
        meta(r).muscle(i).threshold = meta(r).StimAmps(tmp(1));
        meta(r).muscle(i).threshIdx = tmp(1);
        end
    end
end
%% Plot Onsets
msub = [1 2 3 4 5 8]; % Remove Ham/TA and SO/MG for LSP02b
% % esub = [1 2 3 4]; %[5 6 7 8] %[9 10 11 12] %[33 34 35];
figure; maximize; tiledlayout(6,4);
for m = msub
    for e = 1:4
      switch e
            case 1
                esub = [1 2 3 4]; %[5 6 7 8] %[9 10 11 12] %[33 34 35];
            case 2
                esub = [5 6 7 8]; %[9 10 11 12] %[33 34 35];
            case 3
                esub = [9 10 11 12]; %[33 34 35];
            case 4
                esub = [33 34 35];
      end
        ax(m) = nexttile; 
        for i = esub
            r = find([response.elec] == i);
            bubblechart(response(r).muscle(m+8).onset*1000,response(r).muscle(m).onset*1000, [(response(r).StimAmps)],'MarkerFaceAlpha',0.20);
            hold on;
        end
        title([muscles(m)]);
        xlabel('Residual Limb (ms)');
        ylabel('Intact Limb (ms)');
        ylim([10 35]);
        xlim([10 35]);
% % ylim('auto');xlim('auto');
        pl = line(xlim, ylim, 'Color', '#D3D3D3','LineStyle','--');
        box off;
        axis(ax(m), 'square')
        bubblesize(ax(m),[0.3 18])
        

        if m == msub(1) && e == 4
            blgd = bubblelegend('Onset: Residual vs Intact');
            blgd.Layout.Tile = 'west';
        elseif m == msub(1) && e==1
            legendCell = cellstr(num2str(esub', 'E%-d'));
            lgd = legend(legendCell, 'Location', 'northoutside', 'orientation', 'horizontal');
        elseif m == msub(end) && e==2
            legendCell = cellstr(num2str(esub', 'E%-d'));
            lgd = legend(legendCell, 'Location', 'southoutside', 'orientation', 'horizontal');
        elseif m == msub(1) && e==3
            legendCell = cellstr(num2str(esub', 'E%-d'));
            lgd = legend(legendCell, 'Location', 'northoutside', 'orientation', 'horizontal');
        elseif m == msub(end) && e==4
            legendCell = cellstr(num2str(esub', 'E%-d'));
            lgd = legend(legendCell, 'Location', 'southoutside', 'orientation', 'horizontal');
        end
    end
        
end
sgtitle('Latency of Response Onset');

%% Plot Onsets - Peak Max time
msub = [1 2 3 4 5 8]; % Remove Ham/TA and SO/MG for LSP02b
% % esub = [1 2 3 4]; %[5 6 7 8] %[9 10 11 12] %[33 34 35];
figure; tiledlayout(6,4); maximize; 
for m = msub
    for e = 1:4
      switch e
            case 1
                esub = [1 2 3 4]; %[5 6 7 8] %[9 10 11 12] %[33 34 35];
            case 2
                esub = [5 6 7 8]; %[9 10 11 12] %[33 34 35];
            case 3
                esub = [9 10 11 12]; %[33 34 35];
            case 4
                esub = [33 34 35];
      end
        ax(m) = nexttile; 
        for i = esub
            r = find([response.elec] == i);
            bubblechart(response(r).muscle(m+8).peaktime*1000,response(r).muscle(m).peaktime*1000, [(response(r).StimAmps)],'MarkerFaceAlpha',0.20);
            hold on;
        end
        title([muscles(m)]);
        xlabel('Residual Limb (ms)');
        ylabel('Intact Limb (ms)');
        ylim([20 40]);
        xlim([20 40]);
% % ylim('auto');xlim('auto');
        pl = line(xlim, ylim, 'Color', '#D3D3D3','LineStyle','--');
        box off;
        axis(ax(m), 'square')
        bubblesize(ax(m),[0.3 18])
        

        if m == msub(1) && e == 4
            blgd = bubblelegend('Onset: Residual vs Intact');
            blgd.Layout.Tile = 'west';
        elseif m == msub(1) && e==1
            legendCell = cellstr(num2str(esub', 'E%-d'));
            lgd = legend(legendCell, 'Location', 'northoutside', 'orientation', 'horizontal');
        elseif m == msub(end) && e==2
            legendCell = cellstr(num2str(esub', 'E%-d'));
            lgd = legend(legendCell, 'Location', 'southoutside', 'orientation', 'horizontal');
        elseif m == msub(1) && e==3
            legendCell = cellstr(num2str(esub', 'E%-d'));
            lgd = legend(legendCell, 'Location', 'northoutside', 'orientation', 'horizontal');
        elseif m == msub(end) && e==4
            legendCell = cellstr(num2str(esub', 'E%-d'));
            lgd = legend(legendCell, 'Location', 'southoutside', 'orientation', 'horizontal');
        end
    end
        
end
sgtitle('Latency of Response Max');

%% Plot Onsets - First Peak
msub = [1 2 3 4 5 8]; % Remove Ham/TA and SO/MG for LSP02b
% % esub = [1 2 3 4]; %[5 6 7 8] %[9 10 11 12] %[33 34 35];
figure; tiledlayout(6,4); maximize; 
for m = msub
    for e = 1:4
      switch e
            case 1
                esub = [1 2 3 4]; %[5 6 7 8] %[9 10 11 12] %[33 34 35];
            case 2
                esub = [5 6 7 8]; %[9 10 11 12] %[33 34 35];
            case 3
                esub = [9 10 11 12]; %[33 34 35];
            case 4
                esub = [33 34 35];
      end
        ax(m) = nexttile; 
        for i = esub
            r = find([response.elec] == i);
            bubblechart(response(r).muscle(m+8).peak_onset*1000,response(r).muscle(m).peak_onset*1000, [(response(r).StimAmps)],'MarkerFaceAlpha',0.20);
            hold on;
        end
        title([muscles(m)]);
        xlabel('Residual Limb (ms)');
        ylabel('Intact Limb (ms)');
        ylim([12 32]);
        xlim([12 32]);
% % ylim('auto');xlim('auto');
        pl = line(xlim, ylim, 'Color', '#D3D3D3','LineStyle','--');
        box off;
        axis(ax(m), 'square')
        bubblesize(ax(m),[0.3 18])
        

        if m == msub(1) && e == 4
            blgd = bubblelegend('Onset: Residual vs Intact');
            blgd.Layout.Tile = 'west';
        elseif m == msub(1) && e==1
            legendCell = cellstr(num2str(esub', 'E%-d'));
            lgd = legend(legendCell, 'Location', 'northoutside', 'orientation', 'horizontal');
        elseif m == msub(end) && e==2
            legendCell = cellstr(num2str(esub', 'E%-d'));
            lgd = legend(legendCell, 'Location', 'southoutside', 'orientation', 'horizontal');
        elseif m == msub(1) && e==3
            legendCell = cellstr(num2str(esub', 'E%-d'));
            lgd = legend(legendCell, 'Location', 'northoutside', 'orientation', 'horizontal');
        elseif m == msub(end) && e==4
            legendCell = cellstr(num2str(esub', 'E%-d'));
            lgd = legend(legendCell, 'Location', 'southoutside', 'orientation', 'horizontal');
        end
    end
        
end
sgtitle('Latency of First Peak');

%Try for one electrode


%% Plot Responses
msub = [1 2 3 4 5 8]; % Remove Ham/TA and SO/MG for LSP02b
% % esub = [1 2 3 4]; %[5 6 7 8] %[9 10 11 12] %[33 34 35];
figure; tiledlayout(6,4); maximize; 
for m = msub
    for e = 1:4
      switch e
            case 1
                esub = [1 2 3 4]; %[5 6 7 8] %[9 10 11 12] %[33 34 35];
            case 2
                esub = [5 6 7 8]; %[9 10 11 12] %[33 34 35];
            case 3
                esub = [9 10 11 12]; %[33 34 35];
            case 4
                esub = [33 34 35];
      end
        ax(m) = nexttile; 
        for i = esub
            r = find([response.elec] == i);
            bubblechart(response(r).muscle(m+8).p2pResponse,response(r).muscle(m).p2pResponse, [(response(r).StimAmps)],'MarkerFaceAlpha',0.20);
            hold on;
            ax_max(i) = max(max(response(r).muscle(m+8).p2pResponse),max(response(r).muscle(m).p2pResponse));
        end
        title([muscles(m)]);
        xlabel('Residual Limb (mV)');
        ylabel('Intact Limb (mV)');
        ylim([0 max(ax_max)]);
        xlim([0 max(ax_max)]);
% % ylim('auto');xlim('auto');
        pl = line(xlim, ylim, 'Color', '#D3D3D3','LineStyle','--');
        box off;
        axis(ax(m), 'square')
        bubblesize(ax(m),[0.3 18])
        

        if m == msub(1) && e == 4
            blgd = bubblelegend('Onset: Residual vs Intact');
            blgd.Layout.Tile = 'west';
        elseif m == msub(1) && e==1
            legendCell = cellstr(num2str(esub', 'E%-d'));
            lgd = legend(legendCell, 'Location', 'northoutside', 'orientation', 'horizontal');
        elseif m == msub(end) && e==2
            legendCell = cellstr(num2str(esub', 'E%-d'));
            lgd = legend(legendCell, 'Location', 'southoutside', 'orientation', 'horizontal');
        elseif m == msub(1) && e==3
            legendCell = cellstr(num2str(esub', 'E%-d'));
            lgd = legend(legendCell, 'Location', 'northoutside', 'orientation', 'horizontal');
        elseif m == msub(end) && e==4
            legendCell = cellstr(num2str(esub', 'E%-d'));
            lgd = legend(legendCell, 'Location', 'southoutside', 'orientation', 'horizontal');
        end
    end
        
end
sgtitle('Response Amplitude: P2P');

%Try for one electrode

%% Plot ALL Onsets
%% Plot Onsets and Thresholds
musclesL = {'VM', 'RF', 'VL', 'BF', 'ST', 'TA', 'MG', 'LG',};
% % msubset = [10 12 16 14];
% msubset = [4 12 5 13];
% % elecs = [4 7 9 11 18 13 20 16 23 26];
elecs = [1 2 3 4 5 6 7 8 9 10 11 12 33 34 35]; %LSP02b
msubset = [9:16 1:8];
% C = linspecer(length(msubset));
cmap = [0.0784    0.1529    0.3020; 0.0784    0.1529    0.3020; 0.0784    0.1529    0.3020; 0.9176    0.4235    0.1255; 0.9176    0.4235    0.1255; 0.9529    0.6157    0.4510; 0.3843    0.4745    0.6157; 0.3843    0.4745    0.6157];
C = [cmap; cmap];

for i = 1:length(msubset)
    m = msubset(i);
    onsets = [];
    pkonsets = [];
    mxonsets = [];
    for r = 1:length(elecs)
        e = find([response.elec] == elecs(r));
        tmpOns = [response(e).muscle(m).onset];
        onsets = [onsets tmpOns];
        pkOns = [response(e).muscle(m).peak_onset];
        pkonsets = [pkonsets pkOns];
        mxOns = [response(e).muscle(m).peaktime];
        mxonsets = [mxonsets mxOns];
    end
    onserror(i) = nanstd( onsets ) / sqrt( length(onsets) );
    onset_mean(i) = nanmean( onsets );
    pkonserror(i) = nanstd( pkonsets ) / sqrt( length(pkonsets) );
    pkonset_mean(i) = nanmean( pkonsets );
    mxonserror(i) = nanstd( mxonsets ) / sqrt( length(mxonsets) );
    mxonset_mean(i) = nanmean( mxonsets );
end

%% Plot
hS6 = figure;%maximize;
tiledlayout(3,1);
% % nexttile
ax = nexttile;
superbar(onset_mean*1000','E', onserror*1000','BarFaceColor', C);%permute(C, [3 1 2]));
ylabel('Onset (ms)');
% tmp = [' ', mlab, ' '];
% xticklabels(tmp);
xticks([1:length(msubset)]);
% xticklabels(musclesL(msubset-8));
xticklabels(Muscle_ID(msubset));
xtickangle(ax,45)
% legend(arrayfun(@num2str,elecs,'UniformOutput',false), 'Location', 'eastoutside', 'Orientation', 'vertical');
hold on;

ax2 = nexttile;
superbar(pkonset_mean*1000','E', pkonserror*1000','BarFaceColor', C);%permute(C, [3 1 2]));
ylabel('First Peak (ms)');
% tmp = [' ', mlab, ' '];
% xticklabels(tmp);
xticks([1:length(msubset)]);
% xticklabels(musclesL(msubset-8));
xticklabels(Muscle_ID(msubset));
xtickangle(ax2,45)
% legend(arrayfun(@num2str,elecs,'UniformOutput',false), 'Location', 'eastoutside', 'Orientation', 'vertical');
hold on;

ax3 = nexttile;
superbar(mxonset_mean*1000','E', mxonserror*1000','BarFaceColor', C);%permute(C, [3 1 2]));
ylabel('Response Max (ms)');
% tmp = [' ', mlab, ' '];
% xticklabels(tmp);
xticks([1:length(msubset)]);
% xticklabels(musclesL(msubset-8));
xticklabels(Muscle_ID(msubset));
xtickangle(ax3,45)
% legend(arrayfun(@num2str,elecs,'UniformOutput',false), 'Location', 'eastoutside', 'Orientation', 'vertical');
hold on;
sgtitle([subjectName ' Response Latencies']);