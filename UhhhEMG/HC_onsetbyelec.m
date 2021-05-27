%% Plot Onsets and Thresholds
musclesL = {'VM', 'RF', 'VL', 'BF', 'ST', 'TA', 'MG', 'LG',};
% % msubset = [10 12 16 14];
% msubset = [4 12 5 13];
% % elecs = [4 7 9 11 18 13 20 16 23 26];
% % elecs = [1 2 3 4 5 6 7 8 9 10 11 12 33 34 35]; %LSP02b
% elecs = [9:16]
elecs = [response.elec];
msubset = [9:16 1:8];
% C = linspecer(length(msubset));
cmap = [0.0784    0.1529    0.3020; 0.0784    0.1529    0.3020; 0.0784    0.1529    0.3020; 0.9176    0.4235    0.1255; 0.9176    0.4235    0.1255; 0.9529    0.6157    0.4510; 0.3843    0.4745    0.6157; 0.3843    0.4745    0.6157];
C = [cmap; cmap];

figure; tiledlayout(length(elecs),2, 'TileSpacing','Compact');

for r = 1:length(elecs)
    e = find([response.elec] == elecs(r));
    mxonsets = [];
    for i = 1:length(msubset)
        m = msubset(i);
        
        mxOns = [response(e).muscle(m).peaktime];
        mxonsets = [mxOns];
        mxonserror(i) = nanstd( mxonsets ) / sqrt( length(mxonsets) );
        mxonset_mean(i) = nanmean( mxonsets );
        
    end
        ax(r) = nexttile(1+2*(r-1));
        superbar(mxonset_mean*1000','E', mxonserror*1000','BarFaceColor', C);%permute(C, [3 1 2]));
        ylabel(['e' num2str(elecs(r))]);
        yticks([0 20 50]);
        if r < length(elecs)
        set(gca,'XColor','none')   
        set(gca,'XTickLabel',[]);
        end
        box off
end
% tmp = [' ', mlab, ' '];
% xticklabels(tmp);
xticks([1:length(msubset)]);
% xticklabels(musclesL(msubset-8));
xticklabels(Muscle_ID(msubset));
xtickangle(ax(r),45)
% legend(arrayfun(@num2str,elecs,'UniformOutput',false), 'Location', 'eastoutside', 'Orientation', 'vertical');

linkaxes([ax(:)],'xy');
sgtitle([subjectName ' Response Latencies']);


% %% Plot Scatter Onsets
% musclesL = {'VM', 'RF', 'VL', 'BF', 'ST', 'TA', 'MG', 'LG',};
% % % msubset = [10 12 16 14];
% % msubset = [4 12 5 13];
% % % elecs = [4 7 9 11 18 13 20 16 23 26];
% % % elecs = [1 2 3 4 5 6 7 8 9 10 11 12 33 34 35]; %LSP02b
% % elecs = [9:16]
% elecs = [response.elec];
% msubset = [9:16 1:8];
% % C = linspecer(length(msubset));
% cmap = [0.0784    0.1529    0.3020; 0.0784    0.1529    0.3020; 0.0784    0.1529    0.3020; 0.9176    0.4235    0.1255; 0.9176    0.4235    0.1255; 0.9529    0.6157    0.4510; 0.3843    0.4745    0.6157; 0.3843    0.4745    0.6157];
% C = [cmap; cmap];
% 
% figure; %tiledlayout(length(elecs),1, 'TileSpacing','Compact');
% 
% for r = 1:length(elecs)
%     e = find([response.elec] == elecs(r));
%     mxonsets = [];
%     for i = 1:length(msubset)
%         m = msubset(i);
%         
%         mxOns = [response(e).muscle(m).peaktime];
%         mxonsets = [mxOns];
%         mxonserror(i) = nanstd( mxonsets ) / sqrt( length(mxonsets) );
%         mxonset_mean(i) = nanmean( mxonsets );
%         
%     end
%         %ax(r) = nexttile;
%         scatter([1:16],mxonset_mean*1000'); hold on; 
% % %        xx = violinplot(mxonset_mean*1000', Muscle_ID); hold on;
% %         superbar(mxonset_mean*1000','E', mxonserror*1000','BarFaceColor', C);%permute(C, [3 1 2]));
% %         ylabel(['e' num2str(elecs(r))]);
% % %         if r < length(elecs)
% % %         set(gca,'XColor','none')   
% % %         set(gca,'XTickLabel',[]);
% % %         end
%         box off
% end
% % tmp = [' ', mlab, ' '];
% % xticklabels(tmp);
% xticks([1:length(msubset)]);
% % xticklabels(musclesL(msubset-8));
% xticklabels(Muscle_ID(msubset));
% % % xtickangle(ax(r),45);
% xtickangle(45);
% % legend(arrayfun(@num2str,elecs,'UniformOutput',false), 'Location', 'eastoutside', 'Orientation', 'vertical');
% % % 
% % % linkaxes([ax(:)],'xy');
% % % sgtitle([subjectName ' Response Latencies']);

%% Plot Response Amp
musclesL = {'VM', 'RF', 'VL', 'BF', 'ST', 'TA', 'MG', 'LG',};
Muscle_ID = {'Intact VM', 'Intact RF', 'Intact VL', 'Intact BF', 'Intact ST', ...
'Intact TA', 'Intact SO', 'Intact LG','Residual VM', 'Residual RF','Residual VL',...
'Residual BF', 'Residual ST', 'Residual Ham', 'Residual MG', 'Residual LG'}; %Ham for 2b, TA for 05
% % msubset = [10 12 16 14];
% msubset = [4 12 5 13];
% % elecs = [4 7 9 11 18 13 20 16 23 26];
% % elecs = [1 2 3 4 5 6 7 8 9 10 11 12 33 34 35]; %LSP02b
% elecs = [9:16]
elecs = [response.elec];
msubset = [9:16 1:8];
% C = linspecer(length(msubset));
cmap = [0.0784    0.1529    0.3020; 0.0784    0.1529    0.3020; 0.0784    0.1529    0.3020; 0.9176    0.4235    0.1255; 0.9176    0.4235    0.1255; 0.9529    0.6157    0.4510; 0.3843    0.4745    0.6157; 0.3843    0.4745    0.6157];
C = [cmap; cmap];


% % elecs = [9 16]

% % elecs = [1 12 33]
figure; tiledlayout(length(elecs),1, 'TileSpacing','Compact');
for r = 1:length(elecs)
    e = find([response.elec] == elecs(r));
    mxonsets = [];
    for i = 1:length(msubset)
        m = msubset(i);
        
        mxOns = (response(e).muscle(m).p2pResponse - response(e).mean(m).rmsbase) ./ response(e).mean(m).rmsbase *100; %[response(e).muscle(m).p2pResponse];%
        mxonsets = mxOns;
        mxonserror(i) = nanstd( mxonsets ) / sqrt( length(mxonsets) );
        mxonset_mean(i) = nanmean( mxonsets );
        
    end
% %         ax(r) = nexttile(2*(r));
        ax(r) = nexttile;
        superbar(mxonset_mean','E', mxonserror','BarFaceColor', C);%permute(C, [3 1 2]));
       
        ylabel({['e' num2str(elecs(r))]; '(/deltaEMG %)'});
        if r < length(elecs)
        set(gca,'XColor','none')   
        set(gca,'XTickLabel',[]);
        end
        box off
end
% tmp = [' ', mlab, ' '];
% xticklabels(tmp);
xticks([1:length(msubset)]);
% xticklabels(musclesL(msubset-8));
xticklabels(Muscle_ID(msubset));
xtickangle(ax(r),45)
% legend(arrayfun(@num2str,elecs,'UniformOutput',false), 'Location', 'eastoutside', 'Orientation', 'vertical');

linkaxes([ax(:)],'xy');
sgtitle([subjectName ' Response Intensity']);
% % sgtitle([{subjectName; 'Onset of Peak (ms)       

%% Plot Response Amp BOXes
musclesL = {'VM', 'RF', 'VL', 'BF', 'ST', 'TA', 'MG', 'LG',};
% % Muscle_ID = {'Intact VM', 'Intact RF', 'Intact VL', 'Intact BF', 'Intact ST', ...
% % 'Intact TA', 'Intact SO', 'Intact LG','Residual VM', 'Residual RF','Residual VL',...
% % 'Residual BF', 'Residual ST', 'Residual TA', 'Residual MG', 'Residual LG'}; %Ham for 2b, TA for 05
% % msubset = [10 12 16 14];
% msubset = [4 12 5 13];
% % elecs = [4 7 9 11 18 13 20 16 23 26];
% % elecs = [1 2 3 4 5 6 7 8 9 10 11 12 33 34 35]; %LSP02b
% elecs = [9:16]
elecs = [response.elec];
msubset = [9:16 1:8];
% C = linspecer(length(msubset));
cmap = [0.0784    0.1529    0.3020; 0.0784    0.1529    0.3020; 0.0784    0.1529    0.3020; 0.9176    0.4235    0.1255; 0.9176    0.4235    0.1255; 0.9529    0.6157    0.4510; 0.3843    0.4745    0.6157; 0.3843    0.4745    0.6157];
C = [cmap; cmap];


% % elecs = [9 16]

elecs = [1 12 33]
figure; tiledlayout(length(elecs),1, 'TileSpacing','Compact');
for r = 1:length(elecs)
    e = find([response.elec] == elecs(r));
    Rx = [];

    for i = 1:length(msubset)
        m = msubset(i);
        
        norm_peak = (response(e).muscle(m).p2pResponse - response(e).mean(m).rmsbase) ./ response(e).mean(m).rmsbase *100; %[response(e).muscle(m).p2pResponse];%
        Rx(i,:) = norm_peak;
       
    end
        ax(r) = nexttile;
        boxplot(Rx', 'BoxStyle','filled', 'Colors', C); hold on;
        ylabel({['e' num2str(elecs(r))]; '(\deltaEMG %)'});
        if r < length(elecs)
        set(gca,'XColor','none')   
        set(gca,'XTickLabel',[]);
        end
        box off
end
% tmp = [' ', mlab, ' '];
% xticklabels(tmp);
xticks([1:length(msubset)]);
% xticklabels(musclesL(msubset-8));
xticklabels(Muscle_ID(msubset));
xtickangle(ax(r),45)
% legend(arrayfun(@num2str,elecs,'UniformOutput',false), 'Location', 'eastoutside', 'Orientation', 'vertical');

linkaxes([ax(:)],'xy'); %ylim([0 14000]);
sgtitle(['Subject 1, Response Intensity']);%subjectName
% % sgtitle([{subjectName; 'Onset of Peak (ms)     Intensity (% Baseline)'}]);

%% Try again for All Onsets:
elecs = [response.elec];
msubset = [9:16];
% C = linspecer(length(msubset));
cmap = [0.0784    0.1529    0.3020; 0.0784    0.1529    0.3020; 0.0784    0.1529    0.3020; 0.9176    0.4235    0.1255; 0.9176    0.4235    0.1255; 0.9529    0.6157    0.4510; 0.3843    0.4745    0.6157; 0.3843    0.4745    0.6157];
C = [cmap; cmap];
onsets = [];

for i = 1:length(msubset)
    m = msubset(i);
    ons = [];
    for r = 1:length(elecs)
        e = find([response.elec] == elecs(r));
        tmpOns = [response(e).muscle(m).peaktime]; %.onset
        ons = [ons tmpOns];
    end
    onsets(i,:) = ons;
    mxOn(i,:) = nanmean(ons);
    mxOn_err(i,:) = nanstd( ons ) / sqrt( length(ons) );
end
figure;
xx = boxplot(onsets'*1000, 'BoxStyle','filled', 'Colors', C,'orientation','horizontal','Symbol', '.'); hold on;
%         ylabel({['e' num2str(elecs(r))]; '(\deltaEMG %)'});
        xlim([10 40]);
% %         if r < length(elecs)
% %         set(gca,'XColor','none')   
% %         set(gca,'XTickLabel',[]);
% %         end
xlabel('Response Onset (ms)');
yticks([1:length(msubset)]);

yticklabels(Muscle_ID(msubset));
sgtitle(['Subject 1, Response Latencies']);%subjectName

figure; 
superbar(mxOn*1000','E', mxOn_err*1000','BarFaceColor', C, 'Orientation', 'h');
xlim([0 40]);
xlabel('Onset of Peak (ms)');
% tmp = [' ', mlab, ' '];
% xticklabels(tmp);
yticks([1:length(msubset)]);
% xticklabels(musclesL(msubset-8));
yticklabels(Muscle_ID(msubset));
% % ytickangle(ax3,45)
% legend(arrayfun(@num2str,elecs,'UniformOutput',false), 'Location', 'eastoutside', 'Orientation', 'vertical');
hold on;
sgtitle(['Subject 1, Response Latencies']);%subjectName
