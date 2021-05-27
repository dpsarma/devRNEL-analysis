%% Heatmap at SenseThreshold and 2X Sense Threshold
% Using response/meta from Hotel Cali NewOnset Calc

load('G:\My Drive\UH3_metagen\LSP05\LSP05_response_meta_9-16.mat')
response = meta
subjectName = 'LSP05'; setName = 'RCs';
Muscle_ID = {'Intact VM', 'Intact RF', 'Intact VL', 'Intact BF', 'Intact ST', ...
'Intact TA', 'Intact SO', 'Intact LG','Residual VM', 'Residual RF','Residual VL',...
'Residual BF', 'Residual ST', 'Residual TA', 'Residual MG', 'Residual LG'};
elecs = [response.elec];
msubset = [1:16];

sensethresh = [2000]; %e9 e6

%% Extract Thresholds
for r = 1:length(elecs)
    for i = 1:length(msubset)
        if isempty(find(~isnan(meta(r).muscle(i).peaktime), 1)) %.onset
            meta(r).muscle(i).exist = 0;
            meta(r).muscle(i).threshold = NaN;
            meta(r).muscle(i).threshIdx = NaN;
        else
            meta(r).muscle(i).exist = 1;
        tmp = find(~isnan(meta(r).muscle(i).peaktime)); %.onset
        meta(r).muscle(i).threshold = meta(r).StimAmps(tmp(1));
        meta(r).muscle(i).threshIdx = tmp(1);
        end
    end
end

%% Extract Intensity values. 
TwoXThresh = 2*sensethresh;
e = find(elecs == 16);
resp2TIdx = find(meta(e).StimAmps == TwoXThresh);
respTIdx = find(meta(e).StimAmps == sensethresh);

for m = 1:length(msubset)
    Thresh_peak(m) = meta(e).muscle(m).p2pResponse(respTIdx);
    TwoX_peak(m) = meta(e).muscle(m).p2pResponse(resp2TIdx);
    Thresh_norm(m) = (meta(e).muscle(m).p2pResponse(respTIdx) - meta(e).mean(m).rmsbase(respTIdx)) / (meta(e).mean(m).rmsbase(respTIdx)) *100;
    TwoX_norm(m) = (meta(e).muscle(m).p2pResponse(resp2TIdx) - meta(e).mean(m).rmsbase(resp2TIdx)) / (meta(e).mean(m).rmsbase(resp2TIdx)) *100;
    TwoXvsThresh(m) = (meta(e).muscle(m).p2pResponse(resp2TIdx) - meta(e).muscle(m).p2pResponse(respTIdx)) / (meta(e).muscle(m).p2pResponse(respTIdx))*100;
    ratio_TwoTvsT(m) = (meta(e).muscle(m).p2pResponse(resp2TIdx) / meta(e).muscle(m).p2pResponse(respTIdx));
end
cmap = magma; colormap(cmap);
tiledlayout('flow');

% % nexttile;
% % imagesc(Thresh_peak); caxis([0 inf]);
% % colorbar
% % title('Peak at Thresh'); xticks(1:16); xticklabels(Muscle_ID); xtickangle(45)
% % 
nexttile;
imagesc(TwoX_peak); caxis([0 inf]);
colorbar
title('Peak at 2X'); xticks(1:16); xticklabels(Muscle_ID); xtickangle(45)

nexttile;
imagesc(Thresh_norm); caxis([0 inf]);
colorbar
title('Thresh: % Change of Peak to Base'); xticks(1:16); xticklabels(Muscle_ID); xtickangle(45)
nexttile; 
imagesc(TwoX_norm); caxis([0 inf]);
colorbar
title('2X: % Change of Peak to Base'); xticks(1:16); xticklabels(Muscle_ID); xtickangle(45)
nexttile;
imagesc(TwoXvsThresh); caxis([0 inf]);
colorbar
title('% Change of 2X to Threshold'); xticks(1:16); xticklabels(Muscle_ID); xtickangle(45)
% % nexttile;
% % imagesc(TwoXvsThresh); caxis([0 inf]);
% % colorbar
% % title('Ratio of 2X to Threshold'); xticks(1:16); xticklabels(Muscle_ID); xtickangle(45)


