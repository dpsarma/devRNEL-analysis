
%VL-RF: '14274d' BF-ST: 'ea6c20'  LG-MG: '62799d' TA: 'f39d73'}
subjectName = 'LSP05';
switch subjectName
    case 'LSP05'
        Muscle_ID = {'Right VM', 'Right RF', 'Right VL', 'Right BF', 'Right ST', ...
        'Right TA', 'Right SO', 'Right LG','Left VM', 'Left RF','Left VL',...
        'Left BF', 'Left ST', 'Left TA', 'Left MG', 'Left LG'};
        load('D:\FigRescources\UH3\HotelCali\LSP05\LSP05_response_meta.mat')
        elecs = [1 2 3 4 7 8 9 10 11 12 13 14 15 16 18 19 20 23 24 25 26 27 28 29 30 31];
        response5 = response;
        muscles = {'VM', 'RF', 'VL', 'BF', 'ST', 'SO/TA', 'MG', 'LG',};

        
    
    case 'LSP02b'
        Muscle_ID = {'Right VM', 'Right RF', 'Right VL', 'Right BF', 'Right ST', ...
        'Right TA', 'Right SO', 'Right LG','Left VM','Left RF','Left VL',...
        'Left BF', 'Left ST', 'Left Ham', 'Left MG', 'Left LG'};
        load('D:\FigRescources\UH3\HotelCali\LSP02b\LSP02b_response_meta.mat');
        elecs = [1 2 3 4 5 6 7 8 9 10 11 12 33 34];
        response2b = response;
        muscles = {'VM', 'RF', 'VL', 'BF', 'ST', 'Ham/TA', 'MG', 'LG',};
end
% cmap = [0.0784    0.1529    0.3020; 0.0784    0.1529    0.3020; 0.0784    0.1529    0.3020; 0.9176    0.4235    0.1255; 0.9176    0.4235    0.1255; 0.9529    0.6157    0.4510; 0.3843    0.4745    0.6157; 0.3843    0.4745    0.6157];

%% PCA for AUC
for iE = 1:length(response)

for iE = 1:length(response)
    for iM = 1 : length(Muscle_ID)
        xinfo = response(iE).muscle(iM).p2pResponse;
        ampinfo = [ 0 response(iE).StimAmps];
        aucData(iE,iM) = sum(xinfo.*diff(ampinfo));
        normAuc(iE,iM) = aucData(iE,iM)/max(aucData(:,iM));
    end
    
end

%% Extract Threshold Amps
for r = 1:length(response)
    for m = 1:length(Muscle_ID)
        threshes(r,m) = response(r).muscle(m).threshold;
        max_response(r,m) = max([response(r).muscle(m).p2pResponse]);
        active(r,m) = response(r).muscle(m).exist;
    end
end

C = linspecer(length(Muscle_ID));


for i = 1:length(elecs)
    for m = 1:8
    T_ratio(i,m) = threshes(i,m+8)/threshes(i,m);
    R_ratio(i,m) = max_response(i,m+8)/max_response(i,m);
    AUC_ratio(i,m) = (aucData(i,m+8) - aucData(i,m))/(aucData(i,m+8) + aucData(i,m));
    norm_AUCrat(i,m) = (normAuc(i,m+8) - normAuc(i,m))/(normAuc(i,m+8) + normAuc(i,m));
    end
end

%% Plot Activation

% % x = sum(active)/length(active)
% % figure; bar(x)
% % xticks(1:16)
% % xticklabels(Muscle_ID)

%% Plot Scatter


% % % threshes(threshes == 10000) = NaN;
% % figure;tiledlayout('flow');
% % nexttile;
% % for i = 1:length(Muscle_ID)/2
% % scatter(elecs,threshes(:,i)',[],C(i,:));
% % % %     p = polyfit(elecs,threshes(:,i)',2);
% % % %     h = refcurve(p);
% % % %     h.Color = C(i,:);
% % hold on;
% % end
% % legend(Muscle_ID(1:8),'location','eastoutside');
% % nexttile;
% % for i = 9:length(Muscle_ID)
% % scatter(elecs,threshes(:,i)',[],C(i,:));
% % hold on;
% % % %     p = polyfit(elecs,threshes(:,i)',2);
% % % %     h = refcurve(p);
% % % %     h.Color = C(i,:);
% % end
% % legend(Muscle_ID(9:16),'location','eastoutside');
%% Threshes against each other
% figure; 
% type = {'KF', 'KF', 'KF', 'KE', 'KE', 'AE', 'AF', 'AF'};
% for i=1:length(Muscle_ID)/2
%     bubbleplot(threshes(:,i+8),threshes(:,i),[],T_ratio(:,i), CC, muscles);
%     hold on;
%     xlabel('Residual Limb');
%     ylabel('Intact Limb');
% end
% legend(muscles);

% figure; 
% for i=1:length(elecs)
%     scatter([1:8], T_ratio(i,:));
%     hold on;
% end

% % figure;tiledlayout('flow');maximize;
type = {'KF', 'KF', 'KF', 'KE', 'KE', 'AE', 'AF', 'AF'};
CC = [0.0784    0.1529    0.3020; 0.0784    0.1529    0.3020; 0.0784    0.1529    0.3020; 0.9176    0.4235    0.1255; 0.9176    0.4235    0.1255; 0.9529    0.6157    0.4510; 0.3843    0.4745    0.6157; 0.3843    0.4745    0.6157];
colormap(CC);
C3 = linspecer(length(Muscle_ID)/2);

nexttile; 
for i=1:length(elecs)
% %     bubbleplot(threshes(i,1:8),threshes(i,9:16),[],AUC_ratio(i,:), type,muscles);
    bubblechart(threshes(i,4:6), threshes(i,12:14), normAuc(i,4:6), 'MarkerFaceAlpha',0.20);
    
    hold on;
    ylabel('Residual Limb');
    xlabel('Intact Limb');
    ylim([0,6500]);
    xlim([0,6500]);
end
legend(muscles([4:6]),'Location','eastoutside');
bubblelegend('Activation: Residual vs Intact','Location','northeastoutside')
pl = line(xlim, ylim, 'Color', '#D3D3D3','LineStyle','--');
box off;
title([subjectName ' - Flexors']);
%%
nexttile; 
for i=1:length(elecs)
% %     bubbleplot(threshes(i,1:8),threshes(i,9:16),[],AUC_ratio(i,:), type,muscles);
    bubblechart(threshes(i,[1:3 7:8]), threshes(i,[9:11 15:16]), normAuc(i,[1:3 7:8]), 'MarkerFaceAlpha',0.20);
    
    hold on;
    ylabel('Residual Limb');
    xlabel('Intact Limb');
    ylim([0,6500]);
    xlim([0,6500]);
end
legend(muscles([1:3 7:8]),'Location','eastoutside');
bubblelegend('Activation: Residual vs Intact','Location','northeastoutside')
pl = line([0 6000],[0 6000], 'Color', '#D3D3D3','LineStyle','--');
box off;
title([subjectName ' - Extensors']);
%% Activation Ratio
