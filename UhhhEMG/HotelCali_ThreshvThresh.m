
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
        sensethresh = [4000 4000 3500 3000 3000 2000 2000 2500 2500 2000 ...
    2000 2000 2000 2000 3000 3000 2000 2000 2000 2000 1000 1000 2000 1000 1000 2000]; %LSP05


        
    
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

for iM = 1 : length(Muscle_ID)
    for iE = 1:length(response)
        resp(iE) = nanmax(response(iE).muscle(iM).p2pResponse);
    end
    maxResp(iM) =  max(resp);
end
    

for iE = 1:length(response)
    for iM = 1 : length(Muscle_ID)
        xinfo = response(iE).muscle(iM).p2pResponse/maxResp(iM);
        ampinfo = [ 0 response(iE).StimAmps];
        aucData(iE,iM) = nansum(xinfo.*diff(ampinfo));
        normAuc(iE,iM) = aucData(iE,iM)/max(aucData(:,iM));
    end
    
end

%% Extract Threshold Amps
for r = 1:length(response)
    for m = 1:length(Muscle_ID)
        threshes(r,m) = response(r).muscle(m).threshold;
        max_response(r,m) = nanmax([response(r).muscle(m).p2pResponse]);
        active(r,m) = response(r).muscle(m).exist;
    end
end

C = linspecer(length(Muscle_ID));
aucData(isnan(threshes)) =  0;

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

% % figure;maximize; t = tiledlayout(2,2,'TileSpacing','Compact'); %tiledlayout('flow');
type = {'KF', 'KF', 'KF', 'KE', 'KE', 'AE', 'AF', 'AF'};
CC = [0.0784    0.1529    0.3020; 0.0784    0.1529    0.3020; 0.0784    0.1529    0.3020; 0.9176    0.4235    0.1255; 0.9176    0.4235    0.1255; 0.9529    0.6157    0.4510; 0.3843    0.4745    0.6157; 0.3843    0.4745    0.6157];
colormap(CC);
C3 = linspecer(length(Muscle_ID)/2);

ax1 = nexttile; 
for i=1:length(elecs)
% %     bubbleplot(threshes(i,1:8),threshes(i,9:16),[],AUC_ratio(i,:), type,muscles);
    bubblechart(threshes(i,4:6), threshes(i,12:14), AUC_ratio(i,4:6), 'MarkerFaceAlpha',0.20);
    
    hold on;
    ylabel('Residual Limb');
    xlabel('Intact Limb');
    ylim([0,6500]);
    xlim([0,6500]);
end


pl = line(xlim, ylim, 'Color', '#D3D3D3','LineStyle','--');
box off;
bubblelim(ax1,[-1 1])
axis(ax1, 'square')
title([subjectName ' - Flexors']);

blgd = bubblelegend('Activation: Residual vs Intact');
lgd = legend(muscles([4:6]));
blgd.Layout.Tile = 'west';
lgd.Layout.Tile = 'west';
%%
ax2 = nexttile; 
for i=1:length(elecs)
% %     bubbleplot(threshes(i,1:8),threshes(i,9:16),[],AUC_ratio(i,:), type,muscles);
    bubblechart(threshes(i,[1:3 7:8]), threshes(i,[9:11 15:16]), AUC_ratio(i,[1:3 7:8]), 'MarkerFaceAlpha',0.20);
    
    hold on;
    ylabel('Residual Limb');
    xlabel('Intact Limb');
    ylim([0,6500]);
    xlim([0,6500]);
end

bubblelim(ax2,[-1 1])
pl = line([0 6000],[0 6000], 'Color', '#D3D3D3','LineStyle','--');
box off;
axis(ax2, 'square')
title([subjectName ' - Extensors']);

blgd = bubblelegend('Activation: Residual vs Intact');
lgd = legend(muscles([1:3 7:8]));
blgd.Layout.Tile = 'east';
lgd.Layout.Tile = 'east';
%% Violin Plot By threshold
figure;
tiledlayout(3,1)
nexttile;
xx = violinplot(threshes(:,[1 9 2 10 3 11]), Muscle_ID([1 9 2 10 3 11]));
C4 = linspecer(length(xx));
for i = length(xx)
    xx(i).ViolinColor = C4(i,:);
end
ylabel('Stim Threshold');
clear xx C4;
box off;

nexttile;
xx = violinplot(threshes(:,[4 12 5 13]), Muscle_ID([4 12 5 13]));
C4 = linspecer(length(xx));
for i = length(xx)
    xx(i).ViolinColor = C4(i,:);
end
ylabel('Stim Threshold');
clear xx C4;
box off;

nexttile;
xx = violinplot(threshes(:,[6 14 7 15 8 16]), Muscle_ID([6 14 7 15 8 16]));
C4 = linspecer(length(xx));
for i = length(xx)
    xx(i).ViolinColor = C4(i,:);
end
ylabel('Stim Threshold');
clear C4 xx;
sgtitle(subjectName)
box off;

% % for i = [1 9 2 10 3 11 4 12 5 13 6 14 7 16 8 16]
% %     xx(i).ViolinColor = C(i,:);
% % end

%% Bilaterality
for iE = 1:length(elecs)
    n_stim(iE) = length(response(iE).StimAmps);
    for iM = [1 2 3 4 5 8] %1:length(Muscle_ID)/2
        if response(iE).muscle(iM).exist && isnan(response(iE).muscle(iM).threshIdx)
            response(iE).muscle(iM).exist = 0
        elseif response(iE).muscle(iM+8).exist && isnan(response(iE).muscle(iM+8).threshIdx)
            response(iE).muscle(iM+8).exist = 0
        end
        
        if ~response(iE).muscle(iM).exist && ~response(iE).muscle(iM+8).exist
            n_inactive(iE,iM) = length(response(iE).StimAmps);
            n_intact(iE,iM) = 0;
            n_residual(iE,iM) = 0;
            n_bilateral(iE,iM) = 0;
            disp(['Neither Exist - E' num2str(iE) 'M' num2str(iM) '-' num2str(n_inactive(iE,iM)+n_intact(iE,iM)+n_residual(iE,iM)+n_bilateral(iE,iM))]);
        elseif response(iE).muscle(iM).exist && ~response(iE).muscle(iM+8).exist
            n_inactive(iE,iM) = response(iE).muscle(iM).threshIdx-1;
            n_intact(iE,iM) = length(response(iE).StimAmps) - response(iE).muscle(iM).threshIdx +1;
            n_residual(iE,iM) = 0;
            n_bilateral(iE,iM) = 0;
            disp(['Intact Only - E' num2str(iE) 'M' num2str(iM) '-' num2str(n_inactive(iE,iM)+n_intact(iE,iM)+n_residual(iE,iM)+n_bilateral(iE,iM))]);
        elseif ~response(iE).muscle(iM).exist && response(iE).muscle(iM+8).exist
            n_inactive(iE,iM) = response(iE).muscle(iM+8).threshIdx-1;
            n_intact(iE,iM) = 0;
            n_residual(iE,iM) = length(response(iE).StimAmps) - response(iE).muscle(iM+8).threshIdx +1;
            n_bilateral(iE,iM) = 0;
            disp(['Residual Only - E' num2str(iE) 'M' num2str(iM) '-' num2str(n_inactive(iE,iM)+n_intact(iE,iM)+n_residual(iE,iM)+n_bilateral(iE,iM))]);
        elseif response(iE).muscle(iM).exist && response(iE).muscle(iM+8).exist
            if response(iE).muscle(iM).threshIdx > response(iE).muscle(iM+8).threshIdx 
                n_inactive(iE,iM) = response(iE).muscle(iM+8).threshIdx-1;
                n_intact(iE,iM) = 0;
                n_residual(iE,iM) = response(iE).muscle(iM).threshIdx - response(iE).muscle(iM+8).threshIdx;
                n_bilateral(iE,iM) = length(response(iE).StimAmps) - response(iE).muscle(iM).threshIdx + 1;
                disp(['Residual First - E' num2str(iE) 'M' num2str(iM) '-' num2str(n_inactive(iE,iM)+n_intact(iE,iM)+n_residual(iE,iM)+n_bilateral(iE,iM))]);
            elseif response(iE).muscle(iM).threshIdx <= response(iE).muscle(iM+8).threshIdx
                n_inactive(iE,iM) = response(iE).muscle(iM).threshIdx-1;
                n_intact(iE,iM) = response(iE).muscle(iM+8).threshIdx - response(iE).muscle(iM).threshIdx;
                n_residual(iE,iM) = 0;
                n_bilateral(iE,iM) = length(response(iE).StimAmps) - response(iE).muscle(iM+8).threshIdx + 1;
                disp(['Intact First - E' num2str(iE) 'M' num2str(iM) '-' num2str(n_inactive(iE,iM)+n_intact(iE,iM)+n_residual(iE,iM)+n_bilateral(iE,iM))]);
            end
        else
            disp('Fail');
        end
    end
end
%% Plot Bilaterality
colormap('magma');
% leglabs = {'Inactive', 'Bilateral', 'Residual Only', 'Intact Only'};
leglabs = {'Residual Only', 'Intact Only', 'Bilateral', };
% p_state = [(sum(n_inactive)/sum(n_stim))' (sum(n_bilateral)/sum(n_stim))'  (sum(n_residual)/sum(n_stim))'  (sum(n_intact)/sum(n_stim))'];
p_state = [(sum(n_residual)./(sum(n_stim)-sum(n_inactive)))' (sum(n_intact)./(sum(n_stim)-sum(n_inactive)))' (sum(n_bilateral)./(sum(n_stim)-sum(n_inactive)))'    ];

figure;
bar(p_state([1 2 3 4 5 8],:),'stacked');
legend(leglabs,'Location', 'southoutside','Orientation','horizontal')
xticklabels(muscles([1 2 3 4 5 8]));
box off
ylim([0 1.1])
ylabel('% of Stims')
yticks([0 0.25 0.5 0.75 1.0])
% yticks([
title(subjectName)
            
%% Sensory Thresh vs Normal Thresh

% % sensethresh = [4000 4000 3500 3000 3000 2000 2000 2500 2500 2000 ...
% %     2000 2000 2000 2000 3000 3000 2000 2000 2000 2000 1000 1000 2000 1000 1000 2000]; %LSP05


for r = 1:length(response)
    motorthresh(r) = min([response(r).muscle.threshold]);
    allthresh(r,:) = [response(r).muscle.threshold];
end
figure; maximize; box off; nexttile;
xx = violinplot(allthresh');
for i = 1:length(xx)
    xx(i).ViolinColor = [0.9176    0.4235    0.1255];
end
hold on; hl = plot([1:26],sensethresh, '-o','LineWidth',2,...
    'MarkerEdgeColor','k','MarkerSize',10);
xlabel('Tested Electrodes')
ylabel('Threshold Amplitude (mA)');
xticklabels([response.elec]);
legend(hl, {'Sensory Threshold'});

    %'MarkerFaceColor',[.49 1 .63],...

%%
%%scatter this
nexttile;
C = linspecer(length(Muscle_ID));
for r = 1:16
scatter(elecs,allthresh(:,r),30, C(r,:),'filled');
hold on;
end
hl = scatter(elecs,sensethresh, 50,'k');
legend([Muscle_ID 'Sense Thresh'] , 'Location','eastoutside');
xticks([1:2:32]);

%% Thresh vs Thresh 9 & 16
nexttile;
C = linspecer(length(Muscle_ID)/2);

allthresh(isnan(threshes)) =  6500;
for r = 1:8
scatter(allthresh(7,r+8),allthresh(7,r),100,C(r,:),'^','filled','MarkerFaceAlpha',0.8);
hold on;
scatter(allthresh(14,r+8),allthresh(14,r),200, C(r,:),'v','filled','MarkerFaceAlpha',0.8);
end
hl = scatter(sensethresh([7 14]),sensethresh([7 14]), 120,'k');
legend([muscles 'Sensory Threshold'] , 'Location','eastoutside');
xlim([0 6500]); ylim([0 6500])
ylabel('Intact Limb'); xlabel('Residual Limb');

%% Plot Thresh V Thresh LSP02b
msub = [1 2 3 4 5 8]; % Remove Ham/TA and SO/MG for LSP02b
% % esub = [1 2 3 4]; %[5 6 7 8] %[9 10 11 12] %[33 34 35];
figure; tiledlayout(2,4); maximize; 
for m = 1:2
    switch m
            case 1
                msub1 = [4:5]; 
                msub2 = [12:13];
                titstring = [' - Flexors'];
            case 2
                msub1 = [1:3 7:8];
                msub2 = [9:11 15:16];
                titstring = [' - Extensors'];
    end
    for e = 5%1:4
      switch e
            case 1
                esub = [1 2 3 4]; %[5 6 7 8] %[9 10 11 12] %[33 34 35];
                estring = 'e1-4';
            case 2
                esub = [5 6 7 8]; %[9 10 11 12] %[33 34 35];
                estring = 'e5-8';
            case 3
                esub = [9 10 11 12]; %[33 34 35];
                estring = 'e9-12';
            case 4
                esub = [33 34 35];
                estring = 'e33-35';
            case 5
                esub = [9 16];
                estring = 'e9 vs e16';
      end
        hx(m) = nexttile; 
        for r = esub
            i = find([response.elec] == r);
            bubblechart(threshes(i,msub2), threshes(i,msub1), AUC_ratio(i,msub1), 'MarkerFaceAlpha',0.20);
            hold on;
%             ax_max(i) = max(max(response(r).muscle(m+8).p2pResponse),max(response(r).muscle(m).p2pResponse));
        end
            xlabel('Residual (mA)');
            ylabel('Intact (mA)');
            ylim([0,6500]);
            xlim([0,6500]);
            bubblesize(hx(m),[1 20])
            pl = line(xlim, ylim, 'Color', '#D3D3D3','LineStyle','--');
            box off;
            bubblelim(hx(m),[-1 1])
            axis(hx(m), 'square')
            title([estring titstring]);
            xticks([0 2000 4000 6000]);
            yticks([0 2000 4000 6000]);
            xticklabels([0 2 4 6]);
            yticklabels([0 2 4 6]);
            
% % 
% %             lgd = legend(muscles(msub1));
% %             blgd.Layout.Tile = 'west';
% %             lgd.Layout.Tile = 'west';

        

        if m == 1 && e == 4
            blgd = bubblelegend('Activation: Residual vs Intact');
            blgd.Layout.Tile = 'west';
       elseif m == 1 && e==1
            lgd = legend(muscles(msub1), 'Location', 'northwestoutside');
        elseif m == 2 && e==1
            lgd = legend(muscles(msub1), 'Location', 'southwestoutside');
        end
    end
        
end
sgtitle('Comparison of Threshold: Residual vs Intact');

%% Plot Thresh V Thresh LSP05 9 vs 16
% msub = [1 2 3 4 5 8]; % Remove Ham/TA and SO/MG for LSP02b
% % esub = [1 2 3 4]; %[5 6 7 8] %[9 10 11 12] %[33 34 35];
threshes2 = threshes;
threshes2(isnan(allthresh)) =  6500;
xtmp = find(threshes2(14,1:8) == 6000);
rtmp = [50 0 -50 0 -100]%randi([-250,50],7,1)
for i = [1 2 3 5] %to skip adjusting stimAmp of 1750 and only those at 1500/6000
    ii = i;
% %     if i == 5
% %         ii = 4;
% %     end
    threshes2(14,xtmp(i)) = threshes2(14,xtmp(i)) + rtmp(ii);
end
    
rgbcol = {'#003f5c'; '#2f4b7c'; '#665191'; '#a05195'; '#d45087'; '#f95d6a'; '#ff7c43'; '#ffa600'};
% rgbcol = {'#762a83'; '#9970ab'; '#c2a5cf'; '#e7d4e8'; '#d9f0d3';'#a6dba0'; '#5aae61'; '#1b7837'};
% C = hex2rgb(rgbcol);
D = colormap(flipud(magma(10)));
C = D(2:end,:);
% C = linspecer(16, 'qualitative');
hx = figure; %%tiledlayout(2,1, 'TileSpacing','compact'); maximize; 
         m = 0;
        for r = [9 16]
            m = m+1;
% %             hx(m) = nexttile;
            i = find([response.elec] == r);
            for c = 1:8
% %             bubblechart(threshes(i,c+8), threshes(i,c), (aucData(i,c+8)-aucData(i,c))/(aucData(i,c+8)+aucData(i,c)), C(c,:), 'MarkerFaceAlpha','flat'); %AUC_ratio(i,c)'MarkerEdgeAlpha','flat',
                if r == 9
                    scatter(threshes2(i,c+8), threshes2(i,c), 200, C(c,:), 'filled', 'd');
                elseif r == 16
                    scatter(threshes2(i,c+8), threshes2(i,c), 200, C(c,:), 'filled','h');
                end
            hold on;
            end
% %             hl = scatter(sensethresh([7 14]),[0 0], 120, 'x','k');
% % % %             title(['e' num2str(r)]);
% %            
% %             xlabel('Residual (mA)');
% %             ylabel('Intact (mA)');
% %             ylim([0,7000]);
% %             xlim([0,7000]);
% % % %             bubblesize(hx(m),[10 36])
% %             pl = line([0 6000], [6000, 6000], 'Color', '#D3D3D3','LineStyle','--');
% %             pl = line([6000 6000], [6000, 0], 'Color', '#D3D3D3','LineStyle','--');
% %             pl = line(xlim, ylim, 'Color', '#D3D3D3','LineStyle','--');
% %             box off;
% % % %             bubblelim(hx(m),[-1 1])
% %             axis(hx(m), 'square')
% %             xticks([0 2000 4000 6000]);
% %             yticks([0 2000 4000 6000]);
% %             xticklabels([0 2 4 6]);
% %             yticklabels([0 2 4 6]);
%             blgd = bubblelegend({'Total Recruitment', '(mV*mA)'});%('Activation: Residual vs Intact');
%             blgd.Layout.Tile = 'west';
%  
%             lgd = legend([muscles {'Sensory', 'Threshold'}], 'Location', 'northwestoutside');
 
         end
% %             blgd = bubblelegend({'Degree of Recruitment', 'Residual vs Intact'});%('Activation: Residual vs Intact');
% %             blgd.Layout.Tile = 'west';
% %  
 hl = scatter(sensethresh([7 14]),[0 0], 120, 'x','k');
% %             title(['e' num2str(r)]);
           
            xlabel('Residual (mA)');
            ylabel('Intact (mA)');
            ylim([0,7000]);
            xlim([0,7000]);
% %             bubblesize(hx(m),[10 36])
            pl = line([0 6000], [6000, 6000], 'Color', '#D3D3D3','LineStyle','--');
            pl = line([6000 6000], [6000, 0], 'Color', '#D3D3D3','LineStyle','--');
            pl = line(xlim, ylim, 'Color', '#D3D3D3','LineStyle','--');
            box off;
% %             bubblelim(hx(m),[-1 1])
% %             axis(hx, 'square')
            axis square;
            xticks([0 2000 4000 6000]);
            yticks([0 2000 4000 6000]);
            xticklabels([0 2 4 6]);
            yticklabels([0 2 4 6]);
            lgd = legend([muscles muscles {'Sensory', 'Threshold'}], 'Location', 'northwestoutside');
            sgtitle('Comparison of Threshold: Residual vs Intact');