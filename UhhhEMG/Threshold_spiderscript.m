%% Plotting Thresholds

% % subjectName = 'LSP05'; 
% % reportPath = ['D:\FigRescources\UH3\' subjectName '\emgRecruitment_Summary\']; %['R:\data_generated\human\uh3_stim\' subjectName '\emgRecruitment_Summary\'];
% % 
% % 
% % %Uses the saved Threshold Values for all Electrodes to Plot Spider and
% % %Heatmap
% % load('D:\DATA\UH3 testing\LSP05\data_gen\threshes\1Hz-500us_thresholds_srtd.mat')
% % % % load('D:\FigRescources\UH3\LSP05\emgRecruitment_Summary\thresholds\1Hz-500us_thresholds.mat')
% % musclesR = {'VM', 'RF', 'VL', 'BF', 'ST', 'TA', 'SO', 'LG',};
% % musclesL = {'VM', 'RF', 'VL', 'BF', 'ST', 'TA', 'MG', 'LG',};

%Remapped
musclesR = {'VM', 'VL', 'RF', 'BF', 'ST', 'TA', 'SO', 'LG',};
musclesL = {'VM', 'VL', 'RF', 'BF', 'ST', 'TA', 'MG', 'LG',};


% % %% Plot heatmap for Thresholds
% % threshes = zeros(size(RBinary,2),16);
% % for i = 1:size(RBinary,2)
% %     if isempty(RBinary(i).thresholdAmp)
% %         continue;
% %     end
% %     
% %     threshes(i,:) = RBinary(i).thresholdAmp;
% % end
% % threshes(threshes== 0) = 10000;
% % musclesR = {'VM', 'RF', 'VL', 'BF', 'ST', 'TA', 'SO', 'LG',};
% % musclesL = {'VM', 'RF', 'VL', 'BF', 'ST', 'TA', 'MG', 'LG',};
% % elec_labels = {RBinary.label};
% % cmap = [(flipud(magma(20)));1,1,1];
% % 
% % hSt = figure; maximize;
% % h1 = subplot(121);
% % imagesc(threshes(:,9:16));
% % xticklabels(musclesL);
% % yticks(1:size(RBinary,2));
% % yticklabels(elec_labels);
% % colormap(cmap);
% % caxis([125 6500]);
% % title('Left');
% % % ylabel('Electrodes');
% % colorbar('location','WestOutside', 'AxisLocation','Out',...
% %     'Ticks',[250,500,750,1000,1250,1500,1750,2000,2250,2500,2750,3000,3250,3500,4000,4500,5000,5500,6000]);
% % h2 = subplot(122);
% % imagesc(threshes(:,1:8));
% % xticklabels(musclesR);
% % yticks(1:size(RBinary,2));
% % yticklabels(elec_labels);
% % colormap(cmap);
% % caxis([125 6500]);
% % title('Right');
% % % ylabel('Electrodes');
% % colorbar('location','EastOutside','AxisLocation','Out',...
% %     'Ticks',[250,500,750,1000,1250,1500,1750,2000,2250,2500,2750,3000,3250,3500,4000,4500,5000,5500,6000]);
% % 
% % suptitle('Threshold Heatmap by Electrode')
% % 
% % saveas(hSt,[reportPath '\' strrep(datestr(now),':','_')  'AllElecs-1Hz-500us_ThresholdHeatmap.png']);
% % savefig(hSt,[reportPath '\' strrep(datestr(now),':','_')  'AllElecs-1Hz-500us_ThresholdHeatmap']);
% % close all

% % %% Plot individual Spider plots
% % for ii = 1:size(RBinary,2)
% % RBinary(ii).thresholdAmp(RBinary(ii).thresholdAmp==0) = nan;
% % end
% % 
% % %Plot Individual Electrodes
% % spidermin(1:8) = 1000; spdrmax(1:8) = 6000; spdriter = 5; spdrfill = 'on';
% % for ii = 1%:size(RBinary,2)
% %    hS0 = figure; maximize;
% %     subplot(121); spider_plot([RBinary(ii).thresholdAmp(9:16)/1000;],...
% %         'AxesLabels', musclesL,'FillOption', spdrfill,'FillTransparency', 0.1,...
% %         'AxesLimits',[spidermin; spdrmax]/1000,'AxesInterval', spdriter,...
% %         'Direction', 'clockwise');
% %     title('Residual Limb (left)');
% %     % Legend properties
% %     legend(RBinary(ii).label, 'Location', 'southoutside', 'Orientation', 'horizontal');
% %     subplot(122); spider_plot([RBinary(ii).thresholdAmp(1:8)/1000;],...
% %         'AxesLabels', musclesR,'FillOption', spdrfill,'FillTransparency', 0.1, ...
% %         'AxesLimits', [spidermin; spdrmax]/1000,'AxesInterval', spdriter);
% %     title('Intact Limb (right)');
% %     legend(RBinary(ii).label, 'Location', 'southoutside', 'Orientation', 'horizontal');
% % 
% %     tmptit = sgtitle(['Thresholds for ' RBinary(ii).label ' - 1Hz, 500us']);
% %     set(tmptit, 'Interpreter', 'none');
% % 
% % % %     saveas(hS0,[reportPath '\' strrep(datestr(now),':','_')  '_ThresholdPlot' '_' RBinary(ii).label '.png']);
% % % %     savefig(hS0,[reportPath '\' strrep(datestr(now),':','_')  '_ThresholdPlot' '_' RBinary(ii).label]);
% % % %     close all
% % end


% % %% Plot individual Spider plots
% % for ii = 1:size(RBinary,2)
% % RBinary(ii).thresholdAmp(RBinary(ii).thresholdAmp==0) = nan;
% % end
% % 
% % %Plot Individual Electrodes
% % spidermin(1:8) = 1000; spdrmax(1:8) = 6000; spdriter = 2; spdrfill = 'on';
% % hS0 = figure; maximize;
% % leadpos_lat = (1:2:32);
% % leadpos_med = (2:2:32);
% % 
% % 
% % for ii = 1:size(RBinary,2)
% %     elec = str2num(erase(RBinary(ii).label,'E'));
% %     if  (elec <= 16)
% %         plot_pos = leadpos_lat(str2num(erase(RBinary(ii).label,'E')));
% %         disp('lateral electrode')
% %     elseif (elec > 16)
% %         plot_pos = leadpos_med(str2num(erase(RBinary(ii).label,'E'))-16);
% %         disp('medial electrode');
% %     end
% %     
% %     subplot(16,2,plot_pos);
% %     spider_plot(RBinary(ii).thresholdAmp(9:16)/1000,...
% %         'AxesLabels', musclesL,'FillOption', spdrfill,'FillTransparency', 0.1,...
% %         'AxesLimits',[spidermin; spdrmax]/1000,'AxesInterval', spdriter,...
% %         'Direction', 'clockwise');
% % %     title('Residual Limb (left)');
% %     legend(RBinary(ii).label, 'Location', 'southoutside', 'Orientation', 'horizontal');
% % 
% % % %     subplot(122); spider_plot([RBinary(ii).thresholdAmp(1:8)/1000;],...
% % % %         'AxesLabels', musclesR,'FillOption', spdrfill,'FillTransparency', 0.1, ...
% % % %         'AxesLimits', [spidermin; spdrmax]/1000,'AxesInterval', spdriter);
% % % %     title('Intact Limb (right)');
% % % %     legend(RBinary(ii).label, 'Location', 'northoutside', 'Orientation', 'horizontal');
% % 
% % %     tmptit = sgtitle(['Thresholds for ' RBinary(ii).label ' - 1Hz, 500us']);
% % %     set(tmptit, 'Interpreter', 'none');
% % 
% % 
% % end
% % 
% % % %     saveas(hS0,[reportPath '\' strrep(datestr(now),':','_')  'AllElecs_ThresholdPlot.png']);
% % % %     savefig(hS0,[reportPath '\' strrep(datestr(now),':','_')  'AllElecs_ThresholdPlot']);
% % % %     close all

%% Plot individual Spider plots with inverse
for ii = 1:size(RBinary,2)
RBinary(ii).thresholdAmp(RBinary(ii).thresholdAmp==0) = nan;
end

%Plot Individual Electrodes
spdrmax(1:8) = 1/(250/1000); 
spdrmin(1:8) = 1/(6000/1000); 
spdriter = 5; spdrfill = 'on';
for ii = 1:size(RBinary,2)
% %     if tt == 2
% %        ii = 3;
% %        disp([RBinary(tt).label 'is swapping with' RBinary(ii).label]);
% %     elseif tt == 3
% %         ii = 2;
% %         disp([RBinary(tt).label 'is swapping with' RBinary(ii).label]);
% %     else
% %         ii = tt;
% %         disp(['NoSwap']);
% %     end
    
   hS0 = figure; maximize;
    tempMesh = 1./(RBinary(ii).thresholdAmp/1000);
    tempMesh([2 3]) = tempMesh([3 2]);
    subplot(121); spider_plot([tempMesh(9:16);],...
        'AxesLabels', musclesL,'FillOption', spdrfill,'FillTransparency', 0.1,...
        'AxesLimits',[spdrmin; spdrmax],'AxesInterval', spdriter,...
        'Direction', 'clockwise');
    title('Residual Limb (left)');
    % Legend properties
    legend(RBinary(ii).label, 'Location', 'southoutside', 'Orientation', 'horizontal');
    subplot(122); spider_plot(tempMesh(1:8),...
        'AxesLabels', musclesR,'FillOption', spdrfill,'FillTransparency', 0.1, ...
        'AxesLimits', [spdrmin; spdrmax],'AxesInterval', spdriter);
    title('Intact Limb (right)');
    legend(RBinary(ii).label, 'Location', 'southoutside', 'Orientation', 'horizontal');

    tmptit = sgtitle(['(Inverse) Thresholds in mA for ' RBinary(ii).label ' - 1Hz, 500us']);
    set(tmptit, 'Interpreter', 'none');

    saveas(hS0,[reportPath '\' strrep(datestr(now),':','_')  '_ThresholdInversePlot' '_' RBinary(ii).label '.png']);
    savefig(hS0,[reportPath '\' strrep(datestr(now),':','_')  '_ThresholdInversePlot' '_' RBinary(ii).label]);
    close all
end
