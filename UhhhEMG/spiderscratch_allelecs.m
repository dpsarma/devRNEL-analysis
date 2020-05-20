%% Plot All Electrodes Together
load('D:\DATA\UH3 testing\LSP05\data_gen\threshes\1Hz-500us_thresholds_srtd.mat');


clear leglabs rtemp C
div = 1000;
ii = 1;
for i =  1:length(RBinary) %7:14
    leglabs{ii} =  RBinary(i).label; 
    rtemp(ii,:) = 1./(RBinary(i).thresholdAmp/div);
    ii = ii+1;
end


%Remapped VL switched with RF
rtemp(rtemp==0) = nan;
rtemp(:,[2 3]) = rtemp(:,[3 2]);
musclesR = {'VM', 'VL', 'RF', 'BF', 'ST', 'TA', 'SO', 'LG',};
musclesL = {'VM', 'VL', 'RF', 'BF', 'ST', 'TA', 'MG', 'LG',};

hS1 = figure; %maximize;
spdrfill = 'on';
spidermin(1:8) =(1/(6000/div)); 
spdrmax(1:8) = (1/(250/div)); 
spdriter = 3; 

subset = [21 23];%[15:26]; %
C = linspecer(length(leglabs(subset)))
% colors = linspecer(length(leglabs));
% C = colors(subset,:);

% %     subplot(121); 
    spider_plot(rtemp(subset,9:16),...
        'AxesLabels', musclesL,'AxesInterval', 2,...
        'FillOption', spdrfill,'FillTransparency', 0.1, 'AxesLimits',...
        [spidermin; spdrmax],...
        'AxesInterval', spdriter, 'Direction', 'clockwise',...
        'Color', C);
    title('Residual Limb (medial)');
    % Legend properties
    legend(leglabs(subset), 'Location', 'southwestoutside', 'Orientation', 'vertical');
    
    
% %     subplot(122); spider_plot(rtemp(subset,1:8),...
% %         'AxesLabels', musclesR,'AxesInterval', 2,...
% %         'FillOption', spdrfill,'FillTransparency', 0.1, 'AxesLimits',...
% %         [spidermin; spdrmax],...
% %         'AxesInterval', spdriter,...
% %         'Color', C);
% %     title('Intact Limb (right)');
% %     legend(leglabs(subset), 'Location', 'eastoutside', 'Orientation', 'vertical');
% % 
% %     tmptit = sgtitle('(Inverse) Thresholds (1/mA) for Medial Lead: 1Hz, 500us');
% %     set(tmptit, 'Interpreter', 'none');
    
% % saveas(hS1,[reportPath '\' 'logInverseThresholdPlot_' 'medLead_' strrep(datestr(now),':','_') '_.png']);
% % savefig(hS1,[reportPath '\' 'logInverseThresholdPlot_' 'medLead_' strrep(datestr(now),':','_') ]);
% % saveas(hS1,[reportPath '\' strrep(datestr(now),':','_')  '_ThresholdPlot' '_' setName '.png']);
% % savefig(hS1,[reportPath '\' strrep(datestr(now),':','_')  '_ThresholdPlot' '_' setName]);

% % close all
% % 
% % hS1 = figure;maximize;
% % for mm = 1:8
% %     subplot(4,2,mm) %figure;maximize;%subplot(2,8,mm);
% %     plot(rtemp(:,mm));
% %     ylim([0,6]);
% %     xticks(1:1:length(leglabs));
% %     xticklabels({'1' '2' '3' '4' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '18' '19' '20' '23' '24' '25' '26' '27' '28' '29' '30' '31'});
% %     yticks(0:1:6)
% %     ylabel('Threshold(mA)');
% %     xlabel('electrodes');
% %     title(Muscle_ID(mm));
% % end
% % suptitle('Threshold by Electrode (Intact Limb)');
% % saveas(hS1,[reportPath '\' strrep(datestr(now),':','_')  '_ThresholdbyElecbyMusc' '_Intact' '.png']);
% % % % savefig(hS1,[reportPath '\' strrep(datestr(now),':','_')  '_ThresholdPlot' '_' setName]);

