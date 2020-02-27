%% Spider Plot Scratch 2
% To make both Independent Spider Plots & to make one super plot of stim.

for lcall = 1:4

    switch lcall
        case 1
            % elecs 1-8: day9_set2
            load('D:\FigRescources\UH3\LSP05\emgRecruitment_Summary\Day 9\day9-set2.mat');
            multipolar = 'no';
            anodeElec = [];
            disp('elecs 1-8: day9_set2');
        case 2
            % elecs 9-16: Day8-Set_1_RC_1-6mA_500us_1Hz
            load('D:\FigRescources\UH3\LSP05\emgRecruitment_Summary\Day 8\Day8-Set_1_RC_1-6mA_500us_1Hz.mat');
            multipolar = 'no';
            anodeElec = [];
            disp('elecs 9-16: Day8-Set_1_RC_1-6mA_500us_1Hz');
        case 3
            % elecs 18-24: Day8-Set_2_RC_1-6mA_500us_1Hz
            load('D:\FigRescources\UH3\LSP05\emgRecruitment_Summary\Day 8\Day8-Set_2_RC_1-6mA_500us_1Hz.mat');
            multipolar = 'no';
            anodeElec = [];
            disp('elecs 18-24: Day8-Set_2_RC_1-6mA_500us_1Hz');
        case 4
            % elecs 25-32: day9_set1
            load('D:\FigRescources\UH3\LSP05\emgRecruitment_Summary\Day 9\day9-set1.mat');
            disp('elecs 25-32: day9_set1');
        otherwise
            % multipolar: day9_set3
            load('D:\FigRescources\UH3\LSP05\emgRecruitment_Summary\Day 9\day9-set3.mat');
            disp('multipolar: day9_set3');
    end

    elecs = unique(spinalElec);
    % % RBinary = zeros(length(elecs),length(Muscle_ID));

    for ii = 1:length(elecs)
        eIdx = find(spinalElec == elecs(ii));

       if ~isempty(anodeElec)
           RBinary(elecs(ii)).label = ['E' num2str(elecs(ii)) '::' num2str(anodeElec(eIdx(1)))]; 
       else
           RBinary(elecs(ii)).label = ['E' num2str(elecs(ii))];
       end
        
        % Sort Stim Amps by Electrode
        [sortedStim,I] = sort(stimAmp(eIdx));
        for mm = 1:length(Muscle_ID)
                for ll =  1:length(eIdx)
                    XS(ll,:) = max(trial(eIdx(I(ll))).baseline(mm));
                    YS(ll,:) = max(trial(eIdx(I(ll))).responseRMS(mm));
                end
                mBase = mean(XS);
                stdBase = std(XS);
                YS = YS - mBase;
                for ss = 1:length(XS)
                    if YS(ss) > 2*stdBase
                        Rstim(ss) = 1;
                    else
                        Rstim(ss) = 0;
                    end
                end
                  tmp =  find(Rstim);
                  if isempty(tmp) 
                      RBinary(elecs(ii)).thresholdAmp(mm) = 0;
                  else
                      RBinary(elecs(ii)).thresholdAmp(mm) = sortedStim(tmp(1));
                  end
                 ResponsetoStim{elecs(ii),mm} = Rstim;
                 clear XS YS Rstim
        end
    end
    disp('Finished Processing');
    clearvars -except ResponsetoStim RBinary
end



threshes = zeros(size(RBinary,2),16);
for i = 1:size(RBinary,2)
    if isempty(RBinary(i).thresholdAmp)
        continue;
    end
    
    threshes(i,:) = RBinary(i).thresholdAmp;
end
threshes(threshes== 0) = 10000;
musclesR = {'VM', 'RF', 'VL', 'BF', 'ST', 'TA', 'SO', 'LG',};
musclesL = {'VM', 'RF', 'VL', 'BF', 'ST', 'TA', 'MG', 'LG',};
elec_labels = {RBinary.label};
cmap = [(magma(20));1,1,1];

figure;maximize;
h1 = subplot(121);
imagesc(threshes(:,9:16));
xticklabels(musclesL);
yticks(1:size(RBinary,2));
yticklabels(elec_labels);
colormap(cmap);
caxis([125 6500]);
title('Left');
% ylabel('Electrodes');
colorbar('location','WestOutside', 'AxisLocation','Out',...
    'Ticks',[250,500,750,1000,1250,1500,1750,2000,2250,2500,2750,3000,3250,3500,4000,4500,5000,5500,6000]);
h2 = subplot(122);
imagesc(threshes(:,1:8));
xticklabels(musclesR);
yticks(1:size(RBinary,2));
yticklabels(elec_labels);
colormap(cmap);
caxis([125 6500]);
title('Right');
% ylabel('Electrodes');
colorbar('location','EastOutside','AxisLocation','Out',...
    'Ticks',[250,500,750,1000,1250,1500,1750,2000,2250,2500,2750,3000,3250,3500,4000,4500,5000,5500,6000]);

suptitle('Threshold Heatmap by Electrode')
