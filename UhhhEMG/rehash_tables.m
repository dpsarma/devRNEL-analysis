%% Threshold Tables
load('D:\FigRescources\UH3\LSP05\rehash\responses_500us_1Hz_allElecsRC.mat')

for i = 1:length(response2)
    threshTable(i,:) = [response2(i).muscle(:).threshold];
end
% threshTable(:,[2 3]) = threshTable(:,[3 2]);
musclesL = {'VM', 'VL', 'RF', 'BF', 'ST', 'TA', 'MG', 'LG',};
elecs = [response2.elec];


T = table(elecs', threshTable(:,1), threshTable(:,3), threshTable(:,2), threshTable(:,4), threshTable(:,5), threshTable(:,6), threshTable(:,7), threshTable(:,8),'VariableNames',['Elecs', musclesL]);



%% Activity Count
% Count the range of during each muscle is active
for i = 1:length(response2)
    for m = 9:16
        rCount(i,(m-8)) = length(response2(i).StimAmps(response2(i).muscle(m).threshIdx:end))
        rRatio(i,(m-8)) =  rCount(i,(m-8))/length(response2(i).StimAmps);
    end
end


%% Plot again all together
tiledlayout(4,2)
C = linspecer(8)
nexttile;
for m = 1:3
    plot(threshTable(:,m),'LineWidth',4,'Color',C(m,:)); 
    ylim([0 6000]);
    xticks(1:2:length(elecs))
    xticklabels(elecs(1:2:end));
    hold on;
end
title('Knee Flexors: Threshold');
legend(musclesL(1:3))
ylabel('Threshold Amp');
    
nexttile;
for m = 1:3
    plot(rRatio(:,m),'LineStyle',':','LineWidth',4,'Color',C(m,:)); 
    ylim([0 1]); yticks(0:.25:1); hold on;
        xticks(1:2:length(elecs))
    xticklabels(elecs(1:2:end));
end
title('Knee Flexors: %Responding');
legend(musclesL(1:3),'Location','southeast')
ylabel('Ratio of Tested Amps');
    
nexttile;
for m = 4:5 
    plot(threshTable(:,m),'LineWidth',4,'Color',C(m,:)); 
    title(musclesL(m));
    ylim([0 6000]);
        xticks(1:2:length(elecs))
    xticklabels(elecs(1:2:end))
    hold on;    
end
title('Knee Extensors: Threshold');
legend(musclesL(4:5))
ylabel('Threshold Amp');
    
nexttile;
for m = 4:5

    plot(rRatio(:,m),'LineStyle',':','LineWidth',4,'Color',C(m,:)); 
    ylim([0 1]); yticks(0:.25:1);
    hold on;
        xticks(1:2:length(elecs))
    xticklabels(elecs(1:2:end));
end
title('Knee Extensors: %Responding');
legend(musclesL(4:5),'Location','southeast')
 ylabel('Ratio of Tested Amps');   

nexttile;
for m = 6
    plot(threshTable(:,m),'LineWidth',4,'Color',C(m,:)); 
    title(musclesL(m));
    ylim([0 6000]);
    xticks(1:2:length(elecs))
    xticklabels(elecs(1:2:end))
    hold on;
end
title('Ankle Flexor: Threshold');
legend(musclesL(6))
ylabel('Threshold Amp');
    
nexttile;
for m = 6

    plot(rRatio(:,m),'LineStyle',':','LineWidth',4,'Color',C(m,:)); 
    ylim([0 1]); yticks(0:.25:1); hold on;
        xticks(1:2:length(elecs))
    xticklabels(elecs(1:2:end));
end
title('Ankle Flexor: %Responding');
legend(musclesL(6),'Location','southeast')
ylabel('Ratio of Tested Amps');

nexttile;
for m = 7:8
    plot(threshTable(:,m),'LineWidth',4,'Color',C(m,:)); 
    title(musclesL(m));
    ylim([0 6000]);
    xticks(1:2:length(elecs));
    xticklabels(elecs(1:2:end));
    hold on;
end
title('Ankle Extensors: Threshold');
legend(musclesL(7:8))
xlabel('electrode');
ylabel('Threshold Amp');
    
nexttile;
for m = 7:8

    plot(rRatio(:,m),'LineStyle',':','LineWidth',4,'Color',C(m,:)); 
    ylim([0 1]); yticks(0:.25:1);
    hold on;
        xticks(1:2:length(elecs))
    xticklabels(elecs(1:2:end));
end
title('Ankle Extensors: %Responding');
legend(musclesL(7:8),'Location','southeast')
 xlabel('electrode');  
ylabel('Ratio of Tested Amps');

% nexttile; 
% ylabel('Threshold Amp');ylim([0 6000])
% xlabel('electrode'); xlim([0 length(elecs)]);
% yyaxis right
% ylabel('Ratio of Active Range'); ylim([0 1]);

sgtitle('Threshold Trends by Muscle for All Elecs');

figure; maximize;
tiledlayout(length(elecs)/2,4)
C = linspecer(8);

for i = 14:26
    nexttile;
    for m = [9 11 10]
        plot(response2(i).StimAmps/1000, [response2(i).muscle(m).p2pResponse],...
            'LineWidth',4,'Color',C(m-8,:)); 
        hold on;
    end
    title([num2str(elecs(i)) ' - Knee Flexors']);
    legend(musclesL(1:3),'Location','northwest')
    ylabel('p2p EMG (uV)');
        xlim([0 6])
        
    nexttile;
    for m = [12 13]
        plot(response2(i).StimAmps/1000, [response2(i).muscle(m).p2pResponse],...
            'LineWidth',4,'Color',C(m-8,:)); 
        hold on;
    end
    title([num2str(elecs(i)) ' - Knee Extensors']);
    legend(musclesL(4:5),'Location','northwest')
    ylabel('p2p EMG (uV)');
        xlim([0 6])
        
    nexttile;
    for m = [14]
        plot(response2(i).StimAmps/1000, [response2(i).muscle(m).p2pResponse],...
            'LineWidth',4,'Color',C(m-8,:)); 
        hold on;
    end
    title([num2str(elecs(i)) ' - Ankle Flexor']);
    legend(musclesL(6),'Location','northwest')
    ylabel('p2p EMG (uV)');
        xlim([0 6])
        
    nexttile;
    for m = [15:16]
        plot(response2(i).StimAmps/1000, [response2(i).muscle(m).p2pResponse],...
            'LineWidth',4,'Color',C(m-8,:)); 
        hold on;
    end
    title([num2str(elecs(i)) ' - Ankle Extensors']);
    legend(musclesL(7:8),'Location','northwest')
    ylabel('p2p EMG (uV)');
    xlim([0 6])
end