%% UH3 - LSP02 anaylsis (4b) - MDF Analysis
%----------------------------------------
% tag that is missing from all EMG trials is "Sitting" vs. "Standing"

subMusc = [12 13 4 5];

caseElecs = [2];%[2 4]
C2 = linspecer(5);
freqstring = 1;

E2_F5(4,13).pkVar = NaN;
E2_F5(4,13).timeOnsVar = NaN;
E2_F5(4,13).maxMeanpeak = NaN;
E2_F5(4,13).timeMeanpeak = NaN;


for m = subMusc
    
    for pw =  1:length(caseElecs)
        whichElec = caseElecs(pw);
        switch whichElec
            case 2
                means = arrayfun(@(c) E2_F1(m,c).P2P, (1:length((E2_F1(m,:)))));
                stims = arrayfun(@(c) E2_F1(m,c).Stims, (1:length((E2_F1(m,:)))));
                pkVar = arrayfun(@(c) E2_F1(m,c).P2Perr, (1:length((E2_F1(m,:)))));
                timeOns = arrayfun(@(c) E2_F1(m,c).timeMeanpeak, (1:length((E2_F1(m,:))))) + dtPost;
                timeOnsVar = arrayfun(@(c) E2_F1(m,c).timeOnsVar, (1:length((E2_F1(m,:)))));
                
                [e2f1(m).max i] = max(means);
                e2f1(m).stims = stims(i);
                e2f1(m).pkvar = pkVar(i);
                e2f1(m).timons = timeOns(i);
                e2f1(m).timeonsvar = std(timeOns);
              
                if m == 4
                    e2f1(m).max = e2f1(m).max/10;
                    e2f1(m).timeonsvar = 0.85*e2f1(m).timeonsvar;
                end
                
                clear means stims pkVar timeOns timeOnsVar
                
                means = arrayfun(@(c) E2_F2(m,c).P2P, (1:length((E2_F2(m,:)))));
                stims = arrayfun(@(c) E2_F2(m,c).Stims, (1:length((E2_F2(m,:)))));
                pkVar = arrayfun(@(c) E2_F2(m,c).P2Perr, (1:length((E2_F2(m,:)))));
                timeOns = arrayfun(@(c) E2_F2(m,c).timeMeanpeak, (1:length((E2_F2(m,:))))) + dtPost;
                timeOnsVar = arrayfun(@(c) E2_F2(m,c).timeOnsVar, (1:length((E2_F2(m,:)))));
                                
                [e2f2(m).max i] = max(means);
                if m == 4
                    e2f2(m).max = e2f2(m).max/10;
                end
                e2f2(m).stims = stims(i);
                e2f2(m).pkvar = pkVar(i);
                e2f2(m).timons = timeOns(i);
                e2f2(m).timeonsvar = std(timeOns);
                 
                clear means stims pkVar timeOns timeOnsVar
                
                means = arrayfun(@(c) E2_F5(m,c).P2P, (1:length((E2_F5(m,:)))));
                stims = arrayfun(@(c) E2_F5(m,c).Stims, (1:length((E2_F5(m,:)))));
                pkVar = arrayfun(@(c) E2_F5(m,c).P2Perr, (1:length((E2_F5(m,:)))));
                timeOns = arrayfun(@(c) E2_F5(m,c).timeMeanpeak, (1:length((E2_F5(m,:))))) + dtPost;
                timeOnsVar = arrayfun(@(c) E2_F5(m,c).timeOnsVar, (1:length((E2_F5(m,:)))));
                
                [e2f5(m).max i] = max(means);
                e2f5(m).stims = stims(i);
                e2f5(m).pkvar = pkVar(i);
                e2f5(m).timons = timeOns(i);
                e2f5(m).timeonsvar = std(timeOns);
                if m == 4
                    e2f5(m).max = e2f5(m).max/10;
                    e2f5(m).timeonsvar = timeOnsVar(i);
                end
                                  
                
                 clear means stims pkVar timeOns timeOnsVar
                
               
            
            
            case 4
                means = arrayfun(@(c) E4_F1(m,c).P2P, (1:length((E4_F1(m,:)))));
                stims = arrayfun(@(c) E4_F1(m,c).Stims, (1:length((E4_F1(m,:)))));
                pkVar = arrayfun(@(c) E4_F1(m,c).P2Perr, (1:length((E4_F1(m,:)))));
                timeOns = arrayfun(@(c) E4_F1(m,c).timeMeanpeak, (1:length((E4_F1(m,:))))) + dtPost;
                timeOnsVar = arrayfun(@(c) E4_F1(m,c).timeOnsVar, (1:length((E4_F1(m,:)))));
                
                [e4f1(m).max i] = max(means);
                e4f1(m).stims = stims(i);
                e4f1(m).pkvar = pkVar(i);
                e4f1(m).timons = timeOns(i);
                e4f1(m).timeonsvar = std(timeOns);
                
                 clear means stims pkVar timeOns timeOnsVar
                 
                means = arrayfun(@(c) E4_F2(m,c).P2P, (1:length((E4_F2(m,:)))));
                stims = arrayfun(@(c) E4_F2(m,c).Stims, (1:length((E4_F2(m,:)))));
                pkVar = arrayfun(@(c) E4_F2(m,c).P2Perr, (1:length((E4_F2(m,:)))));
                timeOns = arrayfun(@(c) E4_F2(m,c).timeMeanpeak, (1:length((E4_F2(m,:))))) + dtPost;
                timeOnsVar = arrayfun(@(c) E4_F2(m,c).timeOnsVar, (1:length((E4_F2(m,:)))));
                
                [e4f2(m).max i] = max(means);
                e4f2(m).stims = stims(i);
                e4f2(m).pkvar = pkVar(i);
                e4f2(m).timons = timeOns(i);
                e4f2(m).timeonsvar = std(timeOns);
                
                 clear means stims pkVar timeOns timeOnsVar
                 
                means = arrayfun(@(c) E4_F5(m,c).P2P, (1:length((E4_F5(m,:)))));
                stims = arrayfun(@(c) E4_F5(m,c).Stims, (1:length((E4_F5(m,:)))));
                pkVar = arrayfun(@(c) E4_F5(m,c).P2Perr, (1:length((E4_F5(m,:)))));
                timeOns = arrayfun(@(c) E4_F5(m,c).timeMeanpeak, (1:length((E4_F5(m,:))))) + dtPost;
                timeOnsVar = arrayfun(@(c) E4_F5(m,c).timeOnsVar, (1:length((E4_F5(m,:)))));
                
                [e4f5(m).max i] = max(means);
                e4f5(m).stims = stims(i);
                e4f5(m).pkvar = pkVar(i);
                e4f5(m).timons = timeOns(i);
                e4f5(m).timeonsvar = std(timeOns);
                
                 clear means stims pkVar timeOns timeOnsVar
        end
            
%             [x, ia, ic] = unique([stims{pw}]);
%             y = [means{pw}];
%             err = [pkVar{pw}];
%             tim = [timeOns{pw}];
%             timerr = [timeOnsVar{pw}];
%             y = y(ia);
%             err = err(ia);
%             tim = tim(ia);
%             timerr = timerr(ia);
%             
%             
%             
%             subplot(211)
%             br1 = bar(x,y,0.5,'FaceColor', C2(pw,:)); hold on;
%             er1 = errorbar(x, y, err,'k', 'LineWidth', 2,'Capsize', 10);
%             er1.LineStyle = 'none';
%             title('Average Peak Response by Amplitude');
% 
%             subplot(212)
%             br2 = bar(x,tim); hold on;
%             er2 = errorbar(x, tim,timerr,'k', 'LineWidth', 2,'Capsize', 10);
%             er2.LineStyle = 'none';
%             title('Average Response Onset by Amplitude');


    end
%     suptitle([char(musclabels(m)) ' - ' num2str(freqstring) 'Hz - Electrode' num2str(whichElec)]);
%     legend('Electrode 2','','Electrode 4','');
end


% % e2f1;e2f2;e2f5;
% % e4f1;e4f2;e4f5;




% % for b = 1:length(subMusc)
% %     m = subMusc(b);
% %     power(:,b) = [e2f1(m).max  e2f2(m).max e2f5(m).max e4f1(m).max e4f2(m).max e4f5(m).max];
% %     onset(:,b) = [e2f1(m).timons  e2f2(m).timons  e2f5(m).timons e4f1(m).timons e4f2(m).timons e4f5(m).timons];
% %     powerErr(:,b) = [e2f1(m).pkvar e2f2(m).pkvar e2f5(m).pkvar e4f1(m).pkvar e4f2(m).pkvar e4f5(m).pkvar];
% %     onsErr(:,b) = [e2f1(m).timeonsvar  e2f2(m).timeonsvar  e2f5(m).timeonsvar e4f1(m).timeonsvar e4f2(m).timeonsvar e4f5(m).timeonsvar];
% % %     stimLvl(:,b) = [e2f1(m).stims e4f1(m).stims e2f2(m).stims e4f2(m).stims e2f5(m).stims e4f5(m).stims];
% % end

for b = 1:length(subMusc)
    m = subMusc(b);
    power(:,b) = [e2f1(m).max  e2f2(m).max e2f5(m).max];
    onset(:,b) = [e2f1(m).timons  e2f2(m).timons  e2f5(m).timons];
    powerErr(:,b) = [e2f1(m).pkvar e2f2(m).pkvar e2f5(m).pkvar];
    onsErr(:,b) = [e2f1(m).timeonsvar  e2f2(m).timeonsvar  e2f5(m).timeonsvar];
end


hF(1) = figure(1); maximize;
% clist = categorical(musclabels(subMusc))
clist = 1:4;
hBar = bar(clist,power'); 
for k1 = 1:size(power,1)
    clist(k1,:) = bsxfun(@plus, hBar(1).XData, [hBar(k1).XOffset]');
    ydt(k1,:) = hBar(k1).YData;
end
hold on;

er1 = errorbar(clist,ydt,powerErr,'k', 'LineWidth', 2,'Capsize', 10);
for k2 = 1:size(er1,2)
    er1(k2).LineStyle = 'none';
end

set(gca, 'XTickLabel', musclabels(subMusc))
legend('Electrode 2, 1Hz', 'Electrode 2, 2Hz', 'Electrode 2, 5Hz','Location','best');
ylabel('Peak to Peak Response (uV)');
set(gca,'FontSize',20);

hF(2) = figure(2); maximize;
% clist = categorical(musclabels(subMusc))
clist2 = 1:4;
hBar2 = bar(clist2,onset'); 
for k1 = 1:size(onset,1)
    clist2(k1,:) = bsxfun(@plus, hBar2(1).XData, [hBar2(k1).XOffset]');
    ydt2(k1,:) = hBar2(k1).YData;
end
hold on;

er2 = errorbar(clist2,ydt2,onsErr,'k', 'LineWidth', 2,'Capsize', 10);
for k2 = 1:size(er2,2)
    er2(k2).LineStyle = 'none';
end

set(gca, 'XTickLabel', musclabels(subMusc))
legend('Electrode 2, 1Hz', 'Electrode 2, 2Hz', 'Electrode 2, 5Hz','Location','best');
ylabel('Response Onset (sec)');
set(gca,'FontSize',20);

clear power powerErr onsErr stimLvl onset

% % %% Electrode 4 separately
% % for b = 1:length(subMusc)
% %     m = subMusc(b);
% %     power(:,b) = [e4f1(m).max e4f2(m).max e4f5(m).max];
% %     onset(:,b) = [e4f1(m).timons e4f2(m).timons e4f5(m).timons];
% %     powerErr(:,b) = [e4f1(m).pkvar e4f2(m).pkvar e4f5(m).pkvar];
% %     onsErr(:,b) = [e4f1(m).timeonsvar e4f2(m).timeonsvar e4f5(m).timeonsvar];
% % %     stimLvl(:,b) = [e2f1(m).stims e4f1(m).stims e2f2(m).stims e4f2(m).stims e2f5(m).stims e4f5(m).stims];
% % end
% % 
% % hF(3) = figure(3); maximize;
% % % clist = categorical(musclabels(subMusc))
% % clist = 1:4;
% % hBar = bar(clist,power'); 
% % for k1 = 1:size(power,1)
% %     clist(k1,:) = bsxfun(@plus, hBar(1).XData, [hBar(k1).XOffset]');
% %     ydt(k1,:) = hBar(k1).YData;
% % end
% % hold on;
% % 
% % er1 = errorbar(clist,ydt,powerErr,'k', 'LineWidth', 2,'Capsize', 10);
% % for k2 = 1:size(er1,2)
% %     er1(k2).LineStyle = 'none';
% % end
% % 
% % set(gca, 'XTickLabel', musclabels(subMusc))
% % legend('Electrode 4, 1Hz', 'Electrode 4, 2Hz', 'Electrode 4, 5Hz','Location','best');
% % ylabel('Peak to Peak Response (uV)');
% % set(gca,'FontSize',20);
% % 
% % hF(4) = figure(4); maximize;
% % % clist = categorical(musclabels(subMusc))
% % clist2 = 1:4;
% % hBar2 = bar(clist2,onset'); 
% % for k1 = 1:size(onset,1)
% %     clist2(k1,:) = bsxfun(@plus, hBar2(1).XData, [hBar2(k1).XOffset]');
% %     ydt2(k1,:) = hBar2(k1).YData;
% % end
% % hold on;
% % 
% % er2 = errorbar(clist2,ydt2,onsErr,'k', 'LineWidth', 2,'Capsize', 10);
% % for k2 = 1:size(er2,2)
% %     er2(k2).LineStyle = 'none';
% % end
% % 
% % set(gca, 'XTickLabel', musclabels(subMusc))
% % legend('Electrode 4, 1Hz', 'Electrode 4, 2Hz', 'Electrode 4, 5Hz','Location','best');
% % ylabel('Response Onset (sec)');
% % set(gca,'FontSize',20);
% % 
% % clear power powerErr onsErr stimLvl onset
% % 
% % figuresave = 'R:\users\dsarma\LSP02b\Figures\';
% % for cnt = 1:length(hF)
% %     savefig(hF(cnt),[figuresave 'FreqResponseSummary_' num2str(cnt)]);
% %     saveas(hF(cnt),[figuresave 'FreqResponseSummary_' num2str(cnt) '.png']);
% % end
% % 
