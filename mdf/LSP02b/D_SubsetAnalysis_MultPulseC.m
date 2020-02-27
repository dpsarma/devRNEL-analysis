%% UH3 - LSP02 anaylsis (3c) - MDF Analysis
%----------------------------------------
% tag that is missing from all EMG trials is "Sitting" vs. "Standing"

subMusc = [12 13 4 5];

caseElecs = [4];%[2 4]
C2 = linspecer(5);
freqstring = 2;

for m = subMusc
%     figure; maximize;
    for pw =  1:length(caseElecs)
        whichElec = caseElecs(pw);
        switch whichElec
            case 2
                means = arrayfun(@(c) E2_PW200(m,c).P2P, (1:length((E2_PW200(m,:)))));
                stims = arrayfun(@(c) E2_PW200(m,c).Stims, (1:length((E2_PW200(m,:)))));
                pkVar = arrayfun(@(c) E2_PW200(m,c).P2Perr, (1:length((E2_PW200(m,:)))));
                timeOns = arrayfun(@(c) E2_PW200(m,c).timeMeanpeak, (1:length((E2_PW200(m,:))))) + dtPost;
                timeOnsVar = arrayfun(@(c) E2_PW200(m,c).timeOnsVar, (1:length((E2_PW200(m,:))))); %%%%% Calculate this from timOns!
                
                [e2pw200(m).max i] = max(means);
                if m == 4
                    e2pw200(m).max = e2pw200(m).max/10;
                end
                e2pw200(m).stims = stims(i);
                e2pw200(m).pkvar = pkVar(i);
                e2pw200(m).timons = timeOns(i);
                e2pw200(m).timeonsvar = std(timeOns);
                
                clear means stims pkVar timeOns timeOnsVar
                
                means = arrayfun(@(c) E2_PW400(m,c).P2P, (1:length((E2_PW400(m,:)))));
                stims = arrayfun(@(c) E2_PW400(m,c).Stims, (1:length((E2_PW400(m,:)))));
                pkVar = arrayfun(@(c) E2_PW400(m,c).P2Perr, (1:length((E2_PW400(m,:)))));
                timeOns = arrayfun(@(c) E2_PW400(m,c).timeMeanpeak, (1:length((E2_PW400(m,:))))) + dtPost;
                timeOnsVar = arrayfun(@(c) E2_PW400(m,c).timeOnsVar, (1:length((E2_PW400(m,:)))));
                                
                [e2pw400(m).max i] = max(means);
                e2pw400(m).stims = stims(i);
                e2pw400(m).pkvar = pkVar(i);
                e2pw400(m).timons = timeOns(i);
                e2pw400(m).timeonsvar = std(timeOns);
                 
                clear means stims pkVar timeOns timeOnsVar
                
                means = arrayfun(@(c) E2_PW600(m,c).P2P, (1:length((E2_PW600(m,:)))));
                stims = arrayfun(@(c) E2_PW600(m,c).Stims, (1:length((E2_PW600(m,:)))));
                pkVar = arrayfun(@(c) E2_PW600(m,c).P2Perr, (1:length((E2_PW600(m,:)))));
                timeOns = arrayfun(@(c) E2_PW600(m,c).timeMeanpeak, (1:length((E2_PW600(m,:))))) + dtPost;
                timeOnsVar = arrayfun(@(c) E2_PW600(m,c).timeOnsVar, (1:length((E2_PW600(m,:)))));
                
                [e2pw600(m).max i] = max(means);
                e2pw600(m).stims = stims(i);
                e2pw600(m).pkvar = pkVar(i);
                e2pw600(m).timons = timeOns(i);
                e2pw600(m).timeonsvar = std(timeOns);
                
                 clear means stims pkVar timeOns timeOnsVar
                 
                means = arrayfun(@(c) E2_PW800(m,c).P2P, (1:length((E2_PW800(m,:)))));
                stims = arrayfun(@(c) E2_PW800(m,c).Stims, (1:length((E2_PW800(m,:)))));
                pkVar = arrayfun(@(c) E2_PW800(m,c).P2Perr, (1:length((E2_PW800(m,:)))));
                timeOns = arrayfun(@(c) E2_PW800(m,c).timeMeanpeak, (1:length((E2_PW800(m,:))))) + dtPost;
                timeOnsVar = arrayfun(@(c) E2_PW800(m,c).timeOnsVar, (1:length((E2_PW800(m,:)))));
                
                [e2pw800(m).max i] = max(means);
                e2pw800(m).stims = stims(i);
                e2pw800(m).pkvar = pkVar(i);
                e2pw800(m).timons = timeOns(i);
                e2pw800(m).timeonsvar = std(timeOns);
                
                 clear means stims pkVar timeOns timeOnsVar
                 
                means = arrayfun(@(c) E2_PW1000(m,c).P2P, (1:length((E2_PW1000(m,:)))));
                stims = arrayfun(@(c) E2_PW1000(m,c).Stims, (1:length((E2_PW1000(m,:)))));
                pkVar = arrayfun(@(c) E2_PW1000(m,c).P2Perr, (1:length((E2_PW1000(m,:)))));
                timeOns = arrayfun(@(c) E2_PW1000(m,c).timeMeanpeak, (1:length((E2_PW1000(m,:))))) + dtPost;
                timeOnsVar = arrayfun(@(c) E2_PW1000(m,c).timeOnsVar, (1:length((E2_PW1000(m,:)))));
                
                [e2pw1000(m).max i] = max(means);
                e2pw1000(m).stims = stims(i);
                e2pw1000(m).pkvar = pkVar(i);
                e2pw1000(m).timons = timeOns(i);
                e2pw1000(m).timeonsvar = std(timeOns);
                
                 clear means stims pkVar timeOns timeOnsVar
                
               
            
            
            case 4
                means = arrayfun(@(c) E4_PW200(m,c).P2P, (1:length((E4_PW200(m,:)))));
                stims = arrayfun(@(c) E4_PW200(m,c).Stims, (1:length((E4_PW200(m,:)))));
                pkVar = arrayfun(@(c) E4_PW200(m,c).P2Perr, (1:length((E4_PW200(m,:)))));
                timeOns = arrayfun(@(c) E4_PW200(m,c).timeMeanpeak, (1:length((E4_PW200(m,:))))) + dtPost;
                timeOnsVar = arrayfun(@(c) E4_PW200(m,c).timeOnsVar, (1:length((E4_PW200(m,:)))));
                
                [e4pw200(m).max i] = max(means);
                e4pw200(m).stims = stims(i);
                e4pw200(m).pkvar = pkVar(i);
                e4pw200(m).timons = timeOns(i);
                e4pw200(m).timeonsvar = std(timeOns);
                
                 clear means stims pkVar timeOns timeOnsVar
                 
                means = arrayfun(@(c) E4_PW400(m,c).P2P, (1:length((E4_PW400(m,:)))));
                stims = arrayfun(@(c) E4_PW400(m,c).Stims, (1:length((E4_PW400(m,:)))));
                pkVar = arrayfun(@(c) E4_PW400(m,c).P2Perr, (1:length((E4_PW400(m,:)))));
                timeOns = arrayfun(@(c) E4_PW400(m,c).timeMeanpeak, (1:length((E4_PW400(m,:))))) + dtPost;
                timeOnsVar = arrayfun(@(c) E4_PW400(m,c).timeOnsVar, (1:length((E4_PW400(m,:)))));
                
                [e4pw400(m).max i] = max(means);
                e4pw400(m).stims = stims(i);
                e4pw400(m).pkvar = pkVar(i);
                e4pw400(m).timons = timeOns(i);
                e4pw400(m).timeonsvar = std(timeOns);
                
                 clear means stims pkVar timeOns timeOnsVar
                 
                means = arrayfun(@(c) E4_PW600(m,c).P2P, (1:length((E4_PW600(m,:)))));
                stims = arrayfun(@(c) E4_PW600(m,c).Stims, (1:length((E4_PW600(m,:)))));
                pkVar = arrayfun(@(c) E4_PW600(m,c).P2Perr, (1:length((E4_PW600(m,:)))));
                timeOns = arrayfun(@(c) E4_PW600(m,c).timeMeanpeak, (1:length((E4_PW600(m,:))))) + dtPost;
                timeOnsVar = arrayfun(@(c) E4_PW600(m,c).timeOnsVar, (1:length((E4_PW600(m,:)))));
                
                [e4pw600(m).max i] = max(means);
                e4pw600(m).stims = stims(i);
                e4pw600(m).pkvar = pkVar(i);
                e4pw600(m).timons = timeOns(i);
                e4pw600(m).timeonsvar = std(timeOns);
                
                 clear means stims pkVar timeOns timeOnsVar
                 
                 means = arrayfun(@(c) E4_PW800(m,c).P2P, (1:length((E4_PW800(m,:)))));
                stims = arrayfun(@(c) E4_PW800(m,c).Stims, (1:length((E4_PW800(m,:)))));
                pkVar = arrayfun(@(c) E4_PW800(m,c).P2Perr, (1:length((E4_PW800(m,:)))));
                timeOns = arrayfun(@(c) E4_PW800(m,c).timeMeanpeak, (1:length((E4_PW800(m,:))))) + dtPost;
                timeOnsVar = arrayfun(@(c) E4_PW800(m,c).timeOnsVar, (1:length((E4_PW800(m,:)))));
                
                [e4pw800(m).max i] = max(means);
                e4pw800(m).stims = stims(i);
                e4pw800(m).pkvar = pkVar(i);
                e4pw800(m).timons = timeOns(i);
                e4pw800(m).timeonsvar = std(timeOns);
                
                 clear means stims pkVar timeOns timeOnsVar
                 
                means = arrayfun(@(c) E4_PW1000(m,c).P2P, (1:length((E4_PW1000(m,:)))));
                stims = arrayfun(@(c) E4_PW1000(m,c).Stims, (1:length((E4_PW1000(m,:)))));
                pkVar = arrayfun(@(c) E4_PW1000(m,c).P2Perr, (1:length((E4_PW1000(m,:)))));
                timeOns = arrayfun(@(c) E4_PW1000(m,c).timeMeanpeak, (1:length((E4_PW1000(m,:))))) + dtPost;
                timeOnsVar = arrayfun(@(c) E4_PW1000(m,c).timeOnsVar, (1:length((E4_PW1000(m,:)))));
                
                [e4pw1000(m).max i] = max(means);
                e4pw1000(m).stims = stims(i);
                e4pw1000(m).pkvar = pkVar(i);
                e4pw1000(m).timons = timeOns(i);
                e4pw1000(m).timeonsvar = std(timeOns);
                
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

% % %% Plotting Electrode 2
% % e2pw200;e2pw400;e2pw600;e2pw800;e2pw1000;
% % e4pw200;e4pw400;e4pw600;e4pw800;e4pw1000;
% % 
% % 
% % % % for b = 1:length(subMusc)
% % % %     m = subMusc(b);
% % % %     power(:,b) = [e2pw200(m).max  e2pw400(m).max e2pw600(m).max e2pw800(m).max e2pw1000(m).max e4pw200(m).max e4pw400(m).max e4pw600(m).max e4pw800(m).max e4pw1000(m).max];
% % % %     onset(:,b) = [e2pw200(m).timons  e2pw400(m).timons  e2pw600(m).timons e2pw800(m).timons e2pw1000(m).timons e4pw200(m).timons e4pw400(m).timons e4pw600(m).timons e4pw800(m).timons e4pw1000(m).timons];
% % % %     powerErr(:,b) = [e2pw200(m).pkvar e2pw400(m).pkvar e2pw600(m).pkvar e2pw800(m).pkvar e2pw1000(m).pkvar e4pw200(m).pkvar e4pw400(m).pkvar e4pw600(m).pkvar e4pw800(m).pkvar e4pw1000(m).pkvar];
% % % %     onsErr(:,b) = [e2pw200(m).timeonsvar  e2pw400(m).timeonsvar  e2pw600(m).timeonsvar e2pw800(m).timeonsvar e2pw1000(m).timeonsvar e4pw200(m).timeonsvar e4pw400(m).timeonsvar e4pw600(m).timeonsvar e4pw800(m).timeonsvar e4pw1000(m).timeonsvar];
% % % % %     stimLvl(:,b) = [e2pw200(m).stims e4pw200(m).stims e2pw400(m).stims e4pw400(m).stims e2pw600(m).stims e4pw600(m).stims];
% % % % end
% % 
% % e2pw200(4).max = 0.75*e2pw200(m).max;
% % e2pw200(12).max = 0.75*e2pw200(m).max;
% % e2pw200(13).max = 0.75*e2pw200(m).max;
% % e4pw200(12).max = 0.5*e2pw200(m).max;
% % e2pw200(5).timons = 0.75*e2pw200(m).timons;
% % 
% % for b = 1:length(subMusc)
% %     m = subMusc(b);      
% %     power(:,b) = [e2pw200(m).max  e2pw400(m).max e2pw600(m).max e2pw800(m).max e2pw1000(m).max];
% %     onset(:,b) = [e2pw200(m).timons  e2pw400(m).timons  e2pw600(m).timons e2pw800(m).timons e2pw1000(m).timons];
% %     powerErr(:,b) = [e2pw200(m).pkvar e2pw400(m).pkvar e2pw600(m).pkvar e2pw800(m).pkvar e2pw1000(m).pkvar];
% %     onsErr(:,b) = [e2pw200(m).timeonsvar  e2pw400(m).timeonsvar  e2pw600(m).timeonsvar e2pw800(m).timeonsvar e2pw1000(m).timeonsvar];
% % %     stimLvl(:,b) = [e2pw200(m).stims e4pw200(m).stims e2pw400(m).stims e4pw400(m).stims e2pw600(m).stims e4pw600(m).stims];
% % end
% % 
% % hF(1) = figure(1); maximize;
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
% % legend('Electrode 2, 200us', 'Electrode 2, 400us', 'Electrode 2, 600us', 'Electrode 2, 800us', 'Electrode 2, 1000us','Location', 'best');
% % ylabel('Peak to Peak Response (mV)');
% % set(gca,'FontSize',20);
% % 
% % hF(2) = figure(2); maximize;
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
% % legend('Electrode 2, 200us', 'Electrode 2, 400us', 'Electrode 2, 600us', 'Electrode 2, 800us', 'Electrode 2, 1000us', 'Location', 'best');
% % ylabel('Response Onset (sec)');
% % set(gca,'FontSize',20);
% % 
% % clear power powerErr onsErr stimLvl onset

%% Electrode 4
for b = 1:length(subMusc)
    m = subMusc(b);
    power(:,b) = [e4pw200(m).max e4pw400(m).max e4pw600(m).max e4pw800(m).max e4pw1000(m).max];
    onset(:,b) = [e4pw200(m).timons e4pw400(m).timons e4pw600(m).timons e4pw800(m).timons e4pw1000(m).timons];
    powerErr(:,b) = [e4pw200(m).pkvar e4pw400(m).pkvar e4pw600(m).pkvar e4pw800(m).pkvar e4pw1000(m).pkvar];
    onsErr(:,b) = [e4pw200(m).timeonsvar e4pw400(m).timeonsvar e4pw600(m).timeonsvar e4pw800(m).timeonsvar e4pw1000(m).timeonsvar];
%     stimLvl(:,b) = [e2pw200(m).stims e4pw200(m).stims e2pw400(m).stims e4pw400(m).stims e2pw600(m).stims e4pw600(m).stims];
end

hF(3) = figure(3); maximize;
% clist = categorical(musclabels(subMusc))
clist = 1:4;
hBar3 = bar(clist,power'); 
for k1 = 1:size(power,1)
    clist(k1,:) = bsxfun(@plus, hBar3(1).XData, [hBar3(k1).XOffset]');
    ydt(k1,:) = hBar3(k1).YData;
end
hold on;

er3 = errorbar(clist,ydt,powerErr,'k', 'LineWidth', 2,'Capsize', 10);
for k2 = 1:size(er3,2)
    er3(k2).LineStyle = 'none';
end

set(gca, 'XTickLabel', musclabels(subMusc))
legend('Electrode 4, 200us', 'Electrode 4, 400us', 'Electrode 4, 600us', 'Electrode 4, 800us', 'Electrode 4, 1000us','Location', 'best');
ylabel('Peak to Peak Response (mV)');
set(gca,'FontSize',20);

hF(4) = figure(4); maximize;
% clist = categorical(musclabels(subMusc))
clist2 = 1:4;
hBar4 = bar(clist2,onset'); 
for k1 = 1:size(onset,1)
    clist2(k1,:) = bsxfun(@plus, hBar4(1).XData, [hBar4(k1).XOffset]');
    ydt2(k1,:) = hBar4(k1).YData;
end
hold on;

er4 = errorbar(clist2,ydt2,onsErr,'k', 'LineWidth', 2,'Capsize', 10);
for k2 = 1:size(er4,2)
    er4(k2).LineStyle = 'none';
end

set(gca, 'XTickLabel', musclabels(subMusc))
legend('Electrode 4, 200us', 'Electrode 4, 400us', 'Electrode 4, 600us', 'Electrode 4, 800us', 'Electrode 4, 1000us','Location', 'best');
ylabel('Response Onset (sec)');
set(gca,'FontSize',20);

clear power powerErr onsErr stimLvl onset

figuresave = 'R:\users\dsarma\LSP02b\Figures\';
for cnt = 1:length(hF)
    savefig(hF(cnt),[figuresave 'PulseResponseSummary_' num2str(cnt)]);
    saveas(hF(cnt),[figuresave 'PulseResponseSummary_' num2str(cnt) '.png']);
end