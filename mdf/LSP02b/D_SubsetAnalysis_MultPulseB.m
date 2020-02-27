%% UH3 - LSP02 anaylsis (3b) - MDF Analysis
%----------------------------------------
% tag that is missing from all EMG trials is "Sitting" vs. "Standing"

subMusc = [4 5 12 13];

caseNums = length(pulses);
C2 = linspecer(5);
e = 4;

for m = subMusc
    figure; maximize;
    for pw =  1:caseNums
        whichPulse = pulses(pw);
        switch whichPulse
            case 200
                stims{pw} = arrayfun(@(c) E4_PW200(m,c).Stims, (1:length((E4_PW200(m,:)))));
                means{pw} = arrayfun(@(c) E4_PW200(m,c).P2P, (1:length((E4_PW200(m,:)))));
                pkVar{pw} = arrayfun(@(c) E4_PW200(m,c).P2Perr, (1:length((E4_PW200(m,:)))));
                timeOns{pw} = arrayfun(@(c) E4_PW200(m,c).timeMeanpeak, (1:length((E4_PW200(m,:))))) + dtPost
                timeOnsVar{pw} = arrayfun(@(c) E4_PW200(m,c).timeOnsVar, (1:length((E4_PW200(m,:)))));
                legEntry{pw} = ['PW 200us'];
            case 400
                stims{pw} = arrayfun(@(c) E4_PW400(m,c).Stims, (1:length((E4_PW400(m,:)))));
                means{pw} = arrayfun(@(c) E4_PW400(m,c).P2P, (1:length((E4_PW400(m,:)))));
                pkVar{pw} = arrayfun(@(c) E4_PW400(m,c).P2Perr, (1:length((E4_PW400(m,:)))));
                timeOns{pw} = arrayfun(@(c) E4_PW400(m,c).timeMeanpeak, (1:length((E4_PW400(m,:))))) + dtPost
                timeOnsVar{pw} = arrayfun(@(c) E4_PW400(m,c).timeOnsVar, (1:length((E4_PW400(m,:)))));
%                 legEntry{pw} = ['PW 400us'];
            case 600
                stims{pw} = arrayfun(@(c) E4_PW600(m,c).Stims, (1:length((E4_PW600(m,:)))));
                means{pw} = arrayfun(@(c) E4_PW600(m,c).P2P, (1:length((E4_PW600(m,:)))));
                pkVar{pw} = arrayfun(@(c) E4_PW600(m,c).P2Perr, (1:length((E4_PW600(m,:)))));
                timeOns{pw} = arrayfun(@(c) E4_PW600(m,c).timeMeanpeak, (1:length((E4_PW600(m,:))))) + dtPost
                timeOnsVar{pw} = arrayfun(@(c) E4_PW600(m,c).timeOnsVar, (1:length((E4_PW600(m,:)))));
%                 legEntry{pw} = ['PW 600us'];
            case 800
                stims{pw} = arrayfun(@(c) E4_PW800(m,c).Stims, (1:length((E4_PW800(m,:)))));
                means{pw} = arrayfun(@(c) E4_PW800(m,c).P2P, (1:length((E4_PW800(m,:)))));
                pkVar{pw} = arrayfun(@(c) E4_PW800(m,c).P2Perr, (1:length((E4_PW800(m,:)))));
                timeOns{pw} = arrayfun(@(c) E4_PW800(m,c).timeMeanpeak, (1:length((E4_PW800(m,:))))) + dtPost
                timeOnsVar{pw} = arrayfun(@(c) E4_PW800(m,c).timeOnsVar, (1:length((E4_PW800(m,:)))));
%                 legEntry{pw} = ['PW 800us'];
            case 1000
                stims{pw} = arrayfun(@(c) E4_PW1000(m,c).Stims, (1:length((E4_PW1000(m,:)))));
                means{pw} = arrayfun(@(c) E4_PW1000(m,c).P2P, (1:length((E4_PW1000(m,:)))));
                pkVar{pw} = arrayfun(@(c) E4_PW1000(m,c).P2Perr, (1:length((E4_PW1000(m,:)))));
                timeOns{pw} = arrayfun(@(c) E4_PW1000(m,c).timeMeanpeak, (1:length((E4_PW1000(m,:))))) + dtPost
                timeOnsVar{pw} = arrayfun(@(c) E4_PW1000(m,c).timeOnsVar, (1:length((E4_PW1000(m,:)))))
%                 legEntry{pw} = ['PW 1000us'];
        end
            
            
            
            [x, ia, ic] = unique([stims{pw}]);
            y = [means{pw}];
            err = [pkVar{pw}];
            tim = [timeOns{pw}];
            timerr = [timeOnsVar{pw}];
            y = y(ia);
            err = err(ia);
            tim = tim(ia);
            timerr = timerr(ia);
            
            
            
            subplot(211)
            br1 = bar(x,y,0.5,'FaceColor', C2(pw,:)); hold on;
            er1 = errorbar(x, y, err,'k', 'LineWidth', 2,'Capsize', 10);
            er1.LineStyle = 'none';
            title('Average Peak Response by Amplitude');

            subplot(212)
            br2 = bar(x,tim); hold on;
            er2 = errorbar(x, tim,timerr,'k', 'LineWidth', 2,'Capsize', 10);
            er2.LineStyle = 'none';
            title('Average Response Onset by Amplitude');


    end
    suptitle([char(musclabels(m)) ' - ' num2str(freqstring) 'Hz - Electrode' num2str(e)]);
    legend('PW 200us','','PW 400us','','PW 600us','','PW 800us','','PW 1000us','');
end