%% Plot Pulse Width and Normal RC Together

    %Sort Stim for Normal RCs
    [sortedStim,I] = sort(stimAmpRC);
    %Sort PulseWidths for PW RC
    [sortedPulse,I2] = sort(pulseWidthPW);
    
    sortedCharge_RC =  sortedStim.*pulseWidthRC(I);  
    sortedCharge_PW =  sortedPulse.*stimAmpPW(I2); 
    

    %% Plotting Recruitment Curves
    responseRC = [trialRC(:).responseRMS]';
    responsePW = [trialPW(:).responseRMS]';
    

    hF1 = figure; maximize;
    
    
    
%     musclesR = {'VM', 'RF', 'VL', 'BF', 'ST', 'TA', 'SO', 'LG',};
%     musclesL = {'VM', 'RF', 'VL', 'BF', 'ST', 'TA', 'MG', 'LG',};
    
        for m = 1:length(Muscle_ID)
            
            for ll =  1:length(trialRC)
                XRC(ll,:) = trialRC(ll).baseline(m);              
            end
            mBaseRC = mean(XRC);
            for pp =  1:length(trialPW)
                XPW(pp,:) = trialPW(pp).baseline(m);              
            end
            mBasePW = mean(XPW);
            
            hFs(m) = subplot(2,8,m)
            yyA = smooth(sortedCharge_RC,responseRC(I,m)-mBaseRC,'sgolay'); %/max(response(:,m)) to normalize
            plot(sortedCharge_RC,yyA,'LineWidth',3);
            hold on;
            yyB = smooth(sortedCharge_PW,responsePW(I2,m)-mBasePW,'sgolay'); %/max(response(:,m)) to normalize
            plot(sortedCharge_PW,yyB,'LineWidth',3);
            hold on;
            ylabel('(post-Stim) Mean RMS  (uV)'); %'(post-Stim) Mean RMS  (uV)'
            xlabel('Charge Delivered (uA*us)');
            title([Muscle_ID(m)]);
            legend('Ampl Mod','PW Mod');
        end
        linkaxes([hFs],'xy');

sgtitle('Comparison of RC as Charge Delivered PW vs Amp, Day 12 E 16');


saveas(hF1,[reportPath '\' strrep(datestr(now),':','_')  '_Charge Balanced RCs-Day12.png']);
savefig(hF1,[reportPath '\' strrep(datestr(now),':','_') '_Charge Balanced RCs-Day12']);
    
    


