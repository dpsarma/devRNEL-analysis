Muscle_ID = {'Right VM', 'Right RF', 'Right VL', 'Right BF', 'Right ST', ...
    'Right TA', 'Right SO', 'Right LG','Left VM', 'Left RF', 'Left VL',...
    'Left BF', 'Left ST', 'Left TA', 'Left MG', 'Left LG'};

% % %% Plotting Stim Averages
elecs = [9:16];
for ii = 1:length(elecs)
% %     eIdx = find(spinalElec == elecs(ii));
% %     [sortedStim,I] = sort(stimAmp(eIdx));
% %     [sortedPulse,I2] = sort(pulseWidth(eIdx));
% %     [sortedFreq,I3] = sort(stimFrequency(eIdx));
% % 

% %     results(ii).sortedStim = sortedStim;
    if length(results(ii).sortedStim) > 10
        iter = floor(length(results(ii).sortedStim)/10);
    else 
        iter = 2;
    end
% %     
% %     

% %     
% %     h(ii,1) = figure; maximize;
% %     h(ii,2) = figure; maximize;
    h(ii,3) = figure; maximize;
    for c = 9:length(Muscle_ID)
        disp(['Processing ' Muscle_ID(c)]);
% %         for fnum = 1:iter:length(eIdx)
% %             results(ii).muscle(c).meanTrace(fnum,:) = abs(trial(eIdx(I(fnum))).meanTrace{c});
% %             index = findchangepts(abs(trial(eIdx(I(fnum))).meanTrace{c}), 'MaxNumChanges',6,'Statistic','rms'); %1 is stim, 2 is start, 3 is complex start, 4 is complex end, 5 is response end, 6 is return to baseline
% %             if length(index) < 6
% %                 results(ii).muscle(c).onset(fnum) = nan;
% %                 results(ii).muscle(c).offset(fnum) = nan;
% %             else
% %                 results(ii).muscle(c).onset(fnum) = trial(eIdx(I(fnum))).epochTimeVec(index(2));
% % %                 onset_2(fnum) = trial(eIdx(I(fnum))).epochTimeVec(index(3));
% %                 results(ii).muscle(c).offset(fnum) = trial(eIdx(I(fnum))).epochTimeVec(index(5));
% %             end
% %             
% %             results(ii).muscle(c).p2p(fnum) = peak2peak(abs(trial(eIdx(I(fnum))).meanTrace{c}));
% %             results(ii).muscle(c).peak = max(abs(trial(eIdx(I(fnum))).meanTrace{c}));
% %             disp([num2str(fnum) ' of ' num2str(length(eIdx))]);
% %         end
        
% %         disp('Plotting Peak');
% %         figure(h(ii,1)); subplot(8,2,c); plot(results(ii).muscle(c).p2p); hold on; %plot(peak);
% %         title(Muscle_ID(c)); 
% %         xticks([1:iter+4:length(stimAmp)]); xticklabels({sortedStim(1:iter+4:end)/1000});
% %         xlabel('stim amp (mA)'); ylabel('emg (uV)');
% %         
% %          disp('Plotting Times');
% %         figure(h(ii,2)); subplot(8,2,c); plot(results(ii).muscle(c).onset*1000); 
% %         hold on;  plot( results(ii).muscle(c).offset*1000); legend('onset','offset');
% %         xticks([1:iter+4:length(stimAmp)]); xticklabels({sortedStim(1:iter+4:end)/1000});
% %         ylim([0,60]); 
% %         xlabel('stim amp (mA)'); ylabel('t post-stim (ms)');
% %         title(Muscle_ID(c));
        
         disp('Plotting Waterfall');
        figure(h(ii,3)); subplot(4,2,c-8);
        waterfall(results(ii).muscle(c).meanTrace)
        xlim([0,4500])
        xticks(0:1500:4500)
        xticklabels({-50;0;50;150})
%         xlabel('Time (ms)');
%         zlabel('EMG (uV)');
        yticklabels({results(ii).sortedStim(1:iter+2:end)});
        yticks([1:iter+2:length(results(ii).sortedStim)]);
        xaz =  -17.8667; xzen =  22.5799;
        view(xaz,xzen);
        title(Muscle_ID(c))
    end
% %     figure(h(ii,1)); suptitle(['Electrode ' num2str(elecs(ii)) ': P2P response by Muscle']);
% %      figure(h(ii,2)); suptitle(['Electrode ' num2str(elecs(ii)) ': Response Timing by Muscle']);
      figure(h(ii,3)); suptitle(['Electrode ' num2str(elecs(ii)) ': Response Progression by Muscle']);
end

for ii = 1:length(elecs)
% %     disp('saving peaks');
% %     saveas(h(ii,1),[reportPath '\' strrep(datestr(now),':','_')  '_PeakProgression' '_e' num2str(elecs(ii)) '.png']);
% %     savefig(h(ii,1),[reportPath '\' strrep(datestr(now),':','_')  '_PeakProgression' '_e' num2str(elecs(ii))]);
% %     disp('saving response times');
% %     saveas(h(ii,2),[reportPath '\' strrep(datestr(now),':','_')  '_ResponseTimes' '_e' num2str(elecs(ii)) '.png']);
% %     savefig(h(ii,2),[reportPath '\' strrep(datestr(now),':','_')  '_ResponseTimes' '_e' num2str(elecs(ii))]);
     disp('saving waterfalls')
    saveas(h(ii,3),[reportPath '\'  'Residual_Waterfalls' '_e' num2str(elecs(ii)) '.png']);
    savefig(h(ii,3),[reportPath '\'  'Residual_Waterfalls' '_e' num2str(elecs(ii))]);
    
end
% % 
% % % %     pause(0.1)
% % % %     close all
% %     
% % for ii = 1:length(elecs)
% %     for c = 1:length(Muscle_ID)
% %         mpeak(ii,c) = max(results(ii).muscle(c).p2p);
% %     end 
% % end
% % hp = figure;maximize;
% % for c = 1:length(Muscle_ID)
% %     hs(c) = subplot(2,8,c);
% %     plot(mpeak(:,c))
% %     xticks(1:iter+1:8)
% %     xticklabels({9:iter+1:16})
% %     xlabel('electrodes');
% %     if c == 1 || c == 9
% %         ylabel('P2P EMG (uV)');
% %     end
% %     title(Muscle_ID(c));
% % end
% % linkaxes([hs(1:7),hs(9:end)],'xy');
% % % linkaxes([hs(9:end)],'xy');
% % suptitle('Peak Activation by Electrode');
% %     saveas(hp,[reportPath '\'  '_PeaksByElectrode.png']);
% %     savefig(hp,[reportPath '\'  '_PeaksByElectrode']);