%%%% UH3 EMG - Recruitment Curves - Quick Load Trial Stats - Testing Version %%
%% Basic Plotting for Frequency and Set on Elec 9 and 16
    freqs = unique([trial.stimFrequency]);
    pwidths = unique([trial.pulseWidth]);
    elecs = unique([trial.spinalElec]);

    %Proto code for looping
% % for p = 1:length(pwidths)
% %     check_p = pwidths(p);
% %     disp(['Pulse is ' num2str(check_p)]);
% %     for f = 1:length(freqs)
% %         check_f = freqs(f);
% %         disp(['Freq is ' num2str(check_p)]);
% %         eIdx = find(([trial.stimFrequency] == check_f) & ([trial.pulseWidth] == check_p));
% %         for t = 1:length(eIdx)
% %             n = eIdx(t);
% %             for i = 1:length(msubset)
% %                 m = msubset(i);
% %             end
% %         end    
% %     end
% % end

% Swapped RF & VL so 3<->2 & 11<->10
msubset = [9:16]; %or 1:8 for right/intact
msubset([2 3]) = msubset([3 2]);

% % %% Plot Raw Traces
% % 
% % for p = 1:length(pwidths)
% %     check_p = pwidths(p);
% %     disp(['Pulse is ' num2str(check_p)]);
% %     for f = 1:length(freqs)
% %         check_f = freqs(f);
% %         disp(['Freq is ' num2str(check_p)]);
% %         eIdx = find(([trial.stimFrequency] == check_f) & ([trial.pulseWidth] == check_p));
% %         for t = 1:length(eIdx)
% %             n = eIdx(t);
% %             disp(['Trial is ' num2str(n)]);
% %             h(t) = figure;
% %             set(h(t), 'Visible', 'off');
% %             maximize;
% %             disp(['Plotting Full Trace']);
% %             for i = 1:length(msubset)
% %                 m = msubset(i);
% %                 subplot(2,length(msubset)/2,i);
% %                 plot(trial(n).timeVec, trial(n).emgF(m,:));                 
% %                 hold on;
% %                 vline(trial(n).stims/fs);
% %                 title(Muscle_ID(m));ylabel('EMG (uV)'); xlabel('Time (s)');                 
% %             end
% %             tit = suptitle(['E' num2str(trial(n).spinalElec) '- Full Trace: ' num2str(trial(n).stimFrequency) 'Hz_' ... 
% %                 num2str(trial(n).pulseWidth) 'us_' num2str(trial(n).stimAmp/1000) 'mA']);
% %             set(tit, 'Interpreter', 'none') 
% %             disp(['Saving Full Trace ' num2str(n)]);
% %             saveas(h(t),[setPath '\FullTraces\' setDescrpt '- Full Trace_' num2str(trial(n).stimFrequency) 'Hz_' ... 
% %                 num2str(trial(n).pulseWidth) 'ms_' num2str(trial(n).stimAmp/1000) 'mA.png']); %strrep(datestr(now),':','_')
% %             close;
% %         end    
% %     end
% % end
% % 
% % %% Plot Rectified Traces
% % 
% % for p = 3:length(pwidths)
% %     check_p = pwidths(p);
% %     disp(['Pulse is ' num2str(check_p)]);
% %     for f = 1:length(freqs)
% %         check_f = freqs(f);
% %         disp(['Freq is ' num2str(check_p)]);
% %         eIdx = find(([trial.stimFrequency] == check_f) & ([trial.pulseWidth] == check_p));
% %         for t = 1:length(eIdx)
% %             n = eIdx(t);
% %             disp(['Trial is ' num2str(n)]);
% %             h(t) = figure;
% %             set(h(t), 'Visible', 'off');
% %             maximize;
% %             disp(['Plotting Full Trace']);
% %             for i = 1:length(msubset)
% %                 m = msubset(i);
% %                 subplot(2,length(msubset)/2,i);
% %                 plot(trial(n).timeVec, abs(trial(n).emgF(m,:)));                 
% %                 hold on;
% %                 vline(trial(n).stims/fs);
% %                 title(Muscle_ID(m));ylabel('EMG (uV)'); xlabel('Time (s)');                 
% %             end
% %             tit = suptitle(['E' num2str(trial(n).spinalElec) '- Rectified Full Trace: ' num2str(trial(n).stimFrequency) 'Hz_' ... 
% %                 num2str(trial(n).pulseWidth) 'us_' num2str(trial(n).stimAmp/1000) 'mA']);
% %             set(tit, 'Interpreter', 'none') 
% %             disp(['Saving Full Trace ' num2str(n)]);
% %             saveas(h(t),[setPath '\RectFullTraces\' 'E' num2str(trial(n).spinalElec) '- Rectified Full Trace_' num2str(trial(n).stimFrequency) 'Hz_' ... 
% %                 num2str(trial(n).pulseWidth) 'ms_' num2str(trial(n).stimAmp/1000) 'mA.png']); %strrep(datestr(now),':','_')
% %             close;
% %         end    
% %     end
% % end
% % 
% % %% Plotting Waterfalls
% % for p = 1:length(pwidths)
% %     check_p = pwidths(p);
% %     disp(['Pulse is ' num2str(check_p)]);
% %     for f = 1:length(freqs)
% %         check_f = freqs(f);
% %         disp(['Freq is ' num2str(check_p)]);
% %         eIdx = find(([trial.stimFrequency] == check_f) & ([trial.pulseWidth] == check_p));
% %         [sortedStim,I] = sort([trial(eIdx).stimAmp]);
% %         
% %         if length(sortedStim) > 10
% %             iter = floor(length(sortedStim)/10);
% %         else 
% %             iter = 2;
% %         end
% %     
% %         h(p,f) = figure; maximize;
% %        
% %         for i = 1:length(msubset)
% %             m = msubset(i);
% %                  disp(['Processing ' Muscle_ID(m)]);
% %             for t = 1:iter:length(eIdx)
% %                 n = eIdx(I(t));
% %                 if t == 1
% %                     clear meanTrace
% %                 end
% %                 meanTrace(t,:) = abs(trial(n).meanTrace{m});
% %                 index = findchangepts(meanTrace, 'MaxNumChanges',6,'Statistic','rms'); %1 is stim, 2 is start, 3 is complex start, 4 is complex end, 5 is response end, 6 is return to baseline
% %                 if length(index) < 6
% %                     onset(t) = nan;
% %     %                 results(ii).muscle(c).offset(fnum) = nan;
% %                 else
% %                     onset(t) = trial(n).epochTimeVec(index(3));%2
% %     %                 onset_2(fnum) = trial(n).epochTimeVec(index(3));
% %     %                 results(ii).muscle(c).offset(fnum) = trial(n).epochTimeVec(index(5));
% %                 end
% % 
% %                 disp([num2str(t) ' of ' num2str(length(eIdx))]);
% %             end
% %             
% %             disp('Plotting Waterfall');
% %             figure(h(p,f));
% %             subplot(2,length(msubset)/2,i);
% % % %             figure(h(ii,3)); subplot(4,2,c-8);
% %             waterfall(meanTrace)
% %             xlim([0,6000])
% %             xticks(0:1500:6000)
% %             xticklabels({-50;0;50;150})
% %     %         xlabel('Time (ms)');
% %     %         zlabel('EMG (uV)');
% %             yticklabels(sortedStim(1:iter+2:end));
% %             yticks([1:iter+2:length(sortedStim)]);
% %             xaz =  -17.8667; xzen =  22.5799;
% %             view(xaz,xzen);
% %             title(Muscle_ID(m))
% %              
% %             
% %         end   
% %         figure(h(p,f));
% %         tit = suptitle(['E' num2str(trial(n).spinalElec) '- Waterfalls: ' num2str(check_f) 'Hz_' ... 
% %                  num2str(check_p) 'ms']);
% %         set(tit, 'Interpreter', 'none');
% %         
% %         
% %             disp(['Saving Full Trace ' num2str(n)]);
% %             saveas(h(p,f),[setPath '\Waterfalls\' 'E' num2str(trial(n).spinalElec) '- Waterfalls_' num2str(check_f) 'Hz_' ... 
% %                 num2str(check_p) 'ms.png']); %strrep(datestr(now),':','_')
% %             close(h(p,f));  
% %             
% %     end
% %      
% % end

%% Plot Rectified Traces of just first 4 stims
msubset = [11 12 16];
check_p = .500;
disp(['Pulse is ' num2str(check_p)]);
check_f = [1 2 10 20];
disp(['Freq is ' num2str(check_f)]);
check_amp = 5000;
nPre = .2*fs;
nPost = .2*fs;
hF = figure;   maximize;
for f = 1:length(check_f)
    eIdx = find(([trial.stimFrequency] == check_f(f)) & ([trial.pulseWidth] == check_p) & ([trial.stimAmp] == check_amp));
        disp(['Freq is ' num2str(check_f(f))]);
     %Account for multiple trials

      for t = 1%:length(eIdx)
          if isempty(eIdx)
              continue;
          end
        n = eIdx(t);
        disp(['Trial is ' num2str(n)]);
     
        disp(['Plotting Rect Full Trace']);

        disp(num2str([trial(n).stims(4)]));
        for i = 1:length(msubset)
            m = msubset(i);
            h(3*(f-1)+i) = subplot(length(check_f),length(msubset),3*(f-1)+i);
            timbegin = (trial(n).stims(1)-nPre):(trial(1).stims(4)+nPost);
            timend = (trial(n).stims(4)+nPost);
            plot(trial(n).timeVec(timbegin:timend), abs(trial(n).emgF(m,timbegin:timend)));                 
            hold on;
          vline(trial(n).stims(1:4)/fs,{'r','r','r','r'},{'1','2','3','4'});
            title([char(Muscle_ID(m)) ' - ' num2str(check_f(f)) 'Hz']);ylabel('EMG (uV)'); xlabel('Time (s)');                 
        end
      end 
end
linkaxes([h([1 4 7 10 13])],'xy'); linkaxes([h([2 5 8 11 14])],'y'); linkaxes([h([3 6 9 12 15])],'y');

% % suptitle(['E' num2str(elecs) ': 1st 4 Stims']);
% % saveas(hF,['D:\FigRescources\UH3\LSP05\rehash\PAD\E' num2str(elecs) '_stimtrace.png']);
% % savefig(hF,['D:\FigRescources\UH3\LSP05\rehash\PAD\E' num2str(elecs) '_stimtrace']);

%% Plot Rectified Traces of just first 4 stims for ONE muscle
msubset = [10 12 16];
check_p = .500;
disp(['Pulse is ' num2str(check_p)]);
check_f = [1 2 10 20];
disp(['Freq is ' num2str(check_f)]);
check_amp = 5000;
nPre = .02*fs;
nPost = .05*fs;
for mm = 1:length(msubset)
    hF(mm) = figure;   maximize;
    m = msubset(mm);
    for f = 1:length(check_f)
        eIdx = find(([trial.stimFrequency] == check_f(f)) & ([trial.pulseWidth] == check_p) & ([trial.stimAmp] == check_amp));
            disp(['Freq is ' num2str(check_f(f))]);
         %Account for multiple trials - Only take 1 trial for example

          for t = 1%:length(eIdx)
              if isempty(eIdx)
                  continue;
              end
            n = eIdx(t);
            disp(['Trial is ' num2str(n)]);

            disp(['Plotting Rect Full Trace']);

            disp(num2str([trial(n).stims(4)]));

                for i = 1:4 %No of Stims
                    h(4*(f-1)+i) = subplot(length(check_f),4,4*(f-1)+i);

                    timbegin = (trial(n).stims(i)-nPre);
                    timend = (trial(n).stims(i)+nPost);
                    tim = linspace(-nPre/fs, nPost/fs, nPre + nPost + 1);

                    plot(tim, abs(trial(n).emgF(m,timbegin:timend)));                 
                    hold on;
                   xline(0,'r',num2str(i));
                   xline(trial(n).stimLength/fs,'--r');
                    title([char(Muscle_ID(m)) ' - ' num2str(check_f(f)) 'Hz - stim ' num2str(i)]);
                    ylabel('EMG (uV)'); xlabel('Time (s)');                 
                end
          end
    end 
    linkaxes([h(:)],'xy');
    tit = sgtitle(['E' num2str(elecs) '_' char(Muscle_ID(m)) ': 1st 4 Stims']);
    set(tit,'Interpreter','none');
    saveas(hF(mm),['D:\FigRescources\UH3\LSP05\rehash\PAD\E' num2str(elecs) '_' char(Muscle_ID(m)) '_stimtrace.png']);
    savefig(hF(mm),['D:\FigRescources\UH3\LSP05\rehash\PAD\E' num2str(elecs) '_' char(Muscle_ID(m)) '_stimtrace']);

end

%% Get drop off
close all;

for f = 1:length(check_f)
    eIdx = find(([trial.stimFrequency] == check_f(f)) & ([trial.pulseWidth] == check_p) & ([trial.stimAmp] == check_amp));
        disp(['Freq is ' num2str(check_f(f))]);
     %Account for multiple trials
        if f > 10
            rWindow = 0.03*fs;
        else
            rWindow = .050*fs;
        end
      for t = 1%:length(eIdx)
        n = eIdx(t);
        disp(['Trial is ' num2str(n)]);
     
        disp(['Plotting Rect Full Trace']);

        disp(num2str([trial(n).stims(4)]));
        for i = 1:length(msubset)
            m = msubset(i);
            for s = 1:4
                timbegin = trial(n).stims(s)+(artifactBuffer*fs);
%                 response(i,f,s) = peak2peak(trial(n).emgF(m,timbegin:timbegin+rWindow));
                [peakR(i,f,s),onset(i,f,s)] = max(abs(trial(n).emgF(m,timbegin:timbegin+rWindow)));
                onset(i,f,s) = (onset(i,f,s)/fs + artifactBuffer)*1000;
                if onset(i,f,s) > 40
                    onset(i,f,s) = nan;
                    disp('No Response');
                end
            end                       
        end
    end    
end


%% Plot drop off
hF2 = figure; maximize;
C = linspecer(length(check_f));
 for i = 1:length(msubset)
            m = msubset(i);
            hs(i) = subplot(2,length(msubset),i);
            for f = 1:length(check_f)
                ref = peakR(i,f,1);
                for s = 1:4
                    tmp(f,s) = peakR(i,f,s)/ref;
                end
            end

            plot(tmp','LineWidth',2);
            hold on;
            ylabel('% of Initial Activation');
            title([Muscle_ID(m),{'Tracking PAD'}]);
            xlabel('stim order');
%             legend(arrayfun(@num2str,check_f,'UniformOutput',false),'Location','EastOutside','Orientation','vertical')
            xticks([1:4]);
            yticks([0:0.25:1]);

            
            
            subplot(2,length(msubset),3+i);
            tmp = reshape(onset(i,:,:),[length(check_f),4]);
            plot(tmp','LineWidth',2);  
            hold on;
            ylabel('Onset time (ms)')
%             title(['Response Onset (' char(Muscle_ID(m)) ')']);
            xlabel('stim order');
            xticks([1:4]);
            ylim([0,35]);
            yline(nanmean(mean(tmp)),'k','LineWidth',2);
            yline(nanmean(mean(tmp))+mean(nanstd(tmp)),'--k');
            yline(nanmean(mean(tmp))-mean(nanstd(tmp)),'--k');
            legtmp = arrayfun(@num2str,check_f,'UniformOutput',false);
            legtmp{end+1} = 'mean';
            legend(legtmp,'Location','NorthOutside','Orientation','horizontal')
            
 end
 set(0,'defaultAxesFontSize',12)
 linkaxes([hs(:)],'y');
 sgtitle(['Tracking Post-activation response (E' num2str(elecs) ',5mA)']);
 linkaxes([hs],'xy');
 saveas(hF2,['D:\FigRescources\UH3\LSP05\rehash\VL_BF_LG\E' num2str(elecs) '_PAD.png']);
 savefig(hF2,['D:\FigRescources\UH3\LSP05\rehash\VL_BF_LG\E' num2str(elecs) '_PAD']);
