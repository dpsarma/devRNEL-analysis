%% Plot Traces of just first 2 stims for subset of muscles
        Muscle_ID = {'Right VM', 'Right RF', 'Right VL', 'Right BF', 'Right ST', ...
        'Right TA', 'Right SO', 'Right LG','Left VM', 'Left RF','Left VL',...
        'Left BF', 'Left ST', 'Left TA', 'Left SO', 'Left LG'};
elecs = 1;
msubset = [10 12 16 14];
check_p = .200;
disp(['Pulse is ' num2str(check_p)]);
check_f = [1 2 10 20];
disp(['Freq is ' num2str(check_f)]);
check_amp = 4250;
fs = 7.5e3
nPre = .02*fs;
nPost = .05*fs;
    
for f = 1:length(check_f)

    artifactBuffer = 0.005;
    
    switch check_f(f)
        case 1
            dtPre = 50e-3; % msec before stim onset 
            dtPost = 100e-3; % msec after stim onset 
            rmsWindow = 0.05; % post stim RMS Window in s, so 0.05 is 50ms          
        case 2
            dtPre = 50e-3; % msec before stim onset 
            dtPost = 100e-3; % msec after stim onset 
            rmsWindow = 0.05; % post stim RMS Window in s, so 0.05 is 50ms
        case 5
            dtPre = 50e-3; % msec before stim onset 
            dtPost = 100e-3; % msec after stim onset 
            rmsWindow = 0.05; % post stim RMS Window in s, so 0.05 is 50ms
        case 10
            dtPre = 50e-3; % msec before stim onset 
            dtPost = 95e-3; % msec after stim onset 
            rmsWindow = 0.05; % post stim RMS Window in s, so 0.05 is 50ms
        case 20
            dtPre = 50e-3; % msec before stim onset 
            dtPost = 45e-3; % msec after stim onset 
            rmsWindow = 0.030; % post stim RMS Window in s, so 0.05 is 50ms
        case 200
            dtPre = 50e-3; % msec before stim onset 
            dtPost = 4e-3; % msec after stim onset 
            rmsWindow = 0.050; % post stim RMS Window in s, so 0.05 is 50ms
        otherwise 
            dtPre = 10e-3; % msec before stim onset 
            dtPost = 20e-3; % msec after stim onset 
            rmsWindow = 0.005; % post stim RMS Window in s, so 0.05 is 50ms;
    end
    eIdx = find(([trial.stimFrequency] == check_f(f)) & ([trial.pulseWidth] == check_p) & ([trial.stimAmp] == check_amp));
        disp(['Freq is ' num2str(check_f(f))]);
     %Account for multiple trials - Only take 1 trial for example
    hF(f) = figure;   %maximize;
    for mm = 1:length(msubset)
        m = msubset(mm);
        h(mm) = subplot(length(msubset),1, mm)

              for t = 1%:length(eIdx)
                  if isempty(eIdx)
                      continue;
                  end
                n = eIdx(t);
                disp(['Trial is ' num2str(n)]);

                disp(['Plotting Rect Full Trace']);

                disp(num2str([trial(n).stims(4)]));

    %                 for i = 1:4 %No of Stims
    %                     h(4*(f-1)+i) = subplot(length(check_f),4,4*(f-1)+i);

                        timbegin = (trial(n).stims(1)-nPre);
                        timend = (trial(n).stims(2)+nPost);
                        tim = linspace(-nPre/fs, nPost/fs, timend-timbegin + 1);

                        plot((trial(n).timeVec(timbegin:timend)), (trial(n).emgF(m,timbegin:timend)));                 
                        hold on;
                       xline(trial(n).stims(1)/fs,'r','1');
                       xline(trial(n).stims(2)/fs,'--r','2');
                       
                       xmin = min(trial(n).timeVec(timbegin:timend));
                       xmax = max(trial(n).timeVec(timbegin:timend));
                      xlim([xmin,xmax])
%                       xticklabels(0:0.04:.24)
    %                     title([char(Muscle_ID(m)) ' - ' num2str(check_f(f)) 'Hz - stim ' num2str(i)]);
                        ylabel({char(Muscle_ID(m));'EMG (uV)'});% ylim([-1200,1200]);
                        box(h(mm),'off');

    %                 end
              end

    end 
    xlabel('Time (s)')
    linkaxes([h(:)],'xy');
    tit = sgtitle(['E' num2str(elecs) '_' char(Muscle_ID(m)) ': 1st 4 Stims']);
    set(tit,'Interpreter','none');
    set(0,'defaultAxesFontSize',12)
% %     saveas(hF(f),['D:\FigRescources\UH3\LSP05\rehash\PAD\E' num2str(elecs) '_PAD_' num2str(check_f(f)) '_stimtrace.png']);
% %     saveas(hF(f),['D:\FigRescources\UH3\LSP05\rehash\PAD\E' num2str(elecs) '_PAD_' num2str(check_f(f)) '_stimtrace.svg']);
% %     % %     saveas(hF(mm),['D:\FigRescources\UH3\LSP05\rehash\PAD\E' num2str(elecs) '_' char(Muscle_ID(m)) '_stimtrace.png']);
% %     savefig(hF(mm),['D:\FigRescources\UH3\LSP05\rehash\PAD\E' num2str(elecs) '_' char(Muscle_ID(m)) '_stimtrace']);

end


%% Get drop off for 20Hz
% % close all;

for f = 1:length(check_f)
    eIdx = find(([trial.stimFrequency] == check_f(f)) & ([trial.pulseWidth] == check_p) & ([trial.stimAmp] == check_amp));
        disp(['Freq is ' num2str(check_f(f))]);
     %Account for multiple trials
        if f > 10
            rWindow = 0.03*fs;
        else
            rWindow = .050*fs;
        end
      for t = 1:length(eIdx)
        n = eIdx(t);
        disp(['Trial is ' num2str(n)]);
     
        disp(['Plotting Rect Full Trace']);

        disp(num2str([trial(n).stims(4)]));
        for i = 1:length(msubset)
            m = msubset(i);
            for s = 1:4
                timbegin = trial(n).stims(s)+(artifactBuffer*fs);
                idx = s+2*(t-1);
%                 response(i,f,s) = peak2peak(trial(n).emgF(m,timbegin:timbegin+rWindow));
                [peakR(i,f,idx),onset(i,f,idx)] = max(abs(trial(n).emgF(m,timbegin:timbegin+rWindow)));
                onset(i,f,idx) = (onset(i,f,idx)/fs + artifactBuffer)*1000;
                if onset(i,f,idx) > 40
                    onset(i,f,idx) = nan;
                    disp('No Response');
                end
            end                       
        end
    end    
end