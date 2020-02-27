%% UH3 - LSP02 anaylsis (2) - MDF Recruitment Curve RMS
%----------------------------------------
% tag that is missing from all EMG trials is "Sitting" vs. "Standing"

%% Build Recruitment Curves / Muscle / electrode
figuresave = 'D:\FigRescources\UH3\LSP02\BRAIN meeting 4-19\';
% % e = 10; f = 1; m = 13; tnum = 10;

for m = mnum%1:lengthmusclabels
    hRC(m) = figure;
    hold on;
    for f = fnum; %1:size(datestrings,1)
        for e = elec; %1:length(electrodenums)
            for tnum = 1:length(uuid_list{e,f})

            allStims(tnum) = hash{m,e,f,tnum}.stimAmp;
            allfreqs(tnum) = hash{m,e,f,tnum}.stimfreq;
            allpulses(tnum) = hash{m,e,f,tnum}.stimpulse;
            
            allRMS(tnum) = hash{m,e,f,tnum}.meanRMS;
            allRMSerr(tnum) = hash{m,e,f,tnum}.errRMS;

            
            end
            lgndentry{e} = char(electrodestrings(e));
            [sortedStim,I] = sort(allStims);
            yyA = smooth(sortedStim,allRMS(I),'sgolay');
            errorbar(sortedStim,yyA,allRMSerr(I),'LineWidth',2);
            
        end
    end
    ylabel('(post-Stim) Mean RMS  (uV)')
    xlabel('Stimulation Amplitude (mA)')
    title([char(musclabels(m)) ' ' dayIDs(f,:) 'Amplitude Recruitment Curves'])
    legend(lgndentry(~cellfun('isempty',lgndentry')),'Location', 'northwest');
    savefig(hRC(m),[figuresave char(musclabels(m)) '_RMS_Recruitment_' num2str(freqstring) 'Hz_' num2str(pulsestring) 'us']);
    saveas(hRC(m),[figuresave char(musclabels(m)) '_RMS_Recruitment' num2str(freqstring) 'Hz_' num2str(pulsestring) 'us.png']);
end
