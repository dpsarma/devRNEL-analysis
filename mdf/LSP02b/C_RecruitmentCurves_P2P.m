%% UH3 - LSP02 anaylsis (2b) - MDF Recruitment Curve P2P
%----------------------------------------
% tag that is missing from all EMG trials is "Sitting" vs. "Standing"

%% Build Recruitment Curves / Muscle / electrode
figuresave = 'D:\FigRescources\UH3\LSP02\BRAIN meeting 4-19\';
elec = 9; f = 5; m = 13; tnum = 10;

for m = mnum%1:lengthmusclabels
    hRC(m) = figure;
    hold on;
    for f = fnum; %1:size(datestrings,1)
        for e = elec; %1:length(electrodenums)
            for tnum = 1:length(uuid_list{e,f})

            allStims(tnum) = hash{e,f,tnum}.stimAmp;
            allfreqs(tnum) = hash{e,f,tnum}.stimfreq;
            allpulses(tnum) = hash{e,f,tnum}.stimpulse;
            
            allP2P(tnum) = hash{e,f,tnum}.meanP2P;
            allP2Perr(tnum) = hash{e,f,tnum}.errP2P;
            
            end
            lgndentry{e} = char(electrodestrings(e));
            [sortedStim,I] = sort(allStims);
            yyA = smooth(sortedStim,allP2P(I),'sgolay');
            errorbar(sortedStim,yyA,allP2Perr(I),'LineWidth',2);
            
        end
    end
    ylabel('(post-Stim) Mean P2P  (uV)')
    xlabel('Stimulation Amplitude (mA)')
    title([char(musclabels(m)) ' ' dayIDs(f,:) 'Amplitude Recruitment Curves'])
    legend(lgndentry(~cellfun('isempty',lgndentry')),'Location', 'northwest');
    savefig(hRC(m),[figuresave char(musclabels(m)) '_P2P_Recruitment_' num2str(freqstring) 'Hz_' num2str(pulsestring) 'us']);
    saveas(hRC(m),[figuresave char(musclabels(m)) '_P2P_Recruitment' num2str(freqstring) 'Hz_' num2str(pulsestring) 'us.png']);
end
