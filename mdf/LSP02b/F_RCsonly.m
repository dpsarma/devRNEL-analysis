% elec = [9 10 11 12];
% fnum = 1;
omdfc = mdf.init('auto');

summary =  mdf.load('mdf_type','summary','subject','LSP02b');
musclabels = summary.emgLabel;

electrodenums = [1:48];
datestrings = ['2018-10-08'; '2018-10-09'; '2018-10-11'; '2018-10-16'; '2018-10-19'; '2018-10-24'];
dayIDs = ['Day 02'; 'Day 03'; 'Day 04'; 'Day 07'; 'Day 09'; 'Day 12'];

datasave = 'R:\users\dsarma\LSP02b\matData\mdftest\';
figuresave = 'R:\users\dsarma\LSP02b\Figures\mdftest_RCs\';

for i=1:length(electrodenums)
    if (electrodenums(i)<10)
        electrodestrings{i} = ['Unipolar: 0' num2str(electrodenums(i))];
    else
        electrodestrings{i} = ['Unipolar: ' num2str(electrodenums(i))];
    end 
end

filenames = {'lBF_AllDays_AllElectrodes_1Hz_200us'; 'lST_AllDays_AllElectrodes_1Hz_200us';...
'lLG_AllDays_AllElectrodes_1Hz_200us'; 'lRF_AllDays_AllElectrodes_1Hz_200us';...
'lSO_AllDays_AllElectrodes_1Hz_200us'; 'lVL_AllDays_AllElectrodes_1Hz_200us';...
'lVM_AllDays_AllElectrodes_1Hz_200us'; 'rLG_AllDays_AllElectrodes_1Hz_200us';...
'rBF_AllDays_AllElectrodes_1Hz_200us'; 'rST_AllDays_AllElectrodes_1Hz_200us';...
'rTA_AllDays_AllElectrodes_1Hz_200us'; 'rVL_AllDays_AllElectrodes_1Hz_200us';...
'rRF_AllDays_AllElectrodes_1Hz_200us'; 'rSO_AllDays_AllElectrodes_1Hz_200us'; ...
'rVM_AllDays_AllElectrodes_1Hz_200us'};

% filenames = {'lBF_AllDays_AllElectrodes_1Hz_200us'; 'lST_AllDays_AllElectrodes_1Hz_200us';...
% 'rBF_AllDays_AllElectrodes_1Hz_200us'; 'rST_AllDays_AllElectrodes_1Hz_200us'};


% filenames = {'lBF_AllDays_AllElectrodes_1Hz_200us';};


for k2 = 1:length(filenames)
    C = strsplit(char(filenames(k2)),'_');
    mnum =  find(strcmp(musclabels,C(1)));
    freqstring = str2double(string(regexp(C(4),'[0-9]','match')));
    pulsestring = 200; %str2double(string(erase(C(5),'us.mat')));
    
    
    load(char(filenames(k2)));

%% Plot RC (useless)
    for m = mnum%1:lengthmusclabels
        hRC(m) = figure; maximize; plotcnt = 1;
        hold on;
        for f = 1:size(hash,2) %fnum
            X = vertcat(hash{:,f,:});
            maxP2P = max([X.meanP2P]);
            for e = 1:size(hash,1) %elec
   
                S = vertcat(hash{e,f,:});
                if isempty(S)
                   disp([dayIDs(f,:) ', ' char(electrodestrings(e)) ', ' num2str(freqstring) 'Hz, ' num2str(pulsestring) 'us']);
                   disp('Moving to Next Electrode')
                    continue;
                end
                allP2P = [S.meanP2P];
                allP2Perr = [S.errP2P];
                
                allStims = [S.stimAmp];
                allfreqs  = [S.stimfreq];
                allpulses  = [S.stimpulse];
                
                plotcnt = plotcnt+1;
                lgndentry{plotcnt} = [dayIDs(f,:) ', ' char(electrodestrings(e)) ', ' num2str(freqstring) 'Hz, ' num2str(pulsestring) 'us'];
                
                  if e==11
                    allP2P(allStims == 6) = [];
                    allP2Perr(allStims == 6) = [];
                    allStims(allStims == 6) = [];
                end
                
                [sortedStim,I] = sort(allStims);
                
                allP2Perr = allP2Perr/maxP2P;
                allP2P = allP2P/maxP2P;
                
%                 stimAlls{plotcnt} = sortedStim;
%                 P2PAll{plotcnt} = allP2P(I);
%                 P2PerrAll{plotcnt} = allP2Perr(I);
                
                yyA = smooth(sortedStim,allP2P(I),'sgolay');
                plot(sortedStim, yyA, 'LineWidth',2);
%                 errorbar(sortedStim,yyA,allP2Perr(I),'LineWidth',6);
                
            end
        end
        ylabel('(post-Stim) Mean P2P  (uV)')
        xlabel('Stimulation Amplitude (mA)')
        title([ char(musclabels(m)) ' - Amplitude Recruitment Curves - '  dayIDs(f,:)], 'Interpreter', 'none')
        legend(lgndentry(~cellfun('isempty',lgndentry')),'Location', 'eastoutside');
        ylim([0 1.1]);
        savefig(hRC(m),[figuresave char(musclabels(m)) '_P2P_Recruitment_' num2str(freqstring) 'Hz_' num2str(pulsestring) 'us']);
%         saveas(hRC(m),[figuresave char(musclabels(m)) '_P2P_Recruitment' num2str(freqstring) 'Hz_' num2str(pulsestring) 'us.png']);
    end
end 