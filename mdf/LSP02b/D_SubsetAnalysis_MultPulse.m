%% UH3 - LSP02 anaylsis (3a) - MDF Analysis
%----------------------------------------
% tag that is missing from all EMG trials is "Sitting" vs. "Standing"

%% Load for subset of parameters

datasave = 'R:\users\dsarma\LSP02b\matData\mdftest\';
figuresave = 'D:\FigRescources\UH3\LSP02\BRAIN meeting 4-19\';

summary =  mdf.load('mdf_type','summary');
musclabels = summary.emgLabel;

datestrings = ['2018-10-08'; '2018-10-09'; '2018-10-11'; '2018-10-16'; '2018-10-19'; '2018-10-24'];
dayIDs = ['Day 02'; 'Day 03'; 'Day 04'; 'Day 07'; 'Day 09'; 'Day 12'];
freqs = [1 2 5];
pulses = [200 400 600 800 1000];

electrodenums = [1:48];
for i=1:length(electrodenums)
            if (electrodenums(i)<10)
                electrodestrings{i} = ['Unipolar: 0' num2str(electrodenums(i))];
            else
                electrodestrings{i} = ['Unipolar: ' num2str(electrodenums(i))];
            end 
end


listing = dir('R:\users\dsarma\LSP02b\matData\mdftest\Day 09');
listing(1:2) = [];
for k=1:length(listing)
    filenames{k} = listing(k).name;
    
end

% Peak Window
gWind = 50;
fs =  30e3;
dtPre = 50e-3; % msec before stim onset
dtPost= 5e-3; %msec after stim onset
nSampsPre = floor((dtPre-5e-3)*fs);
nSampsPost = floor((dtPost+dtPre)*fs);
respSamps = floor((45e-3)*fs);

for k2 = 1:length(filenames)
    
    C = strsplit(char(filenames(k2)),'_');

    mnum =  find(strcmp(musclabels,C(1)));
    fnum = find(strcmp(dayIDs,C(2)));
    elec = str2double(string(regexp(C(3),'[0-9]','match')));
    freqstring = str2double(string(regexp(C(4),'[0-9]','match')));
    pulsestring = str2double(string(erase(C(5),'us.mat')));

    disp([num2str(freqstring) 'Hz - ' num2str(pulsestring) 'us']);
    if (freqstring == 1 || freqstring == 5)
        disp('Moving On');
        continue;
    end
    
    load(char(filenames(k2)));
    for m = mnum%1:lengthmusclabels
            for f = fnum; %1:size(datestrings,1)
                for e = elec; %1:length(electrodenums)
                    for tnum = 1:length(hash(mnum,elec,fnum,:))
                        if isempty(hash{m,e,f,tnum})
                            continue;
                        end
    
    
                        switch elec
                            case 2
                                if pulsestring == 200
                                     E2_PW200(m,tnum).Stims = hash{m,e,f,tnum}.stimAmp;
                                     E2_PW200(m,tnum).P2P = hash{m,e,f,tnum}.meanP2P;
                                     E2_PW200(m,tnum).P2Perr = hash{m,e,f,tnum}.errP2P;
                                     E2_PW200(m,tnum).meantrace = mean(hash{m,e,f,tnum}.epochsrect);
                                     E2_PW200(m,tnum).epochtime = hash{m,e,f,tnum}.epochtime;


                                     for i = 1:size(hash{m,e,f,tnum}.epochsrect,1)
                                         [pks,locs,widths,proms] = findpeaks(smoothdata(hash{m,e,f,tnum}.epochsrect(nSampsPost:nSampsPost+respSamps)','gaussian',gWind)',...
                                        E2_PW200(m,tnum).epochtime(nSampsPost:nSampsPost+respSamps), ...
                                        'Annotate','extents', 'WidthReference','halfheight'); 
                                        [maxpeak(i), ipeak] = max(pks);
                                        timeonset(i) = locs(ipeak);
                                     end
                                        E2_PW200(m,tnum).pkVar = std(maxpeak)/sqrt(length(maxpeak));
                                        E2_PW200(m,tnum).timeOnsVar = std(timeonset)/sqrt(length(timeonset));

                                     [pks,locs,widths,proms] = findpeaks(smoothdata(E2_PW200(m,tnum).meantrace(nSampsPost:nSampsPost+respSamps)','gaussian',gWind)',...
                                        E2_PW200(m,tnum).epochtime(nSampsPost:nSampsPost+respSamps), ...
                                        'Annotate','extents', 'WidthReference','halfheight'); 
                                     [E2_PW200(m,tnum).maxMeanpeak, ipeak] = max(pks);
                                     E2_PW200(m,tnum).timeMeanpeak = locs(ipeak);


                                elseif pulsestring == 400
                                     E2_PW400(m,tnum).Stims = hash{m,e,f,tnum}.stimAmp;
                                     E2_PW400(m,tnum).P2P = hash{m,e,f,tnum}.meanP2P;
                                     E2_PW400(m,tnum).P2Perr = hash{m,e,f,tnum}.errP2P;
                                     E2_PW400(m,tnum).meantrace = mean(hash{m,e,f,tnum}.epochsrect);
                                     E2_PW400(m,tnum).epochtime = hash{m,e,f,tnum}.epochtime;


                                     for i = 1:size(hash{m,e,f,tnum}.epochsrect,1)
                                         [pks,locs,widths,proms] = findpeaks(smoothdata(hash{m,e,f,tnum}.epochsrect(nSampsPost:nSampsPost+respSamps)','gaussian',gWind)',...
                                        E2_PW400(m,tnum).epochtime(nSampsPost:nSampsPost+respSamps), ...
                                        'Annotate','extents', 'WidthReference','halfheight'); 
                                        [maxpeak(i), ipeak] = max(pks);
                                        timeonset(i) = locs(ipeak);
                                     end
                                        E2_PW400(m,tnum).pkVar = std(maxpeak)/sqrt(length(maxpeak));
                                        E2_PW400(m,tnum).timeOnsVar = std(timeonset)/sqrt(length(timeonset));

                                     [pks,locs,widths,proms] = findpeaks(smoothdata(E2_PW400(m,tnum).meantrace(nSampsPost:nSampsPost+respSamps)','gaussian',gWind)',...
                                        E2_PW400(m,tnum).epochtime(nSampsPost:nSampsPost+respSamps), ...
                                        'Annotate','extents', 'WidthReference','halfheight'); 
                                     [E2_PW400(m,tnum).maxMeanpeak, ipeak] = max(pks);
                                     E2_PW400(m,tnum).timeMeanpeak = locs(ipeak);

                                elseif pulsestring == 600
                                                     E2_PW600(m,tnum).Stims = hash{m,e,f,tnum}.stimAmp;
                                     E2_PW600(m,tnum).P2P = hash{m,e,f,tnum}.meanP2P;
                                     E2_PW600(m,tnum).P2Perr = hash{m,e,f,tnum}.errP2P;
                                     E2_PW600(m,tnum).meantrace = mean(hash{m,e,f,tnum}.epochsrect);
                                     E2_PW600(m,tnum).epochtime = hash{m,e,f,tnum}.epochtime;


                                     for i = 1:size(hash{m,e,f,tnum}.epochsrect,1)
                                         [pks,locs,widths,proms] = findpeaks(smoothdata(hash{m,e,f,tnum}.epochsrect(nSampsPost:nSampsPost+respSamps)','gaussian',gWind)',...
                                        E2_PW600(m,tnum).epochtime(nSampsPost:nSampsPost+respSamps), ...
                                        'Annotate','extents', 'WidthReference','halfheight'); 
                                        [maxpeak(i), ipeak] = max(pks);
                                        timeonset(i) = locs(ipeak);
                                     end
                                        E2_PW600(m,tnum).pkVar = std(maxpeak)/sqrt(length(maxpeak));
                                        E2_PW600(m,tnum).timeOnsVar = std(timeonset)/sqrt(length(timeonset));

                                     [pks,locs,widths,proms] = findpeaks(smoothdata(E2_PW600(m,tnum).meantrace(nSampsPost:nSampsPost+respSamps)','gaussian',gWind)',...
                                        E2_PW600(m,tnum).epochtime(nSampsPost:nSampsPost+respSamps), ...
                                        'Annotate','extents', 'WidthReference','halfheight'); 
                                     [E2_PW600(m,tnum).maxMeanpeak, ipeak] = max(pks);
                                     E2_PW600(m,tnum).timeMeanpeak = locs(ipeak);

                                elseif pulsestring == 800

                                     E2_PW800(m,tnum).Stims = hash{m,e,f,tnum}.stimAmp;
                                     E2_PW800(m,tnum).P2P = hash{m,e,f,tnum}.meanP2P;
                                     E2_PW800(m,tnum).P2Perr = hash{m,e,f,tnum}.errP2P;
                                     E2_PW800(m,tnum).meantrace = mean(hash{m,e,f,tnum}.epochsrect);
                                     E2_PW800(m,tnum).epochtime = hash{m,e,f,tnum}.epochtime;


                                     for i = 1:size(hash{m,e,f,tnum}.epochsrect,1)
                                         [pks,locs,widths,proms] = findpeaks(smoothdata(hash{m,e,f,tnum}.epochsrect(nSampsPost:nSampsPost+respSamps)','gaussian',gWind)',...
                                        E2_PW800(m,tnum).epochtime(nSampsPost:nSampsPost+respSamps), ...
                                        'Annotate','extents', 'WidthReference','halfheight'); 
                                        [maxpeak(i), ipeak] = max(pks);
                                        timeonset(i) = locs(ipeak);
                                     end
                                        E2_PW800(m,tnum).pkVar = std(maxpeak)/sqrt(length(maxpeak));
                                        E2_PW800(m,tnum).timeOnsVar = std(timeonset)/sqrt(length(timeonset));

                                     [pks,locs,widths,proms] = findpeaks(smoothdata(E2_PW800(m,tnum).meantrace(nSampsPost:nSampsPost+respSamps)','gaussian',gWind)',...
                                        E2_PW800(m,tnum).epochtime(nSampsPost:nSampsPost+respSamps), ...
                                        'Annotate','extents', 'WidthReference','halfheight'); 
                                     [E2_PW800(m,tnum).maxMeanpeak, ipeak] = max(pks);
                                     E2_PW800(m,tnum).timeMeanpeak = locs(ipeak);

                               elseif pulsestring == 1000

                                     E2_PW1000(m,tnum).Stims = hash{m,e,f,tnum}.stimAmp;
                                     E2_PW1000(m,tnum).P2P = hash{m,e,f,tnum}.meanP2P;
                                     E2_PW1000(m,tnum).P2Perr = hash{m,e,f,tnum}.errP2P;
                                     E2_PW1000(m,tnum).meantrace = mean(hash{m,e,f,tnum}.epochsrect);
                                     E2_PW1000(m,tnum).epochtime = hash{m,e,f,tnum}.epochtime;


                                     for i = 1:size(hash{m,e,f,tnum}.epochsrect,1)
                                         [pks,locs,widths,proms] = findpeaks(smoothdata(hash{m,e,f,tnum}.epochsrect(nSampsPost:nSampsPost+respSamps)','gaussian',gWind)',...
                                        E2_PW1000(m,tnum).epochtime(nSampsPost:nSampsPost+respSamps), ...
                                        'Annotate','extents', 'WidthReference','halfheight'); 
                                        [maxpeak(i), ipeak] = max(pks);
                                        timeonset(i) = locs(ipeak);
                                     end
                                        E2_PW1000(m,tnum).pkVar = std(maxpeak)/sqrt(length(maxpeak));
                                        E2_PW1000(m,tnum).timeOnsVar = std(timeonset)/sqrt(length(timeonset));

                                     [pks,locs,widths,proms] = findpeaks(smoothdata(E2_PW1000(m,tnum).meantrace(nSampsPost:nSampsPost+respSamps)','gaussian',gWind)',...
                                        E2_PW1000(m,tnum).epochtime(nSampsPost:nSampsPost+respSamps), ...
                                        'Annotate','extents', 'WidthReference','halfheight'); 
                                     [E2_PW1000(m,tnum).maxMeanpeak, ipeak] = max(pks);
                                     E2_PW1000(m,tnum).timeMeanpeak = locs(ipeak);  
                                else
                                    continue;
                                end
                                
                            case 4
                                if pulsestring == 200
                                     E4_PW200(m,tnum).Stims = hash{m,e,f,tnum}.stimAmp;
                                     E4_PW200(m,tnum).P2P = hash{m,e,f,tnum}.meanP2P;
                                     E4_PW200(m,tnum).P2Perr = hash{m,e,f,tnum}.errP2P;
                                     E4_PW200(m,tnum).meantrace = mean(hash{m,e,f,tnum}.epochsrect);
                                     E4_PW200(m,tnum).epochtime = hash{m,e,f,tnum}.epochtime;


                                     for i = 1:size(hash{m,e,f,tnum}.epochsrect,1)
                                         [pks,locs,widths,proms] = findpeaks(smoothdata(hash{m,e,f,tnum}.epochsrect(nSampsPost:nSampsPost+respSamps)','gaussian',gWind)',...
                                        E4_PW200(m,tnum).epochtime(nSampsPost:nSampsPost+respSamps), ...
                                        'Annotate','extents', 'WidthReference','halfheight'); 
                                        [maxpeak(i), ipeak] = max(pks);
                                        timeonset(i) = locs(ipeak);
                                     end
                                        E4_PW200(m,tnum).pkVar = std(maxpeak)/sqrt(length(maxpeak));
                                        E4_PW200(m,tnum).timeOnsVar = std(timeonset)/sqrt(length(timeonset));

                                     [pks,locs,widths,proms] = findpeaks(smoothdata(E4_PW200(m,tnum).meantrace(nSampsPost:nSampsPost+respSamps)','gaussian',gWind)',...
                                        E4_PW200(m,tnum).epochtime(nSampsPost:nSampsPost+respSamps), ...
                                        'Annotate','extents', 'WidthReference','halfheight'); 
                                     [E4_PW200(m,tnum).maxMeanpeak, ipeak] = max(pks);
                                     E4_PW200(m,tnum).timeMeanpeak = locs(ipeak);


                                elseif pulsestring == 400
                                     E4_PW400(m,tnum).Stims = hash{m,e,f,tnum}.stimAmp;
                                     E4_PW400(m,tnum).P2P = hash{m,e,f,tnum}.meanP2P;
                                     E4_PW400(m,tnum).P2Perr = hash{m,e,f,tnum}.errP2P;
                                     E4_PW400(m,tnum).meantrace = mean(hash{m,e,f,tnum}.epochsrect);
                                     E4_PW400(m,tnum).epochtime = hash{m,e,f,tnum}.epochtime;


                                     for i = 1:size(hash{m,e,f,tnum}.epochsrect,1)
                                         [pks,locs,widths,proms] = findpeaks(smoothdata(hash{m,e,f,tnum}.epochsrect(nSampsPost:nSampsPost+respSamps)','gaussian',gWind)',...
                                        E4_PW400(m,tnum).epochtime(nSampsPost:nSampsPost+respSamps), ...
                                        'Annotate','extents', 'WidthReference','halfheight'); 
                                        [maxpeak(i), ipeak] = max(pks);
                                        timeonset(i) = locs(ipeak);
                                     end
                                        E4_PW400(m,tnum).pkVar = std(maxpeak)/sqrt(length(maxpeak));
                                        E4_PW400(m,tnum).timeOnsVar = std(timeonset)/sqrt(length(timeonset));

                                     [pks,locs,widths,proms] = findpeaks(smoothdata(E4_PW400(m,tnum).meantrace(nSampsPost:nSampsPost+respSamps)','gaussian',gWind)',...
                                        E4_PW400(m,tnum).epochtime(nSampsPost:nSampsPost+respSamps), ...
                                        'Annotate','extents', 'WidthReference','halfheight'); 
                                     [E4_PW400(m,tnum).maxMeanpeak, ipeak] = max(pks);
                                     E4_PW400(m,tnum).timeMeanpeak = locs(ipeak);

                                elseif pulsestring == 600
                                                     E4_PW600(m,tnum).Stims = hash{m,e,f,tnum}.stimAmp;
                                     E4_PW600(m,tnum).P2P = hash{m,e,f,tnum}.meanP2P;
                                     E4_PW600(m,tnum).P2Perr = hash{m,e,f,tnum}.errP2P;
                                     E4_PW600(m,tnum).meantrace = mean(hash{m,e,f,tnum}.epochsrect);
                                     E4_PW600(m,tnum).epochtime = hash{m,e,f,tnum}.epochtime;


                                     for i = 1:size(hash{m,e,f,tnum}.epochsrect,1)
                                         [pks,locs,widths,proms] = findpeaks(smoothdata(hash{m,e,f,tnum}.epochsrect(nSampsPost:nSampsPost+respSamps)','gaussian',gWind)',...
                                        E4_PW600(m,tnum).epochtime(nSampsPost:nSampsPost+respSamps), ...
                                        'Annotate','extents', 'WidthReference','halfheight'); 
                                        [maxpeak(i), ipeak] = max(pks);
                                        timeonset(i) = locs(ipeak);
                                     end
                                        E4_PW600(m,tnum).pkVar = std(maxpeak)/sqrt(length(maxpeak));
                                        E4_PW600(m,tnum).timeOnsVar = std(timeonset)/sqrt(length(timeonset));

                                     [pks,locs,widths,proms] = findpeaks(smoothdata(E4_PW600(m,tnum).meantrace(nSampsPost:nSampsPost+respSamps)','gaussian',gWind)',...
                                        E4_PW600(m,tnum).epochtime(nSampsPost:nSampsPost+respSamps), ...
                                        'Annotate','extents', 'WidthReference','halfheight'); 
                                     [E4_PW600(m,tnum).maxMeanpeak, ipeak] = max(pks);
                                     E4_PW600(m,tnum).timeMeanpeak = locs(ipeak);

                                elseif pulsestring == 800

                                     E4_PW800(m,tnum).Stims = hash{m,e,f,tnum}.stimAmp;
                                     E4_PW800(m,tnum).P2P = hash{m,e,f,tnum}.meanP2P;
                                     E4_PW800(m,tnum).P2Perr = hash{m,e,f,tnum}.errP2P;
                                     E4_PW800(m,tnum).meantrace = mean(hash{m,e,f,tnum}.epochsrect);
                                     E4_PW800(m,tnum).epochtime = hash{m,e,f,tnum}.epochtime;


                                     for i = 1:size(hash{m,e,f,tnum}.epochsrect,1)
                                         [pks,locs,widths,proms] = findpeaks(smoothdata(hash{m,e,f,tnum}.epochsrect(nSampsPost:nSampsPost+respSamps)','gaussian',gWind)',...
                                        E4_PW800(m,tnum).epochtime(nSampsPost:nSampsPost+respSamps), ...
                                        'Annotate','extents', 'WidthReference','halfheight'); 
                                        [maxpeak(i), ipeak] = max(pks);
                                        timeonset(i) = locs(ipeak);
                                     end
                                        E4_PW800(m,tnum).pkVar = std(maxpeak)/sqrt(length(maxpeak));
                                        E4_PW800(m,tnum).timeOnsVar = std(timeonset)/sqrt(length(timeonset));

                                     [pks,locs,widths,proms] = findpeaks(smoothdata(E4_PW800(m,tnum).meantrace(nSampsPost:nSampsPost+respSamps)','gaussian',gWind)',...
                                        E4_PW800(m,tnum).epochtime(nSampsPost:nSampsPost+respSamps), ...
                                        'Annotate','extents', 'WidthReference','halfheight'); 
                                     [E4_PW800(m,tnum).maxMeanpeak, ipeak] = max(pks);
                                     E4_PW800(m,tnum).timeMeanpeak = locs(ipeak);

                               elseif pulsestring == 1000

                                     E4_PW1000(m,tnum).Stims = hash{m,e,f,tnum}.stimAmp;
                                     E4_PW1000(m,tnum).P2P = hash{m,e,f,tnum}.meanP2P;
                                     E4_PW1000(m,tnum).P2Perr = hash{m,e,f,tnum}.errP2P;
                                     E4_PW1000(m,tnum).meantrace = mean(hash{m,e,f,tnum}.epochsrect);
                                     E4_PW1000(m,tnum).epochtime = hash{m,e,f,tnum}.epochtime;


                                     for i = 1:size(hash{m,e,f,tnum}.epochsrect,1)
                                         [pks,locs,widths,proms] = findpeaks(smoothdata(hash{m,e,f,tnum}.epochsrect(nSampsPost:nSampsPost+respSamps)','gaussian',gWind)',...
                                        E4_PW1000(m,tnum).epochtime(nSampsPost:nSampsPost+respSamps), ...
                                        'Annotate','extents', 'WidthReference','halfheight'); 
                                        [maxpeak(i), ipeak] = max(pks);
                                        timeonset(i) = locs(ipeak);
                                     end
                                        E4_PW1000(m,tnum).pkVar = std(maxpeak)/sqrt(length(maxpeak));
                                        E4_PW1000(m,tnum).timeOnsVar = std(timeonset)/sqrt(length(timeonset));

                                     [pks,locs,widths,proms] = findpeaks(smoothdata(E4_PW1000(m,tnum).meantrace(nSampsPost:nSampsPost+respSamps)','gaussian',gWind)',...
                                        E4_PW1000(m,tnum).epochtime(nSampsPost:nSampsPost+respSamps), ...
                                        'Annotate','extents', 'WidthReference','halfheight'); 
                                     [E4_PW1000(m,tnum).maxMeanpeak, ipeak] = max(pks);
                                     E4_PW1000(m,tnum).timeMeanpeak = locs(ipeak);  
                                else
                                    continue;
                                end
                        end
        
        
% % %Plot RC (useless)
% %     for m = mnum%1:lengthmusclabels
% %         hRC(m) = figure;
% %         hold on;
% %         for f = fnum; %1:size(datestrings,1)
% %             for e = elec; %1:length(electrodenums)
% %                 for tnum = 1:length(hash(mnum,elec,fnum,:))
% %                     if isempty(hash{m,e,f,tnum})
% %                         continue;
% %                     end
% % 
% %                 allStims(tnum) = hash{m,e,f,tnum}.stimAmp;
% %                 allfreqs(tnum) = hash{m,e,f,tnum}.stimfreq;
% %                 allpulses(tnum) = hash{m,e,f,tnum}.stimpulse;
% % 
% %                 allP2P(tnum) = hash{m,e,f,tnum}.meanP2P;
% %                 allP2Perr(tnum) = hash{m,e,f,tnum}.errP2P;
% % 
% %                 end
% %                 lgndentry{e} = [char(electrodestrings(e)) ', ' num2str(freqstring) 'Hz, ' num2str(pulsestring) 'us'];
% %                 [sortedStim,I] = sort(allStims);
% %                 yyA = smooth(sortedStim,allP2P(I),'sgolay');
% %                 errorbar(sortedStim,yyA,allP2Perr(I),'LineWidth',2);
% % 
% %             end
% %         end
% %         ylabel('(post-Stim) Mean P2P  (uV)')
% %         xlabel('Stimulation Amplitude (mA)')
% %         title([ char(musclabels(m)) ' - Amplitude Recruitment Curves - '  dayIDs(f,:)], 'Interpreter', 'none')
% %         legend(lgndentry(~cellfun('isempty',lgndentry')),'Location', 'northwest');
% %         savefig(hRC(m),[figuresave char(musclabels(m)) '_P2P_Recruitment_' num2str(freqstring) 'Hz_' num2str(pulsestring) 'us']);
% %         saveas(hRC(m),[figuresave char(musclabels(m)) '_P2P_Recruitment' num2str(freqstring) 'Hz_' num2str(pulsestring) 'us.png']);
% %     end
% %     
    
                    end
                end
            end
    end
end
% 
% subMusc = [4 5 13 14];
% 
% caseNums = length(pulses);
% 
% for m = [4 5 13 14]
% %     for pw =  pulses
% 
%             stims = arrayfun(@(c) E2_PW1000(m,c).Stims, (1:length((E2_PW1000(m,:)))));
%             means = arrayfun(@(c) E2_PW1000(m,c).P2P, (1:length((E2_PW1000(m,:)))));
%             pkVar = arrayfun(@(c) E2_PW1000(14,c).P2Perr, (1:length((E2_PW1000(m,:)))));
%             timeOns = arrayfun(@(c) E2_PW1000(14,c).timeMeanpeak, (1:length((E2_PW1000(m,:))))) + dtPost
%             timeOnsVar = arrayfun(@(c) E2_PW1000(14,c).timeOnsVar, (1:length((E2_PW1000(m,:)))))
%             figure;
%             subplot(211)
%             br1 = bar(stims,means'); hold on;
%             er1 = errorbar(stims, means,pkVar,'k', 'LineWidth', 2,'Capsize', 10);
%             er1.LineStyle = 'none';
%             title('PW1000 - Peak');
%             subplot(212)
%             br2 = bar(stims,timeOns); hold on;
%             er2 = errorbar(stims, timeOns,timeOnsVar,'k', 'LineWidth', 2,'Capsize', 10);
%             er2.LineStyle = 'none';
%             title('PW1000 - Onset');
%             
%     
% %     end
% end