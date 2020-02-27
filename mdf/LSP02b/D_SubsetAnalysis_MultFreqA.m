%% UH3 - LSP02 anaylsis (4a) - MDF Analysis
%----------------------------------------
% tag that is missing from all EMG trials is "Sitting" vs. "Standing"

%% Load for subset of parameters
%omdfc = mdf.init('auto');
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
%     if (freqstring == 2 || freqstring == 5)
%         disp('Moving On');
%         continue;
%     end
    
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
                                if freqstring == 1
                                     E2_F1(m,tnum).Stims = hash{m,e,f,tnum}.stimAmp;
                                     E2_F1(m,tnum).P2P = hash{m,e,f,tnum}.meanP2P;
                                     E2_F1(m,tnum).P2Perr = hash{m,e,f,tnum}.errP2P;
                                     E2_F1(m,tnum).meantrace = mean(hash{m,e,f,tnum}.epochsrect);
                                     E2_F1(m,tnum).epochtime = hash{m,e,f,tnum}.epochtime;


                                     for i = 1:size(hash{m,e,f,tnum}.epochsrect,1)
                                         [pks,locs,widths,proms] = findpeaks(smoothdata(hash{m,e,f,tnum}.epochsrect(nSampsPost:nSampsPost+respSamps)','gaussian',gWind)',...
                                        E2_F1(m,tnum).epochtime(nSampsPost:nSampsPost+respSamps), ...
                                        'Annotate','extents', 'WidthReference','halfheight'); 
                                        [maxpeak(i), ipeak] = max(pks);
                                        timeonset(i) = locs(ipeak);
                                     end
                                        E2_F1(m,tnum).pkVar = std(maxpeak)/sqrt(length(maxpeak));
                                        E2_F1(m,tnum).timeOnsVar = std(timeonset)/sqrt(length(timeonset));

                                     [pks,locs,widths,proms] = findpeaks(smoothdata(E2_F1(m,tnum).meantrace(nSampsPost:nSampsPost+respSamps)','gaussian',gWind)',...
                                        E2_F1(m,tnum).epochtime(nSampsPost:nSampsPost+respSamps), ...
                                        'Annotate','extents', 'WidthReference','halfheight'); 
                                     [E2_F1(m,tnum).maxMeanpeak, ipeak] = max(pks);
                                     E2_F1(m,tnum).timeMeanpeak = locs(ipeak);
                                     clear maxpeak timeonset


                                elseif freqstring == 2
                                     E2_F2(m,tnum).Stims = hash{m,e,f,tnum}.stimAmp;
                                     E2_F2(m,tnum).P2P = hash{m,e,f,tnum}.meanP2P;
                                     E2_F2(m,tnum).P2Perr = hash{m,e,f,tnum}.errP2P;
                                     E2_F2(m,tnum).meantrace = mean(hash{m,e,f,tnum}.epochsrect);
                                     E2_F2(m,tnum).epochtime = hash{m,e,f,tnum}.epochtime;


                                     for i = 1:size(hash{m,e,f,tnum}.epochsrect,1)
                                         [pks,locs,widths,proms] = findpeaks(smoothdata(hash{m,e,f,tnum}.epochsrect(nSampsPost:nSampsPost+respSamps)','gaussian',gWind)',...
                                        E2_F2(m,tnum).epochtime(nSampsPost:nSampsPost+respSamps), ...
                                        'Annotate','extents', 'WidthReference','halfheight'); 
                                        [maxpeak(i), ipeak] = max(pks);
                                        timeonset(i) = locs(ipeak);
                                     end
                                        E2_F2(m,tnum).pkVar = std(maxpeak)/sqrt(length(maxpeak));
                                        E2_F2(m,tnum).timeOnsVar = std(timeonset)/sqrt(length(timeonset));

                                     [pks,locs,widths,proms] = findpeaks(smoothdata(E2_F2(m,tnum).meantrace(nSampsPost:nSampsPost+respSamps)','gaussian',gWind)',...
                                        E2_F2(m,tnum).epochtime(nSampsPost:nSampsPost+respSamps), ...
                                        'Annotate','extents', 'WidthReference','halfheight'); 
                                     [E2_F2(m,tnum).maxMeanpeak, ipeak] = max(pks);
                                     E2_F2(m,tnum).timeMeanpeak = locs(ipeak);
                                     clear maxpeak timeonset

                                elseif freqstring == 5
                                     E2_F5(m,tnum).Stims = hash{m,e,f,tnum}.stimAmp;
                                     E2_F5(m,tnum).P2P = hash{m,e,f,tnum}.meanP2P;
                                     E2_F5(m,tnum).P2Perr = hash{m,e,f,tnum}.errP2P;
                                     E2_F5(m,tnum).meantrace = mean(hash{m,e,f,tnum}.epochsrect);
                                     E2_F5(m,tnum).epochtime = hash{m,e,f,tnum}.epochtime;


                                     for i = 1:size(hash{m,e,f,tnum}.epochsrect,1)
                                         [pks,locs,widths,proms] = findpeaks(smoothdata(hash{m,e,f,tnum}.epochsrect(nSampsPost:nSampsPost+respSamps)','gaussian',gWind)',...
                                        E2_F5(m,tnum).epochtime(nSampsPost:nSampsPost+respSamps), ...
                                        'Annotate','extents', 'WidthReference','halfheight'); 
                                        if isempty(pks)
                                            disp(['No Response: ' musclabels(m)])
                                         continue;
                                        end
                                        [maxpeak(i), ipeak] = max(pks);
                                        timeonset(i) = locs(ipeak);
                                     end
                                        if isempty(pks)
                                            disp(['No Response: ' musclabels(m)])
                                         continue;
                                        end
                                        E2_F5(m,tnum).pkVar = std(maxpeak)/sqrt(length(maxpeak));
                                        E2_F5(m,tnum).timeOnsVar = std(timeonset)/sqrt(length(timeonset));

                                     [pks,locs,widths,proms] = findpeaks(smoothdata(E2_F5(m,tnum).meantrace(nSampsPost:nSampsPost+respSamps)','gaussian',gWind)',...
                                        E2_F5(m,tnum).epochtime(nSampsPost:nSampsPost+respSamps), ...
                                        'Annotate','extents', 'WidthReference','halfheight'); 
                                     [E2_F5(m,tnum).maxMeanpeak, ipeak] = max(pks);
                                     E2_F5(m,tnum).timeMeanpeak = locs(ipeak);
                                     
                                     clear maxpeak timeonset
                                        
                               
                                else
                                    continue;
                                end
                                
                            case 4
                                if freqstring == 1
                                     E4_F1(m,tnum).Stims = hash{m,e,f,tnum}.stimAmp;
                                     E4_F1(m,tnum).P2P = hash{m,e,f,tnum}.meanP2P;
                                     E4_F1(m,tnum).P2Perr = hash{m,e,f,tnum}.errP2P;
                                     E4_F1(m,tnum).meantrace = mean(hash{m,e,f,tnum}.epochsrect);
                                     E4_F1(m,tnum).epochtime = hash{m,e,f,tnum}.epochtime;


                                     for i = 1:size(hash{m,e,f,tnum}.epochsrect,1)
                                         [pks,locs,widths,proms] = findpeaks(smoothdata(hash{m,e,f,tnum}.epochsrect(nSampsPost:nSampsPost+respSamps)','gaussian',gWind)',...
                                        E4_F1(m,tnum).epochtime(nSampsPost:nSampsPost+respSamps), ...
                                        'Annotate','extents', 'WidthReference','halfheight'); 
                                        [maxpeak(i), ipeak] = max(pks);
                                        timeonset(i) = locs(ipeak);
                                     end
                                        E4_F1(m,tnum).pkVar = std(maxpeak)/sqrt(length(maxpeak));
                                        E4_F1(m,tnum).timeOnsVar = std(timeonset)/sqrt(length(timeonset));

                                     [pks,locs,widths,proms] = findpeaks(smoothdata(E4_F1(m,tnum).meantrace(nSampsPost:nSampsPost+respSamps)','gaussian',gWind)',...
                                        E4_F1(m,tnum).epochtime(nSampsPost:nSampsPost+respSamps), ...
                                        'Annotate','extents', 'WidthReference','halfheight'); 
                                     [E4_F1(m,tnum).maxMeanpeak, ipeak] = max(pks);
                                     E4_F1(m,tnum).timeMeanpeak = locs(ipeak);
                                     clear maxpeak timeonset


                                elseif freqstring == 2
                                     E4_F2(m,tnum).Stims = hash{m,e,f,tnum}.stimAmp;
                                     E4_F2(m,tnum).P2P = hash{m,e,f,tnum}.meanP2P;
                                     E4_F2(m,tnum).P2Perr = hash{m,e,f,tnum}.errP2P;
                                     E4_F2(m,tnum).meantrace = mean(hash{m,e,f,tnum}.epochsrect);
                                     E4_F2(m,tnum).epochtime = hash{m,e,f,tnum}.epochtime;


                                     for i = 1:size(hash{m,e,f,tnum}.epochsrect,1)
                                         [pks,locs,widths,proms] = findpeaks(smoothdata(hash{m,e,f,tnum}.epochsrect(nSampsPost:nSampsPost+respSamps)','gaussian',gWind)',...
                                        E4_F2(m,tnum).epochtime(nSampsPost:nSampsPost+respSamps), ...
                                        'Annotate','extents', 'WidthReference','halfheight'); 
                                        [maxpeak(i), ipeak] = max(pks);
                                        timeonset(i) = locs(ipeak);
                                     end
                                        E4_F2(m,tnum).pkVar = std(maxpeak)/sqrt(length(maxpeak));
                                        E4_F2(m,tnum).timeOnsVar = std(timeonset)/sqrt(length(timeonset));

                                     [pks,locs,widths,proms] = findpeaks(smoothdata(E4_F2(m,tnum).meantrace(nSampsPost:nSampsPost+respSamps)','gaussian',gWind)',...
                                        E4_F2(m,tnum).epochtime(nSampsPost:nSampsPost+respSamps), ...
                                        'Annotate','extents', 'WidthReference','halfheight'); 
                                     [E4_F2(m,tnum).maxMeanpeak, ipeak] = max(pks);
                                     E4_F2(m,tnum).timeMeanpeak = locs(ipeak);
                                     clear maxpeak timeonset

                                elseif freqstring == 5
                                     E4_F5(m,tnum).Stims = hash{m,e,f,tnum}.stimAmp;
                                     E4_F5(m,tnum).P2P = hash{m,e,f,tnum}.meanP2P;
                                     E4_F5(m,tnum).P2Perr = hash{m,e,f,tnum}.errP2P;
                                     E4_F5(m,tnum).meantrace = mean(hash{m,e,f,tnum}.epochsrect);
                                     E4_F5(m,tnum).epochtime = hash{m,e,f,tnum}.epochtime;


                                     for i = 1:size(hash{m,e,f,tnum}.epochsrect,1)
                                         [pks,locs,widths,proms] = findpeaks(smoothdata(hash{m,e,f,tnum}.epochsrect(nSampsPost:nSampsPost+respSamps)','gaussian',gWind)',...
                                        E4_F5(m,tnum).epochtime(nSampsPost:nSampsPost+respSamps), ...
                                        'Annotate','extents', 'WidthReference','halfheight'); 
                                        [maxpeak(i), ipeak] = max(pks);
                                        timeonset(i) = locs(ipeak);
                                     end
                                        E4_F5(m,tnum).pkVar = std(maxpeak)/sqrt(length(maxpeak));
                                        E4_F5(m,tnum).timeOnsVar = std(timeonset)/sqrt(length(timeonset));

                                     [pks,locs,widths,proms] = findpeaks(smoothdata(E4_F5(m,tnum).meantrace(nSampsPost:nSampsPost+respSamps)','gaussian',gWind)',...
                                        E4_F5(m,tnum).epochtime(nSampsPost:nSampsPost+respSamps), ...
                                        'Annotate','extents', 'WidthReference','halfheight'); 
                                     [E4_F5(m,tnum).maxMeanpeak, ipeak] = max(pks);
                                     E4_F5(m,tnum).timeMeanpeak = locs(ipeak);
                                     clear maxpeak timeonset

                               else
                                    continue;
                               end
                        end
        
        
                    end
                end
            end
    end
end
