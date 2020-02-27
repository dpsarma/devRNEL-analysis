disp('Beginning Single File Parsing for PSD');
for fNum = 1:size(filenames,2)
    disp('Loading file...')
    load(string(fullfile(pathname,filenames(fNum))));
    procData{fNum}.filename = datafile.filename;
    procData{fNum}.trialInfo = datafile.trialInfo;
    procData{fNum}.probe = datafile.probe;
    
    disp('Pwelching...')
    for probeNum = 1:info.numProbes
       
        for shankNum = 1:size(shanks,1)
           
            [procData{fNum}.probe(probeNum).shank(shankNum).pxx, procData{fNum}.probe(probeNum).shank(shankNum).f] = pwelch(datafile.shankData{probeNum}(shanks(shankNum,:),:)',[],[],[],30000,'power');
        end
    end
    disp(['Finishing last file:' datafile.filename])
    pause(0.1)
    clear datafile;
end


%% Plot ALL PSD's on one figure.
fcount = 1+(info.fstop-info.fstart)
probeNames = ['Probe 1 - '; 'Probe 2 - ';];
shankNames = ['Shank A'; 'Shank B'; 'Shank C'; 'Shank D'];

figure;
maximize;
suptitle('PSD for ALL Bullet-Drum Data');
figCount = 1;
for fNum = info.fstart:info.fstop
    disp(['Plotting PSDs for File: ' num2str(fNum)])
    for probeNum = 1:info.numProbes
        for shankNum = 1:size(info.shanks,1)
            disp(char(join([erase(dataFiles{fNum}.filename,'.ns5') ': ' probeNames(probeNum,:) shankNames(shankNum,:)])));
            subplot(fcount,(info.numProbes*size(info.shanks,1)), figCount);
            pwelch(dataFiles{fNum}.shankData{probeNum}(info.shanks(shankNum,:),:)',[],[],[],30000,'power');
            set(gca, 'XScale', 'log');
            set(gca,'fontsize',7)
%             legend('1','2','3','4','5','6','7','8');
            title(char(join([erase(dataFiles{fNum}.filename,'.ns5') ': ' probeNames(probeNum,:) shankNames(shankNum,:)])),'FontSize',5);
            figCount = figCount+1;
        end
    end
end
               