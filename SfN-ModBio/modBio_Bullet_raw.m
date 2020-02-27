%% ModBio Bullet Raw Look


disp('Loading Bullet Meta File');
load('BulletDrumTrials-dataFile.mat','info')
% [filenames, pathname] = uigetfile('*mat','Select all Files','MultiSelect', 'on');
% save('Bullet_DrumFileNames')

load('Bullet_DrumFileNames')

shankA = [1:8];
shankB = [9:16];
shankC = [17:24];
shankD = [25:32];
shanks = [shankA; shankB; shankC; shankD];
probeNames = ['Probe 1'; 'Probe 2';];
shankNames = ['Shank A'; 'Shank B'; 'Shank C'; 'Shank D'];



disp('Beginning Single File Parsing');
for fNum = 1:size(filenames,2)
    disp('Loading file...')
    load(string(fullfile(pathname,filenames(fNum))));
    disp('Plotting...')
    for probeNum = 1:info.numProbes
        %         handles.stackedRaw(fNum,probeNum) = figure;
        figure;
        maximize;
        suptitle(char(join([string(erase(filenames(fNum),'.mat')) ' - Probe ' num2str(probeNum) ' - raw'])));
        
        for shankNum = 1:size(shanks,1)
            subplot(2,2,shankNum)
            stackedplot(datafile.timeVec, datafile.shankData{probeNum}(shanks(shankNum,:),:),'PlotSpans','equal', 'ShowBaseline',false, 'DataTickMode', 'range');
            title(shankNames(shankNum,:));
        end
        
        suplabel('time(sec)','x');
        suplabel('channel(uV)','y');
    end
    pause(0.1)
    disp('Saving 2 Images...');
    saveallopenfigs_png(char(join(['D:\Figures\SFN modBio\Figures\raw_stacked' string(erase(filenames(fNum),'.mat')) '-Raw-'])));
    clear datafile;
end


disp('Beginning Single File Parsing for PSD');
for fNum = 1:size(filenames,2)
    disp('Loading file...')
    load(string(fullfile(pathname,filenames(fNum))));
    procData{fNum}.filename = datafile.filename;
    procData{fNum}.trialInfo = datafile.trialInfo;
    procData{fNum}.probe = datafile.probe;
    
    disp('Plotting...')
    for probeNum = 1:info.numProbes
        %         handles.stackedRaw(fNum,probeNum) = figure;
        figure;
        maximize;
        suptitle(char(join([string(erase(filenames(fNum),'.mat')) ' - Probe ' num2str(probeNum) ' - PSD'])));
        
        for shankNum = 1:size(shanks,1)
            subplot(2,2,shankNum)
            pwelch(datafile.shankData{probeNum}(shanks(shankNum,:),:)',[],[],[],30000,'power');
            [procData{fNum}.probe(probeNum).shank(shankNum).pxx, procData{fNum}.probe(probeNum).shank(shankNum).f] = pwelch(datafile.shankData{probeNum}(shanks(shankNum,:),:)',[],[],[],30000,'power');
            set(gca, 'XScale', 'log');
            legend('1','2','3','4','5','6','7','8');
            title(shankNames(shankNum,:));
        end
%         maximize;
        
%         suplabel('time(sec)','x');
%         suplabel('channel(uV)','y');
    end
    pause(0.1)
    disp('Saving 2 Images...');
    saveallopenfigs(char(join(['D:\Figures\SFN modBio\Figures\psd' string(erase(filenames(fNum),'.mat')) '-Raw-'])));
    clear datafile;
end


%% Plot just Pxx by Shank
fcount = 1+(info.fstop-info.fstart)

for probeNum = 1:info.numProbes
    
    disp(['Showing Probe Number - ' num2str(probeNum)]);
    
    for shankNum = 1:size(shanks,1)
        disp(['...' shankNames(shankNum,:) '...']);
        
        figure;
        hold on;
        for j=1:fcount
            subplot(round(fcount/2),2,j); %figIdx
            plot(procData{j}.probe(probeNum).shank(shankNum).f, 10*log(procData{j}.probe(probeNum).shank(shankNum).pxx))
            set(gca, 'XScale', 'log');
            xlabel('Frequency (Hz)');
            ylabel ('Power (dB)');
            ylim([-150 150])
            title([procData{j}.filename]); %trialInfo
            disp([procData{j}.filename]);
        end
        maximize;
        suptitle({['Probe ' num2str(probeNum) ': ' shankNames(shankNum,:) ' - PSD for subset "DRUM" dataFiles']});
        
    end
end