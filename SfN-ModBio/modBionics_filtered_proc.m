savepath = 'D:\DATA\SfN2018_ModBio\NHP-Bullet\'; saveName = 'BulletDrumTrial-';

shankA = [1:8]; shankB = [9:16]; shankC = [17:24]; shankD = [25:32];
shanks = [shankA; shankB; shankC; shankD];
probeNames = ['Probe 1'; 'Probe 2';];
shankNames = ['Shank A'; 'Shank B'; 'Shank C'; 'Shank D'];

% Filtering Settings:
fs = 30e3;
lowCut = 100; %lowest frequency to pass 300
highCut = 10e3; %highest frequency to pass 3e3?
Norder = 4;
Wp = [lowCut, highCut]/(.5*fs);
[b,a]=butter(Norder, Wp);

disp('Loading Bullet Meta File');
load('BulletDrumTrials-dataFile.mat','info')
% [filenames, pathname] = uigetfile('*mat','Select all Files','MultiSelect', 'on');
% save('Bullet_DrumFileNames')

load('Bullet_DrumFileNames')

% load('D:\DATA\SfN2018_ModBio\NHP-Bullet\BulletDrumTrial-0073.mat')


disp('Beginning Single File Parsing');
for fNum = 1:size(filenames,2)
    disp('Loading file...')
    load(string(fullfile(pathname,filenames(fNum))));
    disp('Plotting...')
    
    filtData.filename = datafile.filename;
    filtData.trialInfo = datafile.trialInfo;
    filtData.timeVec = datafile.timeVec;

    %% Filter Data
    for probeNum = 1:size(datafile.probe,2)
        filtData.shankData{probeNum} = filtfilt(b,a, datafile.shankData{probeNum}')';
    %     datafile.shankData{probeNum} = xfilt';
    end

    clear datafile
    

    dFileName = string(erase(filtData.filename,'.ns5'));
    dFileNum = string(erase(dFileName, 'datafile'));
    saveFileName = join([savepath saveName dFileNum '-filtered.mat'],'');
    save(saveFileName, 'filtData', '-v7.3');

    %% Plot Single Times Courses
    for probeNum = 1:2
            for shankNum = 1:size(shanks,1)
                hSRaw = figure; maximize;
                plot(filtData.timeVec, filtData.shankData{probeNum}(shanks(shankNum,:),:));
                hold on
                title(char(join([filtData.filename ' - Probe ' num2str(probeNum) ' -  ' shankNames(shankNum,:)])));
                legend({'1'; '2'; '3'; '4'; '5'; '6'; '7'; '8'},'Location','northwest');
                xlabel('time(sec)');
                ylabel('channel(uV)');

                pause(0.1)
                disp('Saving 1 Image...');
                saveas(hSRaw,['D:\Figures\SFN modBio\Figures\\Filtered\timeseries-single\' char(join([erase(filtData.filename,'.ns5') ' - Probe ' num2str(probeNum) ' -  ' shankNames(shankNum,:)])) '-filt.png']);
%                 pause(0.1)
%                 savefig(hSRaw,['D:\Figures\SFN modBio\Figures\\Filtered\timeseries-single\' char(join([erase(filtData.filename,'.ns5') ' - Probe ' num2str(probeNum) ' -  ' shankNames(shankNum,:)]))]);

                close all
            end
    end

    %% Plotting Covariance/RMS/Correlations
    
    disp('Plotting Matrixes...');
    for probeNum = 1:info.numProbes
        x = rms(filtData.shankData{probeNum}');
        sdev = std(filtData.shankData{probeNum}');
        
        hCR = figure; 
        maximize;
        for shankNum = 1:size(shanks,1)
            subplot(4,3,3*shankNum-2); 
            imagesc(cov(filtData.shankData{probeNum}([shanks(shankNum,:)],:)')); 
%             caxis([min(min(minCov)) max(max(maxCov))]); 
            colorbar;
            title([probeNames(probeNum,:) ' - ' shankNames(shankNum,:) ' - Covariance']);

            hBar(shankNum) = subplot(4,3,3*shankNum-1);
            bar(shanks(shankNum,:), x(shanks(shankNum,:)));
            title([probeNames(probeNum,:) ' - ' shankNames(shankNum,:) ' - Average Power']);
            ylabel('RMS (uV)');
            xlabel('Channels');
            
            subplot(4,3,3*shankNum); 
            imagesc(corr(filtData.shankData{probeNum}([shanks(shankNum,:)],:)')); caxis([-1 1]); colorbar;
            title([probeNames(probeNum,:) ' - ' shankNames(shankNum,:) ' - Correlations']);
        end
        suptitle(filenames(fNum));
        linkaxes([hBar],'y');

        saveas(hCR,['D:\Figures\SFN modBio\Figures\Filtered\CovRMSCorr\' char(erase(filenames(fNum),'.mat'))...
            '-' probeNames(probeNum,:) '-filtered-CovRMSCorr.png']);
        savefig(hCR,['D:\Figures\SFN modBio\Figures\Filtered\CovRMSCorr\' char(erase(filenames(fNum),'.mat'))...
            '-' probeNames(probeNum,:) '-filtered-CovRMSCorr']);
        close all;
        
    end
    
    %% Plotting PSDs...
    
    disp('Plotting PSDs...')
    for probeNum = 1:info.numProbes
        %         handles.stackedRaw(fNum,probeNum) = figure;
        figure;
        maximize;
        suptitle(char(join([string(erase(filenames(fNum),'.mat')) ' - Probe ' num2str(probeNum) ' - (filt) PSD'])));
        
        for shankNum = 1:size(shanks,1)
            subplot(2,2,shankNum)
            pwelch(filtData.shankData{probeNum}(shanks(shankNum,:),:)',[],[],[],30000,'power');
            set(gca, 'XScale', 'log');
            legend('1','2','3','4','5','6','7','8','Location', 'northwest');
            title(shankNames(shankNum,:));
        end
%         maximize;
        
%         suplabel('time(sec)','x');
%         suplabel('channel(uV)','y');
    end
    pause(0.1)
    disp('Saving 2 Images...');
    saveallopenfigs(char(join(['D:\Figures\SFN modBio\Figures\Filtered\PSD\' string(erase(filenames(fNum),'.mat')) '-Filt-'])));

    
    

    
end