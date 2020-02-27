

disp('Loading Bullet Meta File');
load('BulletDrumTrials-dataFile.mat','info')
load('Bullet_DrumFileNames')

shankA = [1:8];
shankB = [9:16];
shankC = [17:24];
shankD = [25:32];
shanks = [shankA; shankB; shankC; shankD];
probeNames = ['Probe 1'; 'Probe 2';];
shankNames = ['Shank A'; 'Shank B'; 'Shank C'; 'Shank D'];

disp('Loading Cov matrix File')
load('D:\DATA\SfN2018_ModBio\NHP-Bullet\BulletProcessedData-corr_cov_pxx.mat')

disp('Beginning Single File Parsing');

for fNum = 1:size(filenames,2)
    
    for probeNum = 1:info.numProbes
        
        for shankNum = 1:size(shanks,1)
            xMin(shankNum) = min(min(procData{fNum}.probe(probeNum).shank(shankNum).cov))
            xMax(shankNum) = max(max(procData{fNum}.probe(probeNum).shank(shankNum).cov))
        end
        minCov(fNum, probeNum) = min(xMin)
        maxCov(fNum, probeNum) = max(xMax)
    end
end



for fNum = 1:size(filenames,2)
    disp('Loading file...')
    load(string(fullfile(pathname,filenames(fNum))));
    disp('Plotting...')
    for probeNum = 1:info.numProbes
        x = rms(datafile.shankData{probeNum}')
        sdev = std(datafile.shankData{probeNum}')
        
        hCR = figure; 
        maximize;
        for shankNum = 1:size(shanks,1)
            subplot(4,3,3*shankNum-2); 
            imagesc(procData{fNum}.probe(probeNum).shank(shankNum).cov); 
            caxis([min(min(minCov)) max(max(maxCov))]); colorbar;
            title([probeNames(probeNum,:) ' - ' shankNames(shankNum,:) ' - Covariance']);

            hBar(shankNum) = subplot(4,3,3*shankNum-1);
            bar(shanks(shankNum,:), x(shanks(shankNum,:)));
            title([probeNames(probeNum,:) ' - ' shankNames(shankNum,:) ' - Average Power']);
            ylabel('RMS (uV)');
            xlabel('Channels');
            
            subplot(4,3,3*shankNum); 
            imagesc(procData{fNum}.probe(probeNum).shank(shankNum).corr); caxis([-1 1]); colorbar;
            title([probeNames(probeNum,:) ' - ' shankNames(shankNum,:) ' - Correlations']);
        end
        suptitle(filenames(fNum));
        linkaxes([hBar],'y');
%         savefig(hCR,['D:\Figures\SFN modBio\Figures\' char(erase(filenames(fNum),'.mat'))...
%             '-' probeNames(probeNum,:) '-CovRMSCorr']);
        saveas(hCR,['D:\Figures\SFN modBio\Figures\' char(erase(filenames(fNum),'.mat'))...
            '-' probeNames(probeNum,:) '-CovRMSCorr2.png']);
        close all;
        
    end
end