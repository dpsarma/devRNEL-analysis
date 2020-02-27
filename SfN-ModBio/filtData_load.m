[filenames, pathname] = uigetfile('*mat','Pick files','MultiSelect', 'on');


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


for fNum = 1:size(filenames,2)
    
    disp('Loading file...')
    load(string(fullfile(pathname,filenames(fNum))));
    
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