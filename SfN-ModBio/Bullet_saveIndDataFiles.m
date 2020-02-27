%Create Single Files form those made by the ModBio_Bullet script.

savepath = 'D:\DATA\SfN2018_ModBio\NHP-Bullet\';
saveName = 'BulletDrumTrial-';

disp('Loading Bullet Meta File');
load('BulletDrumTrials-dataFile.mat')

disp('Beginning Single File Parsing');

for i= info.fstart:info.fstop
    datafile = dataFiles{i};
    dFileName = string(erase(dataFiles{i}.filename,'.ns5'));
    dFileNum = string(erase(dFileName, 'datafile'));
    saveFileName = join([savepath saveName dFileNum '.mat'],'');
    disp(['Saving ' dFileName]);
    save(saveFileName, 'datafile', '-v7.3');
    pause (0.1);
    clear datafile
end