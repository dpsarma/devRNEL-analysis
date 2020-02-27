% read_Intan_RHD2000_file


path = 'D:\DATA\AnkleEMG\BGT10\';
[Filenames, Pathname] = uigetfile('D:\DATA\AnkleEMG\BGT10\*.rhd','Pick files','MultiSelect', 'on');

for i = length(Filenames)
    read_IntanRHD_v2(path,Filenames{i});
    schweave(:,:,i) = amplifier_data(1:32,:);
    axion(:,:,i) = amplifier_data(33:64,:);
    forceD(:,:,i) = board_adc_data;
    
end
meanSleeve = mean(schweave,3);
meanHD = mean(axion,3);
meanforce = mean(forceD,3);


% schweave = amplifier_data(1:32,:);
% axion = amplifier_data(33:64,:);
% forceD = board_adc_channels;

% plot(t_amplifier,schweave);
% figure; plot(t_amplifier,axion);

% Filtering Settings:
fs = 30e3;
lowCut = 10; %lowest frequency to pass
highCut = 7500; %highest frequency to pass
Norder = 4;
Wp = [lowCut, highCut]/(.5*fs);
[b,a]=butter(Norder, Wp);
Fsmooth = 10; %Lowpass cutoff frequency for smoothing the rectified signal 10 vs 50
NorderS =  4;
dt = 1/fs;
[bS,aS] = butter(NorderS,Fsmooth/(0.5*fs));


filtSleeve = filtfilt(b,a,meanSleeve')';
filtHD = filtfilt(b,a,meanHD')';
filtSleeve = abs(filtSleeve); filtHD = abs(filtHD);
figure; plot(t_amplifier,filtSleeve);
figure; plot(t_amplifier,filtHD);


movRMS = dsp.MovingRMS(3000);
rmsHD = movRMS(filtHD);
rmsSleeve = movRMS(filtSleeve);

figure; plot(t_amplifier,rmsHD);
figure; plot(t_amplifier,rmsSleeve);

startI = find(t_amplifier==0);
endI = find(t_amplifier==5);
 
figure; plot(t_amplifier(startI:end),rmsHD(:,startI:end));


smSleeve = filtfilt(bS,aS,filtSleeve')';
smHD = filtfilt(bS,aS,filtHD')';

figure; subplot(212); plot(t_amplifier,smHD);title('Plantarflexion: Biocircuits/Axion Array (Posterior)'); xlabel('Time(sec)'); ylabel('Smoothed/Rect EMG (uV)'); xlim([-1,6]);
vline([2 4],{'r', 'k'}, {'Cue End', 'Relax'});
subplot(211); plot(t_amplifier,smSleeve);title('Plantarflexion: miniSchweave Array (Anterior)'); xlabel('Time(sec)'); ylabel('Smoothed/Rect EMG (uV)'); xlim([-1,6]);
vline([2 4],{'r', 'k'}, {'Cue End', 'Relax'});

% % % create the video writer with 1 fps
% % writerObj = VideoWriter([path,'heatmap.avi']);
% % writerObj.FrameRate = 120;
% % % set the seconds per image
% % 
% % % open the video writer
% % open(writerObj);
% % % write the frames to the video
% % for i = startI:fs/writerObj.FrameRate:endI
% %     figure(1)
% %     subplot(211)
% %     imagesc(reshape(smSleeve(:,i),[4,8])); yticks(0:4); caxis([0,mean(mean(smSleeve(:,startI:endI)))+8*std(std(smSleeve(:,startI:endI)))]); colorbar;
% %     title('Schweave 32 - Anterior');
% %     subplot(212)
% %     imagesc(reshape(smHD(1:30,i),[5,6])); yticks(0:5); caxis([0,mean(mean(smHD(:,startI:endI)))+8*std(std(smHD(:,startI:endI)))]); colorbar; %mean(mean(smHD(:,startI:endI)))+2*std(std(smHD(:,startI:endI)))
% %     title('Axion 30 - Posterior');
% %     suptitle(['Plantarflexion: ' num2str(t_amplifier(i)) 'sec'])
% %     F(i) = getframe(gcf);
% %     
% %     % convert the image to a frame
% %     frame = F(i) ;    
% %     writeVideo(writerObj, frame);
% % end
% % % close the writer object
% % close(writerObj);





