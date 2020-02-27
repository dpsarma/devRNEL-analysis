% Replace Standing with Standing
addpath('C:\Users\dsarma\Documents\GitHub\mattools\linspecer\');
% addpath('C:\Users\dsarma\Documents\GitHub\mattools\gridLegend_v1.4\');
C2 = linspecer(13);

Muscle_ID = {'Right VM', 'Right RF', 'Right VL', 'Right BF', 'Right ST', ...
    'Right TA', 'Right SO', 'Right LG', 'Left VM', 'Left RF', 'Left VL',...
    'Left BF', 'Left ST', 'Left Ham', 'Left SO', 'Left LG'};
path_Save = 'R:\users\dsarma\';
path_UH3 = 'R:\data_raw\human\uh3_stim\';
subjectName = 'LSP02b';
reportPath = ['R:\users\dsarma\' subjectName '\'];

fs =  30e3;
dtPre = 50e-3; % msec before stim onset
dtPost= 5e-3; %msec after stim onset
nSampsPre = floor((dtPre-5e-3)*fs);
nSampsPost = floor((dtPost+dtPre)*fs);
% minPeak = 5;
gWind = 50;

for m = 1:length(Muscle_ID)
%     hF4(m) = figure; maximize;
    cCount = 1;
    
    disp (['Processing Muscle: ' Muscle_ID(m)]);
    
    for hIdx = 1:70
        if ~contains(hash(hIdx).muscleNames, Muscle_ID(m))
            continue;
        elseif ~isfield(hash(hIdx).standing,'stimAmp') || hIdx == 9;
            disp(['No Standing File: ' num2str(hIdx)]);
            continue;
        else
            peakInfo(hIdx).muscle(m).Name = Muscle_ID(m);
            stimAmps4 = [hash(hIdx).standing.stimAmp];
            [sortedStim4,I] = sort(stimAmps4);
                       
            for ee = 1:length(hash(hIdx).Indices.standing) 
                disp([Muscle_ID(m) '- Set ' hIdx ' - Trial ' ee]);
%             X(ee,:) = hash(hIdx).standing(I(ee)).epochTimeVec;
%             Y(ee,:) = abs(hash(hIdx).standing(I(ee)).meanTrace{m});
            
%             minPeak = round(mean(abs(hash(hIdx).standing(I(ee)).meanTrace{m}(nSampsPre:end)))+std(abs(hash(hIdx).standing(I(ee)).meanTrace{m}(nSampsPre:end))));
                minPeak = (2*std(smoothdata(abs(hash(hIdx).standing(I(ee)).meanTrace{m}(1:nSampsPre))','gaussian',gWind))+mean(smoothdata(abs(hash(hIdx).standing(I(ee)).meanTrace{m}(1:nSampsPre))','gaussian',gWind)));
% %                 minPeak = var(smoothdata(abs(hash(hIdx).standing(I(ee)).meanTrace{m}(1:nSampsPre))','gaussian',gWind));
                [pks,locs,widths,proms] = findpeaks(smoothdata(abs(hash(hIdx).standing(I(ee)).meanTrace{m}(nSampsPost:end))','gaussian',gWind)',...
                    hash(hIdx).standing(I(ee)).epochTimeVec(nSampsPost:end), ...
                    'Annotate','extents', 'WidthReference','halfheight','MinPeakHeight',minPeak); 
                peakInfo(hIdx).muscle(m).standing(ee).pks = pks;
                peakInfo(hIdx).muscle(m).standing(ee).locs = locs+dtPost;
                peakInfo(hIdx).muscle(m).standing(ee).widths = widths;
                peakInfo(hIdx).muscle(m).standing(ee).proms = proms;
                peakInfo(hIdx).muscle(m).standing(ee).stimAmp = sortedStim4(ee);
                peakInfo(hIdx).muscle(m).standing(ee).day = hash(hIdx).day;
                peakInfo(hIdx).muscle(m).standing(ee).elecNum = hash(hIdx).elecNum;
                peakInfo(hIdx).muscle(m).standing(ee).mean = mean(smoothdata(abs(hash(hIdx).standing(I(ee)).meanTrace{m}(nSampsPost:nSampsPost+nSampsPost))','gaussian',gWind)); %Mean in RMS Window
                peakInfo(hIdx).muscle(m).standing(ee).var = var(smoothdata(abs(hash(hIdx).standing(I(ee)).meanTrace{m}(nSampsPost:nSampsPost+nSampsPost))','gaussian',gWind));
                peakInfo(hIdx).muscle(m).standing(ee).std = std(smoothdata(abs(hash(hIdx).standing(I(ee)).meanTrace{m}(nSampsPost:nSampsPost+nSampsPost))','gaussian',gWind));

% %             %% Plot Window
% %             figH = figure;maximize;
% %             h(1) = subplot(211);
% %             plot(hash(hIdx).standing(I(ee)).epochTimeVec,hash(hIdx).standing(I(ee)).meanTrace{m})
% %             hold on;
% %             vline([(0-5e-3) dtPre dtPost],{'r:','r:','r:'},{'PreStim', 'RMSWindow','PostStim'});
% %             plot(hash(hIdx).standing(I(ee)).epochTimeVec(nSampsPost:end),hash(hIdx).standing(I(ee)).meanTrace{m}(nSampsPost:end),'r')
% %             plot(hash(hIdx).standing(I(ee)).epochTimeVec(1:nSampsPre),hash(hIdx).standing(I(ee)).meanTrace{m}(1:nSampsPre),'g')
% %             hline([(2*std(hash(hIdx).standing(I(ee)).meanTrace{m}(1:nSampsPre))+mean(hash(hIdx).standing(I(ee)).meanTrace{m}(1:nSampsPre))) mean(hash(hIdx).standing(I(ee)).meanTrace{m}(1:nSampsPre))],{'k:', 'k:'},{'2STD', 'mean'})
% %             title([Muscle_ID(m) ' wiggles']);
% %             ylabel('Voltage (mV)')
% %             legend('Raw Mean Trace','Response Period', 'PreStim Period');
% % 
% %             h(2) = subplot(212);
% %             plot(hash(hIdx).standing(I(ee)).epochTimeVec,smoothdata(abs(hash(hIdx).standing(I(ee)).meanTrace{m})','gaussian',gWind))
% %             hold on; vline([(0-5e-3) dtPre dtPost],{'r:','r:','r:'},{'PreStim', 'RMSWindow','PostStim'});
% %             plot(hash(hIdx).standing(I(ee)).epochTimeVec(1:nSampsPre),smoothdata(abs(hash(hIdx).standing(I(ee)).meanTrace{m}(1:nSampsPre))','gaussian',gWind),'g')
% %             plot(hash(hIdx).standing(I(ee)).epochTimeVec(nSampsPost:end),smoothdata(abs(hash(hIdx).standing(I(ee)).meanTrace{m}(nSampsPost:end))','gaussian',gWind),'r-')
% %             hline([(2*std(smoothdata(abs(hash(hIdx).standing(I(ee)).meanTrace{m}(1:nSampsPre))','gaussian',gWind))+mean(smoothdata(abs(hash(hIdx).standing(I(ee)).meanTrace{m}(1:nSampsPre))','gaussian',gWind)))...
% %                 mean(smoothdata(abs(hash(hIdx).standing(I(ee)).meanTrace{m}(1:nSampsPre))','gaussian',gWind)) var(smoothdata(abs(hash(hIdx).standing(I(ee)).meanTrace{m}(1:nSampsPre))','gaussian',gWind))],{'k:', 'k:', 'k:'}, {'2STD', 'mean', 'var'})
% %             title([Muscle_ID(m) ' rectified+smoothed wiggles']);
% %             ylabel('Rectified Voltage (mV)')
% %             xlabel('Time (sec)');
% %             linkaxes([h],'x');
% %             h(1).XLim = [-.05 .1];
% %             legend('Rect/Smooth Mean Trace','Response Period', 'PreStim Period');
% %             disp('Saving Peaks to Set Folder')
% %             saveas(figH,[reportPath Muscle_ID{m} '_signalWindow_Day' num2str(hash(hIdx).day) '_Elec' num2str(hash(hIdx).elecNum) '_standing_' num2str(hIdx) '-' num2str(sortedStim4(ee)) 'mA.png']);
% % 
% %             %% Plot Peaks
% %             hSP = figure; %maximize;            
% %             findpeaks(smoothdata(abs(hash(hIdx).standing(I(ee)).meanTrace{m}(nSampsPost:end))','gaussian',gWind)',...
% %                 hash(hIdx).standing(I(ee)).epochTimeVec(nSampsPost:end), ...
% %                 'Annotate','extents', 'WidthReference','halfheight','MinPeakHeight',minPeak); %'SortStr','descend'
% %             text(locs+.02,pks,num2str((1:numel(pks))'))
% %             title(['Signal Peak Widths: ' Muscle_ID(m) ' Day ' num2str(hash(hIdx).day) ' Elec ' num2str(hash(hIdx).elecNum) ' - ' num2str(sortedStim4(ee)) 'mA']);
% %             ylabel('Rectified Voltage (mV)')
% %             xlabel('Time (sec)');
% %             disp('Saving Peaks to Set Folder')
% %             saveas(hSP,[reportPath Muscle_ID{m} '_signalPeaks_Day' num2str(hash(hIdx).day) '_Elec' num2str(hash(hIdx).elecNum) '_standing_' num2str(hIdx) '-' num2str(sortedStim4(ee)) 'mA.png']);
% %             close all;
% %             
            end
            
            
            
        end
    end
end  
%             clearvars -except hash