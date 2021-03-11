% Filtering Settings:
fs = 30e3;
nyq = 0.5*fs;

n = 115; m =12;


%Bandpass
for m = 16
    figure; maximize;
    tiledlayout(6,3);
for bp = 1:18
switch bp
    case 1
        lowCut = 30;%30; %lowest frequency to pass
        highCut = 800;%800; %highest frequency to pass
        Norder = 1;
        notch = 1;
    case 2
        lowCut = 30;%30; %lowest frequency to pass
        highCut = 800;%800; %highest frequency to pass
        Norder = 2;
        notch = 1;

    case 3
        lowCut = 30;%30; %lowest frequency to pass
        highCut = 800;%800; %highest frequency to pass
        Norder = 4;
        notch = 1;
    case 4
        lowCut = 30;%30; %lowest frequency to pass
        highCut = 800;%800; %highest frequency to pass
        Norder = 1;
        notch = 0;        
    case 5        
        lowCut = 30;%30; %lowest frequency to pass
        highCut = 800;%800; %highest frequency to pass
        Norder = 2;
        notch = 0;
    case 6      
        lowCut = 30;%30; %lowest frequency to pass
        highCut = 800;%800; %highest frequency to pass
        Norder = 4;
        notch = 0;
    case 7
        lowCut = 75;%30; %lowest frequency to pass
        highCut = 7500;%800; %highest frequency to pass
        Norder = 1;
        notch = 1;
    case 8
        lowCut = 75;%30; %lowest frequency to pass
        highCut = 7500;%800; %highest frequency to pass
        Norder = 2;
        notch = 1;

    case 9
        lowCut = 75;%30; %lowest frequency to pass
        highCut = 7500;%800; %highest frequency to pass
        Norder = 4;
        notch = 1;
    case 10
        lowCut = 75;%30; %lowest frequency to pass
        highCut = 7500;%800; %highest frequency to pass
        Norder = 1;
        notch = 0;        
    case 11        
        lowCut = 75;%30; %lowest frequency to pass
        highCut = 7500;%800; %highest frequency to pass
        Norder = 2;
        notch = 0;
    case 12      
        lowCut = 75;%30; %lowest frequency to pass
        highCut = 7500;%800; %highest frequency to pass
        Norder = 4;
        notch = 0;
            case 13
        lowCut = 75;%30; %lowest frequency to pass
        highCut = 800;%800; %highest frequency to pass
        Norder = 1;
        notch = 1;
    case 14
        lowCut = 75;%30; %lowest frequency to pass
        highCut = 800;%800; %highest frequency to pass
        Norder = 2;
        notch = 1;

    case 15
        lowCut = 75;%30; %lowest frequency to pass
        highCut = 800;%800; %highest frequency to pass
        Norder = 4;
        notch = 1;
    case 16
        lowCut = 75;%30; %lowest frequency to pass
        highCut = 800;%800; %highest frequency to pass
        Norder = 1;
        notch = 0;        
    case 17        
        lowCut = 75;%30; %lowest frequency to pass
        highCut = 800;%800; %highest frequency to pass
        Norder = 2;
        notch = 0;
    case 18      
        lowCut = 75;%30; %lowest frequency to pass
        highCut = 800;%800; %highest frequency to pass
        Norder = 4;
        notch = 0;
end


    emg = trial(n).emg;
    timeVec = trial(n).timeVec ;
    Wp = [lowCut, highCut]/(nyq);
    [b,a]=butter(Norder, Wp);
    if notch == 1
        d = designfilt('bandstopiir','FilterOrder',2, ...
                   'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
                   'DesignMethod','butter','SampleRate',fs);
        emgF = filtfilt(d,emg')';
        emgF = filtfilt(b,a,emgF')';
        tString = '+Notch';
    else
        emgF = filtfilt(b,a,emg')';
        tString = ' ';
    end
    
    stimBegins = trial(n).stims;
    stimEnd = nSampsPre + trial(n).stimLength + artifactBuffer*fs;  
    %% Epoch data into Windows 
    disp(['Epoching' num2str(fnum)]);
    for i = 1:size(emg,1)
        asd = emgF(i,:);
        emgStim{i} = cell2mat(arrayfun(@(x) asd((stimBegins(x)-nSampsPre):(stimBegins(x)+(trial(fnum).stimLength+nSampsPost-1))), 1:length(stimBegins), 'UniformOutput',false)');
    end
    
    
    meanTrace = mean(emgStim{1,m});
% %     nexttile;
% %     plot(trial(n).epochTimeVec, abs(meanTrace),'LineWidth',2);
% %     hold on;
% %   plot(trial(n).epochTimeVec, abs(trial(n).meanTrace{1, m}),'LineWidth',2);
  
  
% %   axx(1) = nexttile;
% %       plot(trial(n).epochTimeVec, trial(n).meanTrace{1, m},'LineWidth',2);
% %         title('Original Mean Trace: 30:800 N2 (Notch/filtfilt)')
        
  axx(bp) = nexttile;
      plot(trial(n).epochTimeVec, meanTrace,'LineWidth',2);
      title([num2str(lowCut) ':' num2str(highCut) ', N' num2str(Norder) tString  ]);
%     hold on;
  
clear emg
clear emgF
clear emgStim
end
linkaxes([axx],'xy');
sgtitle(Muscle_ID(m));
end