artifactBuffer = 0.02;
rmsWindow = 0.02;

for fnum = 1:length(emgFilenames)
    stimEnd = nSampsPre + trial(fnum).stimLength + artifactBuffer*fs;    
    for m = 1:length(Muscle_ID)
        % Response
        trial(fnum).responseRMS(m,:) = mean(trial(fnum).epochRMS{m}(stimEnd:stimEnd+(rmsWindow*fs)));
        trial(fnum).responseP2P(m,:) = mean(trial(fnum).emgStim{m}(stimEnd:stimEnd+(rmsWindow*fs)));
    end
end

%% Plotting Recruitment Curves
response = [trial(:).responseRMS]';
% response = [trial(:).responseP2P]';

switch length(Muscle_ID)
    case 16
        pRow = 2; pCol = 3;
        mCount = length(Muscle_ID)/2;
        pIdx = 2;
    case 8
        pRow = 1; pCol = 3;
        mCount = length(Muscle_ID);
        pIdx = 1;
    otherwise
        disp('How Would you like to Plot the Recruitment Curves?');
        return;
end

hF3 = figure; %maximize;
[sortedStim,I] = sort(stimAmp);
for pI = 1:pIdx
    
   for m = 1+mCount*(pI-1):mCount+mCount*(pI-1)
       yyA = smooth(sortedStim/1000,response(I,m),'sgolay');
        plot(sortedStim/1000,yyA,'Color',C(m,:),'LineWidth',3);
%         scatter(stimAmp/1000, response(:,m)');
        hold on;
    end
        ylabel('(post-Stim) Mean RMS  (uV)')
        xlabel('Stimulation Amplitude (mA)')
        title(['By Amplitude']);
        
    if pIdx == 1
        legend([Muscle_ID(1:mCount)],'Location','eastoutside', 'Orientation','vertical');
    elseif pIdx == 2
        legend([Muscle_ID((mCount)+1:length(Muscle_ID))],'Location','eastoutside', 'Orientation','vertical');
    end
%     
end
tmptit = suptitle([setDescrpt ' -  Recruitment Curves']);
set(tmptit, 'Interpreter', 'none');
% Pix_SS = get(0,'screensize');
% hF3.Position = Pix_SS;