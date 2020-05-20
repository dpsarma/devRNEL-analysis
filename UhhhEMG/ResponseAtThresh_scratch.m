%% Getting Response @ Threshold.
% Needs stim thresholds and response/stimAmps from Normal Quick Summary
% load('D:\DATA\UH3 testing\LSP05\data_gen\threshes\Sorted_lowest_threshByElec.mat')
for i=1:length(elec_labels)
    tmpLabels{i} = [cell2mat(elec_labels(i)) '_' num2str(lowestThresh(i))];
end

for l_set = 4:4
    switch l_set
        case 1
            disp('Loading Set 1')
            load('D:\DATA\UH3 testing\LSP05\data_gen\Day8-Set_1_RC_1-6mA_500us_1Hz.mat')
        case 2
             disp('Loading Set 2')
            load('D:\DATA\UH3 testing\LSP05\data_gen\Day8-Set_2_RC_1-6mA_500us_1Hz.mat')
        case 3
             disp('Loading Set 3')
            load('D:\DATA\UH3 testing\LSP05\data_gen\Day12-Set_1_E1-8RC_1-6mA_500us_1Hz.mat')
        case 4
             disp('Loading Set 4')
            load('D:\DATA\UH3 testing\LSP05\data_gen\day9-set1.mat')
    end

    for n = 1:length(elec_labels)
    elecs(n) = str2num(erase(cell2mat(elec_labels(n)),'E'));
    end
    response = [trial(:).responseRMS]';

    if ~exist('lowestThresh','var')
        lowestThresh = Y;
    end

    if ~exist('respAtThresh','var')
        disp('respAtThresh does not exist');
        respAtThresh = zeros(length(elecs),length(Muscle_ID));
        normResp = zeros(length(elecs),length(Muscle_ID));
    end

    liveElecs = zeros(length(elecs),1);
    for ii = 1:length(elecs)
        eIdx = find(spinalElec == elecs(ii));
    %     [sortedStim,I] = sort(stimAmp(eIdx));
    %     [sortedPulse,I2] = sort(pulseWidth(eIdx));
    %     [sortedFreq,I3] = sort(stimFrequency(eIdx));

        if ~isempty(eIdx)
            if respAtThresh(ii)==0
               tI = eIdx(stimAmp(eIdx)==lowestThresh(ii));
               respAtThresh(ii,:) = response(tI(end),:)-trial(tI(end)).baseline'; 
               disp('Updating Response');
               liveElecs(ii) = 1;
            else
                disp('Electrode Already Exists');
            end
        else
            disp('Electrodes Not Found in This Set');
        end
    end

    for ll =  1:length(trial)
        bases(ll,:) = trial(ll).baseline;
        for m =  1:length(Muscle_ID)
            YY(ll,m) = response(ll,m)-bases(ll,m);
        end
    end
    mBase = mean(bases);
    maxResp = max(YY);

    normResp(liveElecs==1,:) = round(respAtThresh(liveElecs==1,:)./maxResp,2);
    normResp(normResp < 0) = 0;

    %% Plotting
    h1 = figure;maximize;
    bar(normResp(liveElecs==1,:)');
    legend(tmpLabels(liveElecs==1), 'Location','northeastoutside','Interpreter', 'none');
    xticklabels(Muscle_ID);
    xticks(1:1:16);
    set(gca,'TickLength',[0 .01]);
    ylim([0 1])
    ylabel('% of Max Activation')
    title('Normalized RMS at lowest Threshold Across Muscles');

    respAtThresh(respAtThresh < 0) = 0;
    h2 = figure;maximize;
    bar(respAtThresh(liveElecs==1,:)');
    legend(tmpLabels(liveElecs==1), 'Location','northeastoutside','Interpreter', 'none');
    xticklabels(Muscle_ID);
    xticks(1:1:16);
    set(gca,'TickLength',[0 .01]);
    ylabel('RMS EMG (uV)')
    title('RMS at lowest Threshold Across Muscles');
    
    saveas(h1,['D:\FigRescources\UH3\LSP05\emgRecruitment_Summary\NormResponseAtThresh_' num2str(l_set) '.png']);
    saveas(h2,['D:\FigRescources\UH3\LSP05\emgRecruitment_Summary\ResponseAtThresh_' num2str(l_set) '.png']);

    clearvars -except respAtThresh lowestThresh elec_labels normResp l_set tmpLabels %liveElecs
    %% Think about plotting the highest response as well.
end
    % % RELOAD EVERYTHING AND SAVE THE NORM VARIABLE
        %% Plotting
    h3 = figure;maximize;
    bar(normResp');
    legend(tmpLabels, 'Location','northeastoutside','Interpreter', 'none');
    xticklabels(Muscle_ID);
    xticks(1:1:16);
    set(gca,'TickLength',[0 .01])
    ylim([0 1])
    ylabel('Norm EMG')
    title('Normalized RMS at lowest Threshold Across Muscles');
 respAtThresh(respAtThresh < 0) = 0;
    h4 = figure;maximize;
    bar(respAtThresh');
    legend(tmpLabels, 'Location','northeastoutside','Interpreter', 'none');
    xticks(1:1:16);
    xticklabels(Muscle_ID);
    set(gca,'TickLength',[0 .01])
    ylabel('RMS EMG (uV)')
    title('RMS at lowest Threshold Across Muscles');
   
    
    saveas(h3,['D:\FigRescources\UH3\LSP05\emgRecruitment_Summary\NormResponseAtThresh_All.png']);
    saveas(h4,['D:\FigRescources\UH3\LSP05\emgRecruitment_Summary\ResponseAtThresh_All.png']);
