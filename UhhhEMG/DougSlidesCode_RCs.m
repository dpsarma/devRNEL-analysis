%% RMS Recruitment Curves
for hIdx = 1:5
    stimAmps = [hash(hIdx).standing.stimAmp];
    [sortedStim,I] = sort(stimAmps);
    
    %% Plotting Recruitment Curves
    responseRMS = [hash(hIdx).standing(:).responseRMS]';
    responseRMSErr = [hash(hIdx).standing(:).responseRMSErr]';
    % response = [trial(:).responseP2P]';

    hF3 = figure; %maximize;
    for m = hash(hIdx).Indices.muscles
        yyA = smooth(sortedStim/1000,responseRMS(I,m),'sgolay');
        errorbar(sortedStim/1000,yyA,responseRMSErr(I,m),'Color',C(m,:));%,'LineWidth',3);
        hold on;
    end
    ylabel('(post-Stim) Mean RMS  (uV)')
    xlabel('Stimulation Amplitude (mA)')
%         title(['By Amplitude']);

    legend([hash(hIdx).muscleNames],'Location','northwest', 'Orientation','vertical');


    tmptit = suptitle(['Day ' num2str(hash(hIdx).day) ' - Elec ' num2str(hash(hIdx).elecNum) ' -  Standing - Amplitude RC (RMS)']);
    set(tmptit, 'Interpreter', 'none');
        % % Pix_SS = get(0,'screensize');
%         hF3.Position = [680 49 572 947];

        saveas(hF3,[reportPath 'RecruitmentCurvesRMS_Day'  num2str(hash(hIdx).day) '_Elec' num2str(hash(hIdx).elecNum) '_Standing.png']);
        savefig(hF3,[reportPath 'RecruitmentCurvesRMS_Day'  num2str(hash(hIdx).day) '_Elec' num2str(hash(hIdx).elecNum) '_Standing' ]);


    stimAmps2 = [hash(hIdx).sitting.stimAmp];
    [sortedStim2,I2] = sort(stimAmps2);
    
    %% Plotting Recruitment Curves
    responseRMS2 = [hash(hIdx).sitting(:).responseRMS]';
    responseRMS2Err = [hash(hIdx).sitting(:).responseRMSErr]';
    % response = [trial(:).responseP2P]';

    hF3 = figure; %maximize;
    for m = hash(hIdx).Indices.muscles
        yyA = smooth(sortedStim2/1000,responseRMS2(I2,m),'sgolay');
        errorbar(sortedStim2/1000,yyA,responseRMS2Err(I2,m),'Color',C(m,:));%,'LineWidth',3);
        hold on;
    end
    ylabel('(post-Stim) Mean RMS  (uV)')
    xlabel('Stimulation Amplitude (mA)')
%         title(['By Amplitude']);

    legend([hash(hIdx).muscleNames],'Location','northwest', 'Orientation','vertical');


    tmptit = suptitle(['Day ' num2str(hash(hIdx).day) ' - Elec ' num2str(hash(hIdx).elecNum) ' -  Sitting - Amplitude RC (RMS)']);
    set(tmptit, 'Interpreter', 'none');
        % % Pix_SS = get(0,'screensize');
%         hF3.Position = [680 49 572 947];

       saveas(hF3,[reportPath 'RecruitmentCurvesRMS_Day'  num2str(hash(hIdx).day) '_Elec' num2str(hash(hIdx).elecNum) '_Sitting.png']);
       savefig(hF3,[reportPath 'RecruitmentCurvesRMS_Day'  num2str(hash(hIdx).day) '_Elec' num2str(hash(hIdx).elecNum) '_Sitting' ]);

end

%% P2P Recruitment Curves 
for hIdx = 1:5
    stimAmps = [hash(hIdx).standing.stimAmp];
    [sortedStim,I] = sort(stimAmps);
    
    %% Plotting Recruitment Curves
    responseP2P = [hash(hIdx).standing(:).responseP2P]';
%     responseP2PErr = [hash(hIdx).standing(:).responseP2PErr]';
    % response = [trial(:).responseP2P]';

    hF3 = figure; %maximize;
    for m = hash(hIdx).Indices.muscles
        yyA = smooth(sortedStim/1000,responseP2P(I,m),'sgolay');
        plot(sortedStim/1000,yyA,'Color',C(m,:),'LineWidth',3);
        hold on;
    end
    ylabel('(post-Stim) Mean P2P  (uV)')
    xlabel('Stimulation Amplitude (mA)')
%         title(['By Amplitude']);

    legend([hash(hIdx).muscleNames],'Location','northwest', 'Orientation','vertical');


    tmptit = suptitle(['Day ' num2str(hash(hIdx).day) ' - Elec ' num2str(hash(hIdx).elecNum) ' -  Standing - Amplitude RC (P2P)']);
    set(tmptit, 'Interpreter', 'none');
        % % Pix_SS = get(0,'screensize');
%         hF3.Position = [680 49 572 947];

    saveas(hF3,[reportPath 'RecruitmentCurvesP2P_Day'  num2str(hash(hIdx).day) '_Elec' num2str(hash(hIdx).elecNum) '_Standing.png']);
    savefig(hF3,[reportPath 'RecruitmentCurvesP2P_Day'  num2str(hash(hIdx).day) '_Elec' num2str(hash(hIdx).elecNum) '_Standing' ]);


    stimAmps2 = [hash(hIdx).sitting.stimAmp];
    [sortedStim2,I2] = sort(stimAmps2);
    
    %% Plotting Recruitment Curves
    responseP2P2 = [hash(hIdx).sitting(:).responseP2P]';
%     responseP2P2Err = [hash(hIdx).sitting(:).responseP2PErr]';
    % response = [trial(:).responseP2P]';

    hF3 = figure; %maximize;
    for m = hash(hIdx).Indices.muscles
        yyA = smooth(sortedStim2/1000,responseP2P2(I2,m),'sgolay');
        plot(sortedStim2/1000,yyA,'Color',C(m,:),'LineWidth',3);
        hold on;
    end
    ylabel('(post-Stim) Mean P2P  (uV)')
    xlabel('Stimulation Amplitude (mA)')
%         title(['By Amplitude']);

    legend([hash(hIdx).muscleNames],'Location','northwest', 'Orientation','vertical');


    tmptit = suptitle(['Day ' num2str(hash(hIdx).day) ' - Elec ' num2str(hash(hIdx).elecNum) ' -  Sitting - Amplitude RC (P2P)']);
    set(tmptit, 'Interpreter', 'none');
        % % Pix_SS = get(0,'screensize');
%         hF3.Position = [680 49 572 947];

    saveas(hF3,[reportPath 'RecruitmentCurvesP2P_Day'  num2str(hash(hIdx).day) '_Elec' num2str(hash(hIdx).elecNum) '_Sitting.png']);
    savefig(hF3,[reportPath 'RecruitmentCurvesP2P_Day'  num2str(hash(hIdx).day) '_Elec' num2str(hash(hIdx).elecNum) '_Sitting' ]);

end