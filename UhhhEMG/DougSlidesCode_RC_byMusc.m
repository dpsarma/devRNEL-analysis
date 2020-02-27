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

for m = 1:length(Muscle_ID)
    hF4(m) = figure; %maximize;
    cCount = 1;
    for hIdx = 1:69
        if ~contains(hash(hIdx).muscleNames, Muscle_ID(m))
            continue;
        elseif ~isfield(hash(hIdx).standing,'stimAmp')
            disp(['No Standing Found - ' num2str(hIdx)]);
            continue;
% %         elseif hIdx == 62 %|| hIdx == 52 || hIdx == 53
% %             disp('Skipping Bad Run');
% %             continue;
        else           
            stimAmps3 = [hash(hIdx).standing.stimAmp];
            [sortedStim3,I] = sort(stimAmps3);
            %% Plotting Recruitment Curves
            responseRMS = [hash(hIdx).standing(:).responseRMS]'; %RMS
            responseRMSErr = [hash(hIdx).standing(:).responseRMSErr]'; %RMS
            % response = [trial(:).responseP2P]';

            yyA = smooth(sortedStim3/1000,responseRMS(I,m),'sgolay');
            errorbar(sortedStim3/1000,yyA,responseRMSErr(I,m),'Color',C2(cCount,:),'LineWidth',2);
%             plot(sortedStim3/1000,yyA,'Color',C2(hIdx,:),'LineWidth',3);
            hold on;
            cCount = cCount + 1;
        end
        ssnName{hIdx,:} = ['Day' num2str(hash(hIdx).day) '_Elec' num2str(hash(hIdx).elecNum)];
    end
    
    ylabel('(post-Stim) Mean RMS  (uV)')
    xlabel('Stimulation Amplitude (mA)')
    title([Muscle_ID(m) ' - Standing - Amplitude RCs - RMS']);
    
    if exist('ssnName')
        tmpleg = legend([ssnName(~cellfun('isempty',ssnName))'],'Location','eastoutside', 'Orientation','vertical');
        set(tmpleg,'Interpreter', 'none');
        
    end
    
%     disp('Saving to Set Folder')
%         savefig(hF4(m),[reportPath Muscle_ID{m} '_RMS_Recruitment_standing']);
%         saveas(hF4(m),[reportPath Muscle_ID{m} '_RMS_Recruitment_standing.png']);
end


for m = 1:length(Muscle_ID)
        disp('Saving to Set Folder')
%         hF4(m).CurrentAxes.YLim = [0 500]
        savefig(hF4(m),[reportPath Muscle_ID{m} '_RMS_Recruitment_standing']);
        saveas(hF4(m),[reportPath Muscle_ID{m} '_RMS_Recruitment_standing.png']);
end
close all;

