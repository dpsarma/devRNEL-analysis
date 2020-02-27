Muscle_ID = {'Right VM', 'Right RF', 'Right VL', 'Right BF', 'Right ST', ...
    'Right TA', 'Right SO', 'Right LG', 'Left VM', 'Left RF', 'Left VL',...
    'Left BF', 'Left ST', 'Left Ham', 'Left SO', 'Left LG'};
path_Save = 'R:\users\dsarma\';
path_UH3 = 'R:\data_raw\human\uh3_stim\';
subjectName = 'LSP02b';
reportPath = ['R:\users\dsarma\' subjectName '\'];

% % for hIdx = 1:69
% %     if ~isfield(hash(hIdx).standing,'stimAmp')
% %        continue;
% %     end
% %     stimAmps = [hash(hIdx).standing.stimAmp];
% %     [sortedStim,I] = sort(stimAmps);
% %     
% %     for m = hash(hIdx).Indices.muscles
% %         for ee = 1:length(hash(hIdx).Indices.standing) 
% % 
% %             X(ee,:) = hash(hIdx).standing(I(ee)).epochTimeVec;
% %             Y(ee,:) = abs(hash(hIdx).standing(I(ee)).meanTrace{m});
% %         %     Z(ee,:) = stimAmp;
% %         end
% % 
% %         hW(m) = figure;
% %         maximize;
% %         waterfall(Y,X)
% % %         xticklabels([-0.5 0 0.05 0.1 0.15]);
% %         xticklabels([-50 -33 -16 0 16 33 50 67 83 100]);
% %         t =sortedStim(1):(sortedStim(end)-sortedStim(1))/9:sortedStim(end);
% %         yticklabels([round(t/1000,2)]);
% % %         yticklabels([sortedStim]);
% % 
% %         xlabel('Time (ms)');
% %         ylabel(['Stim Amplitude (mA)']);
% %         zlabel('Voltage(uV)');
% %         view(22.5,26.8);
% %         hW(m).CurrentAxes.FontSize = 20;
% %         tmpsup = suptitle([Muscle_ID{m} ' - Day ' num2str(hash(hIdx).day) ' - Elec ' num2str(hash(hIdx).elecNum) ' - Standing -  Waterfall']);
% %         set(tmpsup, 'FontSize', 30);
% %         
% %         disp('Saving to Set Folder')
% %         savefig(hW(m),[reportPath Muscle_ID{m} '_Waterfalls_Day' num2str(hash(hIdx).day) '_Elec' num2str(hash(hIdx).elecNum) '_Standing_' num2str(hIdx)]);
% %         saveas(hW(m),[reportPath Muscle_ID{m} '_Waterfalls_Day' num2str(hash(hIdx).day) '_Elec' num2str(hash(hIdx).elecNum) '_Standing_' num2str(hIdx) '.png']);
% %         
% %         clear X Y
% %     end
% %     close all;
% % end


for hIdx = 10%11:69
    stimAmps = [hash(hIdx).sitting.stimAmp];
    [sortedStim,I] = sort(stimAmps);
    
    for m = hash(hIdx).Indices.muscles
        for ee = 1:length(hash(hIdx).Indices.sitting) 

            X(ee,:) = hash(hIdx).sitting(I(ee)).epochTimeVec;
            Y(ee,:) = abs(hash(hIdx).sitting(I(ee)).meanTrace{m});
        %     Z(ee,:) = stimAmp;
        end

        hW2(m) = figure;
        maximize;
        waterfall(Y,X)
        xticklabels([-50 -33 -16 0 16 33 50 67 83 100]);
%         xticklabels([-0.5 0 0.05 0.1 0.15]);
        t =sortedStim(1):(sortedStim(end)-sortedStim(1))/9:sortedStim(end);
        yticklabels([round(t/1000,2)]);
%         yticklabels([sortedStim]);

        xlabel('Time (sec)');
        ylabel(['Stim Amplitude (mA)']);
        zlabel('Voltage(uV)');
        view(22.5,26.8);
        hW2(m).CurrentAxes.FontSize = 20;
        tmpsup = suptitle([Muscle_ID{m} ' - Day ' num2str(hash(hIdx).day) ' - Elec ' num2str(hash(hIdx).elecNum) ' - Sitting -  Waterfall']);
        set(tmpsup, 'FontSize', 30);
        
        disp('Saving to Set Folder')
        savefig(hW2(m),[reportPath Muscle_ID{m} '_Waterfalls_Day' num2str(hash(hIdx).day) '_Elec' num2str(hash(hIdx).elecNum) '_Sitting_' num2str(hIdx)]);
        saveas(hW2(m),[reportPath Muscle_ID{m} '_Waterfalls_Day' num2str(hash(hIdx).day) '_Elec' num2str(hash(hIdx).elecNum) '_Sitting_' num2str(hIdx) '.png']);
        
        clear X Y
    end
    close all;
end