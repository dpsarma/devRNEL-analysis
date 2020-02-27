%% UH3 - LSP02 anaylsis (1) - MDF list load
%----------------------------------------
% tag that is missing from all EMG trials is "Sitting" vs. "Standing"


%% Initialization & Data File Prep
datasave = 'R:\users\dsarma\LSP02b\matData\mdftest\';
omdfc = mdf.init('auto');

summary =  mdf.load('mdf_type','summary');
musclabels = summary.emgLabel;

datestrings = ['2018-10-08'; '2018-10-09'; '2018-10-11'; '2018-10-16'; '2018-10-19'; '2018-10-24'];
dayIDs = ['Day 02'; 'Day 03'; 'Day 04'; 'Day 07'; 'Day 09'; 'Day 12'];
freqs = [1 2 5];
pulses = [200 250 400 500 600 750 800 1000];
electrodenums = [1:48];
for i=1:length(electrodenums)
            if (electrodenums(i)<10)
                electrodestrings{i} = ['Unipolar: 0' num2str(electrodenums(i))];
            else
                electrodestrings{i} = ['Unipolar: ' num2str(electrodenums(i))];
            end 
end
        
for ff = 2:length(freqs)
    freqstring = freqs(ff);

    for p = 1:length(pulses)
        pulsestring = pulses(p);        

        %% Create list of data for loading
        %ee = findTrialsPerDay_Elec('2018-10-08','Unipolar: 09')

        for f = 1:size(datestrings,1)
            for e = 1:length(electrodenums)
                % Create load table by date & electrode number
                disp(['Starting: ' dayIDs(f,:) ' ' electrodestrings(e) ' ' num2str(freqstring) 'Hz_' num2str(pulsestring) 'us'])
                uuid_list{e,f} = findTrialsPerAllVars(datestrings(f,:),char(electrodestrings(e)), freqstring, pulsestring);
% %                 if isempty(uuid_list{e,f})
% %                     disp(['No Trials for: ' dayIDs(f,:) ', ' electrodestrings(e)]);
% %                     continue;
% %                 end
            end
        end
        % % % OR create table for ALL EMG Trials by date
        % % uuid_list = findTrialsPerDay('2018-10-08');


        save(char([datasave 'mdf_uuidlist_' num2str(freqstring) 'Hz_' num2str(pulsestring) 'us.mat']),'uuid_list',...
            'summary','musclabels','dayIDs','datestrings','freqstring','pulsestring','-v7.3');
    end
end