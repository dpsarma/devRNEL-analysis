%% UH3 - LSP02 anaylsis (1) - MDF
%----------------------------------------
% tag that is missing from all EMG trials is "Sitting" vs. "Standing"

%% Initialization & Data Prep
%
% IMPORTANT: make sure to change confSel to the name of your DC as it
% appears in the mdf configuration file
%
omdfc = mdf.init(struct( ...
  'confFile' , 'auto', ...
  'confSel'  , '<Data collection name>'));

summary =  mdf.load('mdf_type','summary');
musclabels = summary.emgLabel;

dates = {'2018-10-08'; '2018-10-09'; '2018-10-11'; '2018-10-16'; '2018-10-19'; '2018-10-24'};
%dayIDs = ['Day 02'; 'Day 03'; 'Day 04'; 'Day 07'; 'Day 09'; 'Day 12'];
%
% I would rewrite this block of code as below
%electrodenums = [1:48];
%for i=1:length(electrodenums)
%    if (electrodenums(i)<10)
%        electrodestrings{i} = ['Unipolar: 0' num2str(electrodenums(i))];
%    else
%        electrodestrings{i} = ['Unipolar: ' num2str(electrodenums(i))];
%    end 
%end
electrodes = arrayfun(@(en) sprintf('Unipolar: %02d',en),[1:48]);

            
%% Create list of data for loading
%
% if you do not need uuids by date and electrode, leave them empty or do
% not pass them in
%
groupedUuids = findTrialsUuid(dates,electrodes);

% flatten the structure
temp = arrayfun(@(x) x.uuids',groupedUuids,'Uniform',0)';
uuids = [temp{:}];

%% Loop for each trial
%
for i = 1:length(uuids)
    
    res = dataPrepSingle(uuids{i});
    
end %for
