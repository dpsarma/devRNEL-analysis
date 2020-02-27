function [uuids] = findTrialsUuid(dates,electrodes)
   
    % checks input argument
    if nargin < 1
      dates = {};
    end
    if nargin < 2
      electrodes = {};
    end

    % aggregation pipe
    aggregation = { ...
      '{$match:{"mdf_def.mdf_type":"trial","mdf_metadata.trialType" : "EMG"}}', ...
      '{$project:{"subject":"$mdf_metadata.subject","date":{$substr:["$mdf_metadata.time",0,10]},"electrode":"$mdf_metadata.electrodeLabel","uuid":"$mdf_def.mdf_uuid"}}', ...
      '{$group:{"_id":{"subject":"$subject","date":"$date","electrode":"$electrode"},"uuids":{$push:"$uuid"}}}', ...
      '{$project:{"_id":0,"subject":"$_id.subject","date":"$_id.date","electrode":"$_id.electrode","uuids":"$uuids"}}'...
    };

    % extract only dates and electrodes requested
    if isa(dates,'cell') && length(dates) > 0
      aggregation{end+1} = ['{$match:{"date":{$in:[' strjoin(cellfun(@(x) ['"' x '"'],dates,'UniformOutput',0),', ')  ']}}}'];
    end
    if isa(electrodes,'cell') && length(electrodes) > 0
      aggregation{end+1} = ['{$match:{"electrode":{$in:[' strjoin(cellfun(@(x) ['"' x '"'],electrodes,'UniformOutput',0),', ')  ']}}}']; 
    end

    dbObj = mdfDB.getInstance();
    raw = dbObj.aggregate(aggregation);

    uuids = struct('subject','','date','','electrode','','uuids','');
    for i = 1:length(raw)
       uuids(end+1).subject = raw{i}.subject;
       uuids(end).date      = raw{i}.date;
       uuids(end).electrode = raw{i}.electrode;
       uuids(end).uuids     = raw{i}.uuids;
    end %for
    uuids(1) = [];
end
