function [allUUID] = findTrialsPerDay(dateString)

%     dateString = '2018-10-05';
    import com.mongodb.BasicDBObject
    import com.mongodb.BasicDBList
    
    match1 = BasicDBObject('mdf_def.mdf_type', BasicDBObject('$eq','trial'));    
    match2 = BasicDBObject('mdf_metadata.time', BasicDBObject('$regex',dateString));
    match3 = BasicDBObject('mdf_metadata.trialType', BasicDBObject('$eq','EMG'));
    matchList = BasicDBList(); 
    matchList.add(match1);
    matchList.add(match2);
    matchList.add(match3);
    matchAnd = BasicDBObject('$and', matchList);
    match  = BasicDBObject('$match',  matchAnd); 

    groupFields = BasicDBObject('_id','$mdf_metadata.subject');
    groupFields.put('uuid', BasicDBObject( '$push', '$mdf_def.mdf_uuid'));
    group = BasicDBObject('$group', groupFields);

    pipe = BasicDBObject.parse('{"_id": 1,"objUUID": "$uuid"}');
    project  = BasicDBObject('$project',  pipe); 

    aggrListFinal = BasicDBList();
    aggrListFinal.add(match); 
    aggrListFinal.add(group);
    aggrListFinal.add(project);

    dbInst = mdfDB.getInstance;
    tmp = dbInst.coll.aggregate(aggrListFinal);
    if ~strcmp(char(tmp.results),'[]')
        tmp2 = loadjson(char(tmp.results));
        allUUID = tmp2{1}.objUUID;
    else
        allUUID  = [];
    end
end