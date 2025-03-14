function mergedSubjects = addSubjects(newSubjects,pv)
%% Merge an existing subjects json file with a new table/file containing subject data
arguments
    newSubjects
    pv.NS_ROOT= getenv("NS_ROOT")
    pv.save  (1,1) logical = false
end

if isstring(newSubjects)
    assert(exist(newSubjects,"file"),'Subject file %s does not exist',newSubjects);
    newSubjects = readtable(newSubjects);
    
    if ismember("KLabNumber", newSubjects.Properties.VariableNames)
        % Export from KLab data base.
        newSubjects = renamevars(newSubjects,["KLabNumber" "DateOfBirth" "Gender" "Handedness"],["subject" "dob" "sex" "handedness"]);
        newSubjects = addvars(newSubjects,repmat("Homo sapiens",height(newSubjects),1),'NewVariableNames','species');
        newSubjects.dob = datetime(newSubjects.dob,"Format","uuuu-MM-dd");
        newSubjects  =convertvars(newSubjects,1:width(newSubjects),"string");
        [newSubjects{newSubjects.handedness=="","handedness"}] = deal("Unknown");
    end
end

jsonFilename = fullfile(pv.NS_ROOT,"subject.json");
assert(exist(jsonFilename,"file"),'Subject file %s does not exist',jsonFilename);
currentSubjects = readJson(jsonFilename);
currentSubjects = struct2table(currentSubjects);

%% Check for duplicates
[newDup,oldDup]=ismember(newSubjects.subject,currentSubjects.subject);
if any(newDup)
    % Check full row match
    dupsOld = currentSubjects(oldDup(oldDup~=0),:);
    dupsNew = newSubjects(newDup,:);
    matchingVars = intersect(dupsOld.Properties.VariableNames,dupsNew.Properties.VariableNames);
    mismatch = ~ismember(dupsNew(:,matchingVars),dupsOld(:,matchingVars));
    if any(mismatch)        
        for i=find(mismatch)
            fprintf('Old: ' ) 
            dupsOld(dupsOld.subject ==dupsNew.subject(i),:)
            fprintf('\n New: ' ) 
            dupsNew(i,:)
            fprintf('\n');
        end
        error('Mismatched entries')
    end   
end

% Merge new information on existing subjects
newVars = setdiff(newSubjects.Properties.VariableNames,currentSubjects.Properties.VariableNames);
mergedSubjects = innerjoin(currentSubjects,newSubjects(newDup,:),"Keys","subject","RightVariables",newVars);
% add completely new subjects
mergedSubjects= [mergedSubjects;newSubjects(~newDup,:)];
% And put a default value for the new information in the current subjects
% that were not updated
updateCurrent = currentSubjects(oldDup==0,:);
if ~isempty(updateCurrent)
for i=1:numel(newVars)   
    updateCurrent = addvars(updateCurrent,repmat("",height(updateCurrent),1),'NewVariableNames',newVars{i});
end
mergedSubjects = [mergedSubject;updateCurrent];
end



if pv.save
    json = table2struct(mergedSubjects);
    writeJson(json,jsonFilename);
end

