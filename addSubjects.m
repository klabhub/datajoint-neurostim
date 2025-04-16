function allSubjects = addSubjects(newSubjects,pv)
%% Merge an existing subjects json file with a new table/file containing subject data
arguments
    newSubjects
    pv.NS_ROOT= getenv("NS_ROOT")
    pv.save  (1,1) logical = false
    pv.force (1,1) logical  =false
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
oldSubjects = readJson(jsonFilename);
oldSubjects = struct2table(oldSubjects);

%% Check for duplicates
[newDup,oldDup]=ismember(newSubjects.subject,oldSubjects.subject);
if any(newDup)
    % Check full row match
    dupsOld = oldSubjects(oldDup(oldDup~=0),:);
    dupsNew = newSubjects(newDup,:);
    matchingVars = intersect(dupsOld.Properties.VariableNames,dupsNew.Properties.VariableNames);
    mismatch = ~ismember(dupsNew(:,matchingVars),dupsOld(:,matchingVars));
    if any(mismatch)    
        % Warn if some subjects get new, updated values.
        fprintf('\n New: ' ) 
        dupsNew(mismatch,:)
        fprintf('\n Old: ' ) 
        dupsOld(ismember(dupsOld.subject,dupsNew.subject(mismatch)),:)
        if ~pv.force
            % Don't update 
            error('Mismatched entries')
        else
            fprintf('Force selected: new will replace old.\n')
        end
    end   
end

% Merge new information on existing subjects
newVarsNotInOld = setdiff(newSubjects.Properties.VariableNames,oldSubjects.Properties.VariableNames);
oldVarsNotInNew = setdiff(oldSubjects.Properties.VariableNames,newSubjects.Properties.VariableNames);
% For subjects in old and new; merge information with new overruling old
mergedSubjects = innerjoin(oldSubjects,newSubjects(newDup,:),"Keys","subject","LeftVariables",oldVarsNotInNew,"RightVariables",newSubjects.Properties.VariableNames);
% For subjects in old only ; add default values for any newVarsNotInOld
oldOnly = oldSubjects(ismember(setdiff(oldSubjects.subject,newSubjects.subject),oldSubjects.subject),:);
if ~isempty(oldOnly) && ~isempty(newVarsNotInOld)
    def = repmat("",height(oldOnly),1);
    def = repmat({def},[1 numel(newVarsNotInOld)]);
    oldOnly = addvars(oldOnly,def{:},'NewVariableNames',newVarsNotInOld);    
end
% For Subjects in new only: add default valuse for any oldVarsNotInNew
newOnly = newSubjects(~newDup,:);
if ~isempty(newOnly) && ~isempty(oldVarsNotInNew)
    def = repmat("",height(newOnly),1);
    def = repmat({def},[1 numel(oldVarsNotInNew)]);
    newOnly = addvars(newOnly,def{:},'NewVariableNames',oldVarsNotInNew);    
end

allSubjects = [oldOnly;mergedSubjects;newOnly];

if pv.save
    json = table2struct(allSubjects);
    writeJson(json,jsonFilename);
end

