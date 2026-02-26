%{
# A class representing a Neurostim parameter/property
-> ns.Plugin
property_name : varchar(25)     # The name of the neurostim.parameter
---
property_value = NULL : longblob    # The value(s) of the constant/trial parameter/event 
property_time = NULL : longblob     # Trial time at which the event occured
property_nstime = NULL : longblob     # Neurostim experiment time at which the event occured
property_trial = NULL : longblob     # Trial in which the event occured
property_type : enum('Global','Parameter','Event','ByteStream') # Type of parameter
%}
% Global type have a single value in an experiment
% Parameter type have a single value per trial
% Event type have values that (can) change inside a trial. The time at which they occur is usually most relevant.
% ByteStream type are parameters that could not be stored "as is" in the
% database (e.g. objects) - They are converted to a byte stream and then
% stored.
% See Experiment.get how parameters can be retrieved from the
% database.
%
% BK - April 2022
classdef PluginParameter < dj.Part
    properties (SetAccess = protected)
        master = ns.Plugin;  % Part  table for the plugin
    end


    methods (Access=public)
        function addNew(tbl,key,name,value,type,trialTime,trial,nsTime,replace)
            arguments
                tbl (1,1) ns.PluginParameter
                key (1,1) struct % Key of the parent plugin
                name (1,:) char
                value
                type
                trialTime = []
                trial = []
                nsTime = []
                replace =false
            end

            key.property_name = name;
            if exists(tbl & key) && replace
                delQuick(tbl&key);
            end
            %% Add the data
            key.property_value = value;
            key.property_type = type;
            key.property_time = trialTime;
            key.property_trial = trial;
            key.property_nstime = nsTime;
            encodeAndInsert(tbl,key)

        end
    end

    methods (Access= ?ns.Plugin)
        function make(tbl,key,prms)
            % Called from ns.Plugin to fill the parameter values
            prmNames = fieldnames(prms);
            nrPrms = numel(prmNames);

            if isempty(key) || nrPrms ==0 ;return;end

            key = repmat(key,[1 nrPrms]); % replicate the plugin pkey part
            fprintf('Adding %d parameters for %s plugin\n',nrPrms,key(1).plugin_name);
            for i=1:nrPrms
                thisPrm = prms.(prmNames{i});
                if thisPrm.cntr==1 || ( thisPrm.cntr ==2 && isempty(thisPrm.log{1}))
                    % Single value: global property
                    % Easy: store only this value
                    type = 'Global';
                    value =thisPrm.value;
                    nsTime = thisPrm.time(end);
                    trial = 1;
                    time = -inf;
                elseif thisPrm.changesInTrial
                    % Changes within a trial : event
                    % Store all values, plus the time at which the event
                    % occurred.
                    type = 'Event';
                    [value,trial,time,nsTime] =get(thisPrm,'matrixIfPossible',true,'cellstrAsString',true,'removeFirstEmpty',true); %neurostim.parameter.get
                else
                    % One value per trial : parameter
                    % Some could be single key presses. So filling in across trials is not always right.
                    % So pull all values,just like events.
                    type = 'Parameter';
                    [value,trial,time,nsTime] = get(thisPrm,'matrixIfPossible',true,'cellstrAsString',true,'removeFirstEmpty',true);
                end


                % Some cleanup to store the values in the database
                if iscell(value)
                    %Replace function handles with strings
                    isFun = cellfun(@(x) (isa(x,'function_handle')),value);
                    [value(isFun)] = cellfun(@func2str,value(isFun),'UniformOutput',false);
                end
                % Neurostim should only store changing values, but that
                % did not always work perfectly. Here we detect values that
                % were really constant so that we can store a single value
                % for the experiment (i.e. Global type).
                if iscellstr(value)  || ischar(value) || isstring(value) || isnumeric(value) || islogical(value)
                    if ischar(value)
                        uValue = value;
                    elseif iscellstr(value) || isstring(value)
                        uValue=  unique(value);
                        if isscalar(uValue) && iscellstr(value) %#ok<ISCLSTR>
                            uValue= uValue{1}; % Get rid of cell
                        end
                        uValue = char(uValue);
                    elseif ismatrix(value)
                        uValue =unique(value);
                        % Dont do this: some events have nan values and only store the
                        % time
                        % if isnumeric(value) && all(isnan(value(:)))
                        % A vector with all nans
                        %    uValue =NaN;
                        %end
                    else % it is something >2D
                        uValue = nan(2,1); % Just a flag to skip the next part
                    end
                    if  size(uValue,1)==1
                        % Really only one value
                        value = uValue;
                        type = 'Global';
                    end
                end


                if isstruct(value) && numel(fieldnames(value))==0
                    value = true(size(value));
                end

                key(i).property_name = thisPrm.name;
                key(i).property_value= value;
                key(i).property_type = type;
                key(i).property_time = time;
                key(i).property_trial = trial;
                key(i).property_nstime = nsTime;
            end

            encodeAndInsert(tbl,key);
        end
    end


    methods (Access=protected)
        function encodeAndInsert(tbl,key)
            names = string({key.property_name});
            if numel(unique(upper(names))) ~= numel(key)
                % Two properties match case insensitively; that will cause
                % problems in the MySql database. The one with lower case
                % is renamed to X_lowercase
                isDouble = (names==lower(names)) & sum(upper(names) == upper(names)')==2;
                newName = char(names(isDouble) + "_lowercase");
                fprintf(2,'Properties of %s match case insentitively. %s  renamed to %s\n',key(1).plugin_name,names(isDouble),newName)
                key(isDouble).property_name = newName;
            end
            for i=1:numel(key)
                if islogical(key(i).property_value)
                    key(i).property_value = double(key(i).property_value);
                end
                if ~(isnumeric(key(i).property_value) || isstruct(key(i).property_value) || ischar(key(i).property_value))
                    % Database cannot not store this value. Convert to byte stream, the user can
                    % get the value by using getArrayFromByteStream,at
                    % least in Matlab (see Experiment.get),
                    key(i).property_value = getByteStreamFromArray(key(i).property_value);
                    key(i).property_type = 'ByteStream';
                end
            end
            insert(tbl,key)
        end
    end

    methods (Access= public)

        function G = get(tbl,pv)
            % Retrieve parameter values in a format that is used by
            % ns.Experiment/get() to return to the user
            arguments
                tbl (1,1) ns.PluginParameter
                pv.nwbRoot = []
                pv.prm {mustBeText} = string.empty
                pv.atTrialTime (1,1) double = NaN
                pv.trial (1,:) double = []
            end

            if ~exists(tbl);G=table;return;end

            if ~isempty(pv.prm)
                restriction  = sprintf('property_name ="%s"',pv.prm);
                if (~isempty(pv.trial) || ~isnan(pv.atTrialTime))
                    % A selection is in place; we'll need firstFrame from cic.
                    restriction = [restriction  ' OR property_name="firstframe"'];
                end
                tbl = tbl & restriction;
                if ~exists(tbl)
                    fprintf('The plugin does not have the property %s\n',pv.prm);
                    G=table;
                    return;
                end
            end

            T = fetchtable(tbl ,'*');


            % Group at the experiment level
            [ix,G ] = findgroups(T(:,["subject" "session_date" "starttime"]));
            S =splitapply(@(plg,name,value,time,nstime,trial,type) {ns.PluginParameter.nestedExperimentTable(plg,name,value,time,nstime,trial,type,pv)},T(:,["plugin_name" "property_name" "property_value" "property_time" "property_nstime" "property_trial" "property_type"]),ix);
            plgNames = unique(T.plugin_name);
            nrPlgs =numel(plgNames);
            for plg= 1:nrPlgs        
                G = addvars(G,cell(height(G),1),'NewVariableNames',plgNames(plg));
                for g= 1:height(G)
                    if ismember(plgNames(plg),S{g}.Properties.VariableNames) 
                        G{g,plgNames(plg)} = {S{g}.(plgNames(plg))};
                    end
                end
            end

            if ~isempty(pv.nwbRoot)
                % Only called for nwb export - one experiment/plugin at a time.
                % Determine when cic was constructedby finding the first
                % entry in nstime for the trial property. This defines time
                % zero in NWB.
                assert(nrPlgs==1,"More than one experiment in get for nwbRoot?");
                tpl=fetch(ns.PluginParameter  & (ns.Experiment &tbl) & 'plugin_name="cic"' & 'property_name="trial"' ,'property_nstime');
                timeZero  = tpl.property_nstime(1);
                names =fetchn(tbl,'property_name');

                for j=1:numel(names)
                    name =names{j};
                    data= G{1,plgNames}{1}.(name);
                    if ~(isnumeric(data) || islogical(data))
                        fprintf('Skipping %s in %s (non numeric data)\n',name,plgName);
                    else
                        timestamps = G{1,plgNames}{1}.([name 'NsTime']);
                        if isempty(timestamps)
                            timestamps =0; %Global property-pretend it was set at t=0.
                        else
                            % Align time stamps to the first GetSecs
                            % call,convert to seconds
                            timestamps =(timestamps-timeZero)/1000;
                        end
                        if size(data,1)== size(timestamps,1)
                            description = sprintf('%s property in %s plugin %d entries',name,plgNames,size(data,1));
                            ts =  types.core.TimeSeries('description',description, 'data',data','data_continuity','step','timestamps',timestamps,'data_unit','notspecified');
                            nwbRoot.stimulus_presentation.set(sprintf('%s_%s',plgNames ,name),ts);
                        else
                            fprintf('Skipping %s (mismatched timestamps)\n',name);
                        end
                    end
                end
            end
        end
    end

    methods (Static)
        function T = nestedExperimentTable(plg,name,value,time,nstime,trial,type,pv)
            % This is called from get with splitapply, ensuring that all of
            % the properties correspond to a single experiment. Here we
            % convert these to a nested table format, with one column per
            % plugin,
            [ix,grp] = findgroups(plg); % Find the plugins for this expt.
            S = splitapply(@(name,value,time,nstime,trial,type) {ns.PluginParameter.nestedPluginTable(name,value,time,nstime,trial,type)},name,value,time,nstime,trial,type,ix);
            T= table(S{:},'VariableNames',cellstr(grp));
            if ~isempty(pv.prm) && (~isnan(pv.atTrialTime) || ~isempty(pv.trial))
                if isscalar(unique(plg))
                    plg = "cic";
                    if ~isfield(T.cic,pv.prm);return;end
                else
                    plg = setdiff(plg,"cic");
                end
                T.(plg) = ns.PluginParameter.attrialtime(T.(plg),pv.prm,T.cic.firstframe.data,pv.atTrialTime,pv.trial);
            end
            if ismember("cic",T.Properties.VariableNames)
                if plg~="cic"
                    T = removevars(T,"cic");
                elseif ~isempty(pv.prm) && lower(pv.prm) ~="firstframe" && isfield(T.cic,'firstframe')
                    T.cic = rmfield(T.cic,"firstframe");
                end
            end
        end

        function S = nestedPluginTable(name,value,time,nstime,trial,type)
            % Called from splitapply in nestedExperimentTable to create a
            % struct per plugin.

            if ~iscell(value)
                value= {value};
                trial = {trial};
                time = {time};
                nstime = {nstime};
            end
            for prm = 1:numel(name)
                nm = lower(name(prm)); % Force lower case
                if strcmpi(type,'ByteStream')
                    % Convert from bytestrem to matlab values
                    value ={getArrayFromByteStream(value{prm})};
                end
                if iscellstr(value{prm}) || ischar(value{prm})
                    thisValue = string(value{prm});
                else
                    thisValue = value{prm};
                end
                S.(nm) = struct('data',thisValue,'trialtime', time{prm},'trial',trial{prm},'clocktime',nstime{prm});
            end
        end

        function props = attrialtime(props,propName,firstFrame,time,trial)
            % This is called from nestedExperimentTable to deal with the "feature" that Neurostim
            % stores only changes to a property, hence if a property did not change in
            % trial n, it should have the value it got in trial n-1.
            %
            % INPUT
            % props - struct with properties
            % propName - The name of the property that should be returned.
            % firstFrame -  the ns time of the first frame in each trial
            % time - The time (relative to the trial start) at which the property
            %           should be determined. Use Inf for 'at the end of the trial'.
            % trial - a vector of trial numbers for which to return the information.
            % Defaults to [], which means all trials.
            %
            % OUTPUT
            % props = An updated struct, one for each requested trial.
            %
            % BK - Feb 2026.
            arguments
                props (1,1) struct
                propName (1,1) string
                firstFrame (1,:) double
                time (1,1) double  = NaN % Select a time in the trial
                trial (1,:) double =[]  % Select trials
            end

            propName = lower(propName); % All properties are lower case
            nrTrials  = numel(firstFrame);
            allEventValues = props.(propName).data;
            allEventTrials = props.(propName).trial;
            allEventTimes  = props.(propName).trialtime;
            allEventNsTimes  = props.(propName).clocktime;
            if isnan(time)
                % No time selection (but a trial selection)
                [keepTrial,loc] = ismember(allEventTrials,trial);
                if iscell(allEventValues)
                    thisValue = allEventValues{keepTrial};                    
                else
                    thisValue = allEventValues(keepTrial);
                end
                props.(propName) = struct('data',thisValue,'trialtime', allEventTimes(keepTrial),'trial',allEventTrials(loc(keepTrial)),'clocktime',allEventNsTimes(keepTrial)); 
            else % Time selection : 1 per trial (optionally followed by trial selection)
                % Initialize with nan
                data = cell(nrTrials,1);
                eventTime = nan(nrTrials,1);
                eventNsTime = nan(nrTrials,1);
                [data{:}] = deal(NaN);
                currentTrial =NaN;
                % Loop through all events (= events when the property changed value)
                % For many events this could be made more efficient by skipping successive
                % events in a given trial, but the find needed for this is probably slower
                % than this loop over all events
                for e=1:numel(allEventTrials)
                    if ~isnan(currentTrial) && allEventTrials(e) > currentTrial+1
                        % Fill in skipped trials with the currentValue
                        trgTrials = currentTrial+1:allEventTrials(e)-1;
                        [data{trgTrials}] =deal(currentValue);
                        eventTime(trgTrials) = -inf; % Indicates that this value was set before the start of the trial
                        eventNsTime(trgTrials) = currentNsTime;
                    end
                    % Next trial or same trial, update until atTrialTime reached.
                    currentTrial = allEventTrials(e);
                    currentTime = allEventTimes(e);
                    currentNsTime = allEventNsTimes(e);
                    if iscell(allEventValues)
                        currentValue = allEventValues{e};
                    else
                        currentValue = allEventValues(e);
                    end
                    if currentTime <= time
                        data{currentTrial} = currentValue;
                        eventTime (currentTrial) = currentTime;
                        eventNsTime(currentTrial) = currentNsTime;
                    elseif isnan(eventTime(currentTrial))
                        % The current event occurred after the atTrialTime and
                        % there was no event in the current trial; use the value from the
                        % previous trial
                        if currentTrial >1
                            data{currentTrial} = data{currentTrial-1};
                            eventTime (currentTrial) = -inf;
                            eventNsTime(currentTrial) = eventNsTime(currentTrial-1);
                        end
                    end
                end
                % Fill in to the end from the last value that was stored.
                if currentTrial~=nrTrials
                    [data{currentTrial+1:nrTrials}] =deal(currentValue);
                    eventNsTime(currentTrial+1:nrTrials) = currentNsTime;
                    eventTime(currentTrial+1:nrTrials)  =-inf;
                end

                if ~isempty(trial)
                    trialNrs = find(~isinf(eventNsTime));
                    keepTrial = ismember(trialNrs,trial);
                else
                    keepTrial = true(numel(eventNsTime),1);
                end
                trialsWithTheEvent = intersect(find(~isinf(eventNsTime)),find(keepTrial)); % Trials in which the event actually ocurred
                props.(propName) = struct('data',{data(keepTrial)},'trialtime', eventTime(keepTrial),'trial',trialsWithTheEvent,'clocktime',eventNsTime(keepTrial));     
            end
            % Try to convert to matrix
            if iscell( props.(propName).data )
                if all(cellfun(@numel, props.(propName).data )==1)
                    props.(propName).data  = cell2mat( props.(propName).data );
                end
            elseif iscellstr(props.(propName).data) 
                props.(propName).data = string(props.(propName).data);
            end
        end
    end

end





