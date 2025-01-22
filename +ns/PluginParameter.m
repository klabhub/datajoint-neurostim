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
                trialTime
                trial
                nsTime
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
                    if ischar(uValue) || size(uValue,1)==1
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

        function v = get(tbl,nwbRoot,plgName)
            arguments
                tbl (1,1) ns.PluginParameter
                nwbRoot   = []
                plgName (1,:) char =''
            end
            % Retrieve all properties in the table as a struct
            % Used by ns.Experiment.get
            % Returns a struct with one field per property.
            persistent warnedAlready
            if isempty(warnedAlready);warnedAlready ={};end
            %% First the Global consts.
            [vals,names,nsTimes] = fetchn(tbl & 'property_type=''Global''' ,'property_value','property_name','property_nstime');
            % Create the struct with name/value
            glbl = cell(1,2*numel(names));
            [glbl{1:2:end}] =deal(names{:});
            [glbl{2:2:end}] = deal(vals{:});
            v  = struct(glbl{:});
            % And add the parmNsTime field for globals as well. Not usually
            % needed, except for a 1 trial experiment (to define
            % firstFrameNsTime for alignment)
            names= strcat(names,'NsTime');
            vals = nsTimes;
            glbl = cell(1,2*numel(names));
            [glbl{1:2:end}] =deal(names{:});
            [glbl{2:2:end}] = deal(vals{:});
            v= orderfields(mergestruct(v,struct(glbl{:})));
            %% Now the parms that change
            % Parameters - they do not change within a trial. The
            % output struct will have a vector/cell with one value for
            % each trial

            % Events - these can happen at any time. The struct
            % contains both the values and the times at which they
            % occurred (e.g. v.X and v.XTime)

            %Bytestream - can contain objects, coded as bytes.
            % Decode here.
            [vals,names,times,nsTimes,trials,parmTypes] = fetchn(tbl - 'property_type =''Global''' ,'property_value','property_name','property_time','property_nstime','property_trial','property_type');
            for j=1:numel(names)
                if any(cellfun(@(x) isfield(v,x),{names{j},[names{j} 'Trial'],[names{j} 'Time'],[names{j} 'NsTime']}))
                    oldName = names{j};
                    names{j} = [names{j} '2'];
                    if ~ismember(oldName,warnedAlready)
                        warnedAlready = cat(2,warnedAlready,{oldName});
                        warnNoTrace('%s already defined; renamed to %s',oldName,names{j});
                    end
                    % This happes because block and
                    % blockTrial are both cic properties. Unlikely to
                    % happen anywhere else.
                end


                name = names{j};
                if strcmpi(parmTypes(j),'ByteStream')
                    v.(name) =getArrayFromByteStream(vals{j});
                else
                    v.(name) =vals{j};
                end
                v.([name 'Time']) = times{j};
                v.([name 'NsTime']) = nsTimes{j};
                v.([name 'Trial']) = trials{j};
            end

            if ~isempty(nwbRoot)
                % Determine when cic was constructedby finding the first
                % entry in nstime for the trial property.
                tpl=fetch(ns.PluginParameter  & (ns.Experiment &tbl) & 'plugin_name="cic"' & 'property_name="trial"' ,'property_nstime');
                timeZero  = tpl.property_nstime(1);
                for j=1:numel(names)
                    name =names{j};
                    data= v.(name);
                    if ~(isnumeric(data) || islogical(data))
                        fprintf('Skipping %s in %s (non numeric data)\n',name,plgName);
                    else
                        timestamps = v.([name 'NsTime']);
                        if isempty(timestamps)
                            timestamps =0; %Global property-pretend it was set at t=0.
                        else
                            % Align time stamps to the first GetSecs
                            % call,convert to seconds
                            timestamps =(timestamps-timeZero)/1000;
                        end
                        if size(data,1)== size(timestamps,1)
                            description = sprintf('%s property in %s plugin %d entries',name,plgName,size(data,1));
                            ts =  types.core.TimeSeries('description',description, 'data',data','data_continuity','step','timestamps',timestamps,'data_unit','notspecified');
                            nwbRoot.stimulus_presentation.set(sprintf('%s_%s',plgName ,name),ts);
                        else
                            fprintf('Skipping %s (mismatched timestamps)\n',name);
                        end
                    end
                end
            end
        end
    end


end



