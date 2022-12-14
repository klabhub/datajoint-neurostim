%{
# A class representing a Neurostim parameter/property
-> ns.Plugin
property_name : varchar(25)     # The name of the neurostim.parameter
---
property_value = NULL : longblob    # The value(s) of the constant/trial parameter/event 
property_time = NULL : longblob     # Time at which the event occured
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

    methods (Access= ?ns.Plugin)
        function make(self,key,prm)
            % Called from ns.Plugin to fill the parameter values
            assert(isa(prm,'neurostim.parameter'),'PluginParameter needs a neurostim.parameter object as its input.');

            time =[];
            trial =[];
            if prm.cntr==1
                % Single value: global property
                % Easy: store only this value
                type = 'Global';
                value =prm.value;
            elseif prm.changesInTrial
                % Changes within a trial : event
                % Store all values, plus the time at which the event
                % occurred.
                type = 'Event';
                [value,trial,time] =get(prm,'matrixIfPossible',true);
            else
                % One value per trial : parameter
                % Some could be single key presses. So filling in across trials is not always right.
                % So pull all values,just like events.
                type = 'Parameter';
                [value,trial,time] = get(prm,'matrixIfPossible',true);
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
                if iscellstr(value) %#ok<ISCLSTR>
                    uValue=  unique(value);
                    if size(uValue,1)==1
                        uValue= uValue{1}; % Get rid of cell
                    end
                elseif ismatrix(value)
                    uValue=  unique(value,'rows');
                    if isnumeric(value) && all(isnan(value(:)))
                        % A vector with all nans
                        uValue =NaN;
                    end
                else % it is something >2D
                    uValue = nan(2,1); % Just a flag to skip the next part
                end
                if size(uValue,1)==1
                    % Really only one value
                    value = uValue;
                    type = 'Global';
                    time = [];
                end
            end


            if isstruct(value) && numel(fieldnames(value))==0
                value = true(size(value));
            end


            key.property_name = prm.name;
            key.property_value = value;
            key.property_type = type;
            key.property_time = time;
            key.property_trial = trial;
            try
                self.insert(key);
            catch me
                if contains(me.message,'Duplicate entry')
                    key.property_name = [key.property_name key.property_name];
                    self.insert(key)
                elseif contains(me.message,'Matlab placeholder') || contains(me.message,'unsupported type')
                    % Database could not store this value. Probably some
                    % kind of object. Convert to byte stream, the user can
                    % get the value by using getArrayFromByteStream if
                    % really needed,  at least in
                    % Matlab (see Experiment.get),
                    key.property_value = getByteStreamFromArray(value);
                    key.property_type = 'ByteStream';
                    self.insert(key);
                else
                    rethrow(me)
                end
            end
        end

    end

    methods (Access= public)

        function v = get(tbl)
            % Retrieve all properties in the table as a struct
            % Used by ns.Experiment.get
            % Returns a struct with one field per property.
            v= struct;
            %% First the Global consts.            
            [vals,names] = fetchn(tbl & 'property_type=''Global''' ,'property_value','property_name');
            glbl = cell(1,2*numel(names));
            [glbl{1:2:end}] =deal(names{:});
            [glbl{2:2:end}] = deal(vals{:});
            v  = struct(glbl{:});

            %% Now the parms that change
            % Parameters - they do not change within a trial. The
            % output struct will have a vector/cell with one value for
            % each trial

            % Events - these can happen at any time. The struct
            % contains both the values and the times at which they
            % occurred (e.g. v.X and v.XTime)

            %Bytestream - can contain objects, coded as bytes.
            % Decode here.
            [vals,names,times,trials,types] = fetchn(tbl - 'property_type =''Global''' ,'property_value','property_name','property_time','property_trial','property_type');
            for j=1:numel(names)
                if strcmpi(types(j),'ByteStream')
                    v.(names{j}) =getArrayFromByteStream(vals{j});
                else
                    v.(names{j}) =vals{j};
                end
                v.([names{j} 'Time']) = times{j};
                v.([names{j} 'Trial']) = trials{j};
            end

        end
    end


end



