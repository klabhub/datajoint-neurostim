%{
# A dimension refers to a group of conditions (each of which represents a set of trials in an Experiment).
-> ns.Experiment 
-> ns.DimensionParm 
%}
%
% A 'Dimension' is a set of conditions. Each condition represents a set of
% trials in an experiment and is stored in the part tables DimensionCondition
% and DimensionTrial.
%
% For instance, when mapping an orientation tuning curve with the gabor plugin,
% the dimension could be called 'ori', and the trials in which the orientation
% was 30 degrees are stored in a DimensionCondition tuple with the name
% gabor:orientation:30.
%
% How a dimension is computed is defined in ns.DimensionParm. Each
% DimensionParm row specifies which plugin parameter(s) drive the conditions,
% together with optional restrict/exclude filters, a transform function, and
% other options. The paradigms for which a DimensionParm applies are stored in
% ns.DimensionParmParadigm.
%
%
% To register a new dimension definition, use ns.Dimension.defineParms or
% insert directly into ns.DimensionParm:
%   ns.Dimension.defineParms(expt, 'gabor', 'orientation', 'ori')
%
% The .name field in DimensionCondition creates unique names for the
% conditions. For instance, if one of the orientations was 10, its name
% would become gabor:orientation:10.
%
% The ns.Tuning table uses Dimension to determine tuning curves.
%
% See also ns.DimensionParm, ns.DimensionParmParadigm, ns.DimensionCondition,
%          ns.Tuning
classdef Dimension < dj.Computed & dj.DJInstance

    methods (Access=public)
        function varargout = split(d,T,pv)
            % Split a table T into subtables  corresponding to the conditions in the dimension.
            % This is done by matching the .trials in each condition to the
            % .trial column in the table. For instance, if the dimension has two
            % conditions, the table T will be split int two subtables
            % containing only the trials from each condition.  The second
            % output gives the values that identify the conditions.
            %
            % With labelOnly =true, the table is not split, but a column is added
            % to indicate which condition the row (trial) corresponds to.
            %
            arguments
                d (1,1) ns.Dimension {mustHaveRows(d,1)}
                T (:,:) table
                pv.labelOnly  (1,1) logical = true
            end
            conditions = ns.DimensionCondition & d;
            nrConditions  = count(conditions);
            dimensionName = fetch1(d,'dimension');
            if pv.labelOnly
                % Add a column that labels the condition in the
                % existing table
                TT = addvars(T,nan(height(T),1), 'NewVariableNames', dimensionName);
            else
                % Split the table rows into multiple tables
                TT =cell(1,nrConditions);
                conditionValues = cell(1,nrConditions);
            end
            cCntr=0;
            for c = fetch(conditions,'value','trials')'
                cCntr= cCntr+1;
                trialsThisCondition = ismember(T.trial,c.trials);
                assert(isscalar(c.value),"Condition value should be a scalar");
                valueThisCondition = c.value{1}; % Stored as cell but is a scalar
                if pv.labelOnly
                    [TT{trialsThisCondition,dimensionName}] = deal(valueThisCondition);
                else
                    TT{cCntr} = T(trialsThisCondition,:);
                    conditionValues{cCntr} =valueThisCondition;
                end
            end
            if pv.labelOnly
                varargout = {TT};
            else
                varargout = {TT,conditionValues};
            end
        end
        function [trials,values] = combine(d,dimension,restriction)
            % Factorial combination of multiple dimensions
            % Returns a cell array with trial numbers that correspond ot
            % the factorial combination of the dimension d with the
            % dimensions whose names are specified as dimNames.
            %
            % EXAMPLE
            %   Create a single row table for the dimension current
            %   dim = ns.Dimension & expt & 'dimension=''current'''
            % Now cross this dimension with the 'stimChannel' dimension
            %   trials= combine(dim,'stimChannel')
            % For an experiment with 5 current levels and 2 stimulation
            % channels this will return a 5x2 cell array where element
            % (i,j) contains the trials in which the ith current was
            % applied to the the jth channel.
            % The function also returns the values of each of the
            % conditions, values(3,2) is a cell array with the value
            % assigned to the third condition in the first dimension and
            % the second condition in the second dimension. (i.e. it
            % matches the layout of the trials output).
            arguments
                d (1,1) {mustHaveRows(d,1)}
            end
            arguments (Repeating)
                dimension (1,1) string
                restriction (1,1) string
            end

            nrDims = 1+numel(dimension);
            dimTrials = cell(1,nrDims);
            dimValue = cell(1,nrDims);
            %  Get the info from the main dimension ; this determines the
            %  experiment
            [dimTrials{1},dimValue{1}]  = fetchn(d*ns.DimensionCondition,'trials','value');
            dTpl = fetch(d); % Only considering the same expt as d
            for i=1:nrDims-1
                dTpl.dimension = char(dimension{i});
                % Fetch each dimension
                if restriction{i} == ""
                    thisRestrict= true;
                else
                    thisRestrict = restriction{i};
                end

                [dimTrials{i+1},dimValue{i+1}]  = fetchn((ns.Dimension&dTpl)*(ns.DimensionCondition & thisRestrict),'trials','value');
            end
            % Check what we have
            nrConditions = cellfun(@numel,dimTrials);
            trials =cell(nrConditions);
            values =cell(nrConditions);
            subs = cell(1,nrDims);
            % Use indexing/subscriptint to combine factorially across an arbitrary number of dimensions
            for i=1:prod(nrConditions)
                [subs{:}] =ind2sub(nrConditions,i);
                thisTrials = cellfun(@(d,r)(dimTrials{d}{r}),num2cell(1:nrDims),subs,'uni',false);
                trials{i} = intersect(thisTrials{:});% Keep trials that occur in each dimension
                values{i} = cellfun(@(d,r)(dimValue{d}{r}{1}),num2cell(1:nrDims),subs,'uni',false);
            end
            out= cellfun(@isempty,trials,'uni',true);
            rowOut = all(out,2:nrDims);
            trials(rowOut,:) = [];
            values(rowOut,:)=[];
        end




    end

    methods (Access=protected)
        function makeTuples(self, key)
            % Called by populate for each (experiment, dimension) pair in
            % keySource not yet present in the table.
            % key contains: subject, session_date, starttime, dimension,
            %               paradigm  (from ns.Experiment * ns.DimensionParmParadigm)

            % Fetch the DimensionParm for this dimension
            parms = fetch1(ns.DimensionParm & key, 'parms');
            % Unpack parms
            plg          = cellstr(parms.plg);
            prm          = cellstr(parms.prm);
            nrPlg = numel(plg);
            nrPrm = numel(prm);
            assert(nrPlg==nrPrm, "Please specify one parameter per plugin")

            if strcmp(parms.fun, '')
                fun = @(x)(x);
            else
                fun = str2func(parms.fun);
            end

            pvSEPARATOR = ":"; % Between parm and value
            ppSEPARATOR = "_"; % Between one parm and the next.

            allTrials = (1:fetch1(ns.Experiment & key, 'nrtrials'))';
            nrTrials  = numel(allTrials);
            valStr    = "";
            valTbl    = table;
            for i = 1:nrPlg
                if parms.nameValueOnly
                    prefix = "";
                else
                    if parms.left==0
                        % Use full name
                        prefix = string(plg{i}) + pvSEPARATOR + string(prm{i}) + pvSEPARATOR;
                    else
                        % Use 1:parms.left of the plg name
                        prefix = string(plg{i}(1:parms.left)) + pvSEPARATOR + string(prm{i}(1:parms.left)) + pvSEPARATOR;
                    end
                end

                ret = get(ns.Experiment & key, plg{i}, 'prm', prm{i}, 'atTrialTime', parms.atTrialTime, 'what', ["data" "trial"]);
                if isempty(ret)
                    fprintf('Plugin %s not in use for %s on %s at %s; skipping dimension %s.\n', plg{i}, key.subject, key.session_date, key.starttime, key.dimension);
                    return;
                end
                prmValues = ret.data;
                if isempty(ret.trial) && isscalar(unique(ret.data))
                    prmTrials = (1:nrTrials)';
                else
                    prmTrials = ret.trial;
                end

                if isempty(prmTrials) && isscalar(prmValues)
                    prmValues = repmat(prmValues, [numel(allTrials) 1]);
                elseif isempty(prmValues) || numel(prmTrials) ~= numel(allTrials) || ~all(prmTrials == allTrials)
                    fprintf('Dimension (%s:%s) not in use for %s on %s at %s\n', plg{i}, prm{i}, key.subject, key.session_date, key.starttime);
                    return;
                end
                prmValues = fun(prmValues);
                if iscell(prmValues)
                    nrPrmValues = cellfun(@numel, prmValues);
                    if all(nrPrmValues == 1)
                        strPrmValues = string(prmValues);
                    else
                        strPrmValues = cellfun(@(x) strjoin(string(x), "_"), prmValues);
                    end
                else
                    strPrmValues = string(prmValues);
                end
                if valStr == ""
                    valStr = prefix + strPrmValues(:);
                else
                    valStr = [valStr prefix + strPrmValues]; %#ok<AGROW>
                end
                valTbl = addvars(valTbl, prmValues, 'NewVariableNames', string(plg{i}) + pvSEPARATOR + string(prm{i}));
            end
            if isempty(prmValues)
                return;
            end
            valStr = fillmissing(valStr, "constant", "unknown");
            isStayTrials = true(size(allTrials));
            if ~isempty(parms.restrict)
                for r = 1:3:numel(parms.restrict)
                    restrictValue = get(ns.Experiment & key, parms.restrict{r}, 'prm', parms.restrict{r+1}, 'atTrialTime', parms.atTrialTime, 'what', 'data');
                    isStayTrials = isStayTrials & ismember(restrictValue(:), parms.restrict{r+2});
                end
            end
            for r = 1:3:numel(parms.exclude)
                excludeValue = get(ns.Experiment & key, parms.exclude{r}, 'prm', parms.exclude{r+1}, 'atTrialTime', parms.atTrialTime, 'what', 'data');
                isStayTrials = isStayTrials & ~ismember(excludeValue, parms.exclude{r+2});
            end
            valStr     = valStr(isStayTrials, :);
            stayTrials = find(isStayTrials);
            valTbl     = valTbl(isStayTrials, :);

            [uValStr, ia, ic] = unique(valStr, 'rows');
            valTbl       = valTbl(ia, :);
            nrConditions = size(uValStr, 1);
            %% Create tuples
            tplC = repmat(key, [nrConditions 1]);
            for c = 1:nrConditions
                tplC(c).name      = strjoin(uValStr(c, :), ppSEPARATOR);
                tplC(c).trials    = stayTrials(ic == c);
                tplC(c).value     = table2cell(valTbl(c, :));
            end
            %% Insert (populate wraps this call in a transaction)
            % Expand trials blob into individual rows for DimensionTrial
            tplT = struct([]);
            for c = 1:nrConditions
                trialNums = num2cell(tplC(c).trials);
                tplRow = repmat(rmfield(tplC(c), {'trials', 'value'}), [numel(trialNums) 1]);
                [tplRow.trial] = deal(trialNums{:});
                tplT = [tplT; tplRow]; %#ok<AGROW>
            end
            insert(self, key);
            tplC = makeMymSafe(tplC);
            insert(ns.DimensionCondition, tplC);
            if ~isempty(tplT)
                tplT = makeMymSafe(tplT);
                insert(ns.DimensionTrial, tplT);
            end
            fprintf('Dimension %s: added %d condition tuples.\n', key.dimension, numel(tplC));
        end
    end
    methods (Static)
        function [nrTrialsPerCondition,overlappingTrials] = summary(expt,pv)
            % Summarize the dimensions in an experiment
            % Returns one table with the number of trials for each
            % condition in each dimension and another table with the number
            % of trials that occur in conditions of separate dimensions
            %
            arguments
                expt (1,1) ns.Experiment {mustHaveRows(expt,1)}
                pv =1
            end


            dims = ns.Dimension & expt;
            assert(exists(dims),"No dimensions have been defined for this experiment (%s on %s @%s)",expt{1}.subject, expt{1}.session_date,expt{1}.starttime);

            dimCons = fetchtable(ns.DimensionCondition & dims,'trials');
            % Pivot tabulating hte number of trials per condition
            % in each dimension
            nrTrialsPerCondition = pivot(dimCons,Columns="dimension",Rows= "name",DataVariable="trials", Method = @(x) numel(x{1}));

            % Get unique trials and dimensions
            allUniqueTrials = unique(cat(1, dimCons.trials{:}));
            allUniqueDims = unique(dimCons.dimension);
            % Create a logical map: Rows = Dimensions, Cols = Trials
            % membership(i, j) is true if Dimension i contains Trial j
            membership = false(numel(allUniqueDims), numel(allUniqueTrials));

            for i = 1:numel(allUniqueDims)
                currentTrials = dimCons.trials(strcmp(dimCons.dimension, allUniqueDims{i}));
                currentTrials = cat(1, currentTrials{:});
                membership(i, :) = ismember(allUniqueTrials, currentTrials);
            end

            % Now your intersection matrix M is just matrix multiplication!
            M = double(membership) * double(membership)';

            % Convert the matrix M to a table
            overlappingTrials= array2table(M, 'VariableNames', allUniqueDims, 'RowNames', allUniqueDims);

        end

        function missing = missingDimensionParms()
            % Return unique (paradigm, dimension) pairs present in ns/dimension
            % but not yet in ns/DimensionParm.  Use this to identify which
            % DimensionParm rows must be created before populate() can run.
            %
            % Returns a struct with fields 'paradigm' and 'dimension'
            % (each a cell array of strings).
            %
            % EXAMPLE
            %   m = ns.Dimension.missingDimensionParms()
            %   % m.paradigm and m.dimension list what needs to be inserted
            schema  = ns.getSchema();
            db = schema.dbname;
            missing = dj.conn().query(sprintf([ ...
                'SELECT DISTINCT d.paradigm, d.dimension ' ...
                'FROM `%s`.`ns/dimension` d ' ...
                'LEFT JOIN `%s`.`ns/#dimension_parm` dp ' ...
                '  ON d.dimension = dp.dimension AND d.paradigm = dp.paradigm ' ...
                'WHERE dp.dimension IS NULL ' ...
                'ORDER BY d.paradigm, d.dimension'], db, db));
            missing = struct2table(missing);
        end

        function tpl = defineParms(expt, plg, prm, name, pv)
            % newer (May 26) use of DimensionParm table.
            %
            % Pass the same arguments previously passed to ns.Dimension.define.
            % then add the return value to the ns.DimensionParm table.
            %
            % See also ns.Dimension.define, ns.DimensionParm
            arguments
                expt (1,1) ns.Experiment
                plg (1,:) {mustBeNonzeroLengthText}
                prm (1,:) {mustBeNonzeroLengthText}
                name (1,1) string
                pv.fun (1,1) function_handle = @(x)(x)
                pv.restrict (1,:) cell = {}
                pv.exclude (1,:) cell = {}
                pv.description (1,1) string = ""
                pv.left (1,1) double = 0
                pv.nameValueOnly (1,1) = false
                pv.atTrialTime (1,1) = 0
                pv.useTable (1,1) = false
            end
            assert(~pv.useTable,"useTable is no longer supported in the DimensionParm design of dimensions.");
            % Derive unique paradigms from the supplied experiment table
            paradigm = string(unique(fetchn(proj(expt, 'paradigm'), 'paradigm')));

            % Convert function handle to string; empty string means identity (default)
            if strcmp(func2str(pv.fun), func2str(@(x)(x)))
                funStr = "";
            else
                funStr = string(func2str(pv.fun));
            end

            parms = struct( ...
                'plg',          {string(plg)}, ...
                'prm',          {string(prm)}, ...
                'fun',          funStr, ...
                'restrict',     {pv.restrict}, ...
                'exclude',      {pv.exclude}, ...
                'left',         pv.left, ...
                'nameValueOnly', pv.nameValueOnly, ...
                'atTrialTime',  pv.atTrialTime);

            tpl = struct( ...
                'dimension',   char(name), ...
                'paradigm',    cellstr(paradigm), ...
                'parms',       parms, ...
                'description', char(pv.description));

        end
    end
end