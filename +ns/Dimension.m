%{
# A dimension refers to a group of conditions (each of which represente a sset of trials in an Experiment).
-> ns.Experiment 
dimension        : varchar(32)              # Name for the dimension (e.g,'ori')
---
plugin           : blob             # The plugin(s) that define this dimension
parameter        : blob             # The parameter(s) that define this dimension
description = NULL :varchar(1024)   #   A brief description
%}
%
% A 'Dimension' is a set of conditions. Each condition represents a set of
% trials in an experiment and is represented in the Part table
% DimensionCondition.
%
% For instance, when mapping an orientation tuning curve with the gabor plugin,
% the dimension could be called 'ori', and the trials in which the orientation
% was 30 degrees are stored in a DimensionCondition tuple with the name
% gabor:orientation:30.
%
% The definition of what constitutes a dimension (i.e. the set of
% stimulus parameters) varies per paradigm, hence the Dimension table has to be populated
% by the user (i.e., it is dj.Manual). The ns.Dimension.define function provides the
% tool to do this.
%
% For instance, in the example above:
% ns.Dimension.defin(expt, 'gabor','orientation','ori')
% to add a row to the Dimension table and a set of rows tothe DimensionCondition table.
%
%
% The .name field creates unique names for the conditions. For instance, if one of
% the orientations was 10, its name would become gabor:orientation:10
%
% The ns.Tuning function uses the Dimension table to determine tuning
% curves.
%
% See also ns.Tuning
classdef Dimension < dj.Manual
    methods (Access=public)
        function [trials,values] = combine(d,dimNames)
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
                dimNames (1,1) string
            end
   
            nrDims = 1+numel(dimNames);
            dimTrials = cell(1,nrDims);
            dimValue = cell(1,nrDims);
            %  Get the info from the main dimension ; this determines the
            %  experiment
            [dimTrials{1},dimValue{1}]  = fetchn(d*ns.DimensionCondition,'trials','value');            
            dTpl = fetch(d); % Only considering the same expt as d
            for i=1:nrDims-1
                dTpl.dimension = char(dimNames{i});
                % Fetch each dimension 
                [dimTrials{i+1},dimValue{i+1}]  = fetchn((ns.Dimension&dTpl)*ns.DimensionCondition,'trials','value');
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
        end
    end
    methods (Static)
        function define(expt,plg,prm,name,pv)
            arguments
                expt (1,1) ns.Experiment
                plg (1,:) {mustBeNonzeroLengthText}
                prm (1,:) {mustBeNonzeroLengthText}
                name (1,1) string
                pv.fun (1,1) function_handle = @(x)(x) % Function handle to convert plg.prm into a number that defines the condition (default is to do nothing; the value determines the condition)
                pv.restrict (1,:) cell = {} % Define only in this subset of trials by specifying a set of allowed values for a plugin parameter
                pv.description (1,1) string ="" % A description to add to the table
                pv.replace (1,1) logical = false % Set to true to replace (all) existing conditions from this expt and this dimension.
                pv.left (1,1) double = NaN  % Reduce the names to this number of chars from the left
                pv.nameValueOnly (1,1) = false  % Set to true to define condition names based o the prm values alone (and not their name).

            end
            if ~exists(expt)
                fprintf('Empty experiment table; no dimensions addded.\n')
                return;
            end
            if ischar(plg);plg={plg};end
            if ischar(prm);prm={prm};end
            pvSEPARATOR = ":"; % Between parm and value
            ppSEPARATOR = "_"; % Between one parm and the next.

            nrPlg = numel(plg);
            nrPrm = numel(prm);
            assert(nrPlg==nrPrm,"Please specify one parameter per plugin")
            existing = (ns.Dimension & 'dimension=''' + name + '''')  & expt ;
            if pv.replace
                del(existing)
            end
            % Only process experiments in which this dimension has not already
            % been defined.
            expt = expt - proj(existing);
            exptTpl = fetch(expt);
            nrExpt = numel(exptTpl);
            tic
            fprintf('Defining Dimension %s for %d experiments...',name,nrExpt)
            totalNrTpl=0;
            % Collect values as string (to make a unique condition name)
            % and as actual values (to store as values).
            for e =1:nrExpt
                valStr = "";
                valTbl = table;
                for i=1:nrPlg
                    if pv.nameValueOnly
                        prefix="";
                    else
                        if isnan(pv.left)
                            prefix = string(plg{i})+pvSEPARATOR +string(prm{i})+pvSEPARATOR;
                        else
                            % Reduce prefix
                            prefix = string(plg{i}(1:pv.left)) +pvSEPARATOR +string(prm{i}(1:pv.left)) +pvSEPARATOR;
                        end
                    end
                    prmValues = get(ns.Experiment & exptTpl(e),plg{i},'prm',prm{i},'atTrialTime',0)';
                    if isempty(prmValues)
                        % This experiment did not use the plugin; error in
                        % the condition specification, skip to the next
                        % experiment
                        break;
                    end
                    % Transform the values using the user-specified
                    % function
                    prmValues = pv.fun(prmValues);
                    if iscell(prmValues)
                        nrPrmValues= cellfun(@numel,prmValues);
                        if all(nrPrmValues==1)
                            strPrmValues= string(prmValues);
                        else
                            strPrmValues= cellfun(@(x) strjoin(string(x),"_"),prmValues);
                        end
                    else
                        strPrmValues =string(prmValues);
                    end
                    if valStr==""
                        valStr= prefix+strPrmValues(:); 
                    else
                        valStr = [valStr prefix+strPrmValues];%#ok<AGROW>
                    end
                    valTbl = addvars(valTbl,prmValues,'NewVariableNames',plg{i} + pvSEPARATOR + prm{i});
                end
                if isempty(prmValues)
                    % This experiment did not use the plugin; error in
                    % the condition specification, skip to the next
                    % experiment
                    continue;
                end
                valStr = fillmissing(valStr,"constant","unknown");
                if ~isempty(pv.restrict)
                    % Determine which trials meet the specified condition
                    restrictValue = get(ns.Experiment & exptTpl(e),pv.restrict{1},'prm',pv.restrict{2}, 'atTrialTime',0)';
                    stay = ismember(restrictValue,pv.restrict{3});
                    valStr =valStr(stay,:);
                    stayTrials = find(stay);
                    valTbl= valTbl(stay,:);
                else
                    stayTrials = 1:fetch1(ns.Experiment&exptTpl(e),'nrtrials');
                end
                % Find the unique combinations of parameters
                [uValStr,ia,ic] =unique(valStr,'rows');                
                valTbl = valTbl(ia,:);
                nrConditions = size(uValStr,1);
                %% Create tuples
                tplD = mergestruct(exptTpl(e), ...
                    struct('dimension',name, ...
                    'description',pv.description,...
                    'plugin',{plg}, ...
                    'parameter',{prm}));
                tplC = repmat(exptTpl(e),[nrConditions 1]);
                for c=1:nrConditions
                    tplC(c).name = strjoin(uValStr(c,:),ppSEPARATOR);
                    tplC(c).dimension  = name;
                    tplC(c).trials = stayTrials(ic==c);
                    tplC(c).value = table2cell(valTbl(c,:)); % Must store as cell for mym
                end
                %% Insert the condition tuples for this experiment as one transaction
                C =dj.conn;
                C.startTransaction
                try
                    insert(ns.Dimension,tplD);
                    insert(ns.DimensionCondition,tplC);
                catch me
                    C.cancelTransaction;
                    rethrow(me)
                end
                C.commitTransaction;
                totalNrTpl = totalNrTpl + numel(tplC);
            end
            fprintf('Done in %s s. Added %d condition tuples.\n',seconds(toc),totalNrTpl);
        end
    end
end