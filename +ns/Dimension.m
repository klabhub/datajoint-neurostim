%{
# A dimension refers to a group of conditions (each of which represente a sset of trials in an Experiment).
-> ns.Experiment 
dimension        : varchar(32)              # Name for the dimension (e.g,'ori')
---
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
     methods (Static)
        function define(expt,plg,prm,name,pv)
            arguments
                expt (1,1) ns.Experiment {mustHaveRows}
                plg (1,:) {mustBeNonzeroLengthText}
                prm (1,:) {mustBeNonzeroLengthText}
                name (1,1) string 
                
                pv.restrict (1,:) cell = {} % Define only in this subset of trials by specifying a set of allowed values for a plugin parameter
                pv.description (1,1) string ="" % A description to add to the table
                pv.replace (1,1) logical = false % Set to true to replace (all) existing conditions from this expt and this dimension.           
                pv.left (1,1) double = NaN  % Reduce the names to this number of chars from the left
                pv.nameValueOnly (1,1) = false  % Set to true to define condition names based o the prm values alone (and not their name).
                
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
            for e =1:nrExpt
                val = [];
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
                    val= [val prefix+string(prmValues)]; %#ok<AGROW> 
                end
                if isempty(prmValues)
                        % This experiment did not use the plugin; error in
                        % the condition specification, skip to the next
                        % experiment                        
                        continue;
                end
                val = fillmissing(val,"constant","unknown");
                if ~isempty(pv.restrict)
                    % Determine which trials meet the specified condition
                      restrictValue = get(ns.Experiment & exptTpl(e),pv.restrict{1},'prm',pv.restrict{2}, 'atTrialTime',0)';
                      stay = ismember(restrictValue,pv.restrict{3});
                      val =val(stay,:);
                      stayTrials = find(stay);
                else 
                      stayTrials = 1:fetch1(ns.Experiment&exptTpl(e),'nrtrials');
                end
                [uVal,~,ix] = unique(val,"rows"); % Sort ascending.
                nrConditions = size(uVal,1);                
                %% Create tuples 
                tplD = mergestruct(exptTpl(e),struct('dimension',name,'description',pv.description));                
                tplC = repmat(exptTpl(e),[nrConditions 1]);
                for c=1:nrConditions                    
                    tplC(c).name = strjoin(uVal(c,:),ppSEPARATOR);
                    tplC(c).dimension  = name;                    
                    tplC(c).trials = stayTrials(ix==c);                                                                   
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