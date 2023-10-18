%{
# A condition refers to a set of trials in an Experiment with matching parameters.
-> ns.Experiment 
name        : varchar(255)      # Condition name (plg_prm_value)
---
conditionnr : smallint unsigned  # Condition number 1:nrConditions
group = NULL : varchar(10)        # Group of conditions 
%}
% As the definition of what constitutes a condition (i.e. the set of
% stimulus parameters) varies per paradigm, this table has to be populated 
% by the user (i.e., it is dj.Manual). The defineConditions function provides the
% tool to do this. For instance, in an experiment where the 'frequency'
% parameter of the 'gabor' plugin was the only parameter that varied across
% trials (i.e., the frequency defines the condition), call:
% defineConditions(ns.Condition, expt, 'gabor','frequency') 
% to create a table of conditions and matching trials (in the
% ns.ConditionTrial part table). 
% The .name field creates unique names for the conditions. For instance, if one of
% the frequencies was 10, its name would become gabor_frequency_10
% 
% The group can be used to define multiple groupings of conditions. For
% isntance, one grouping could assign trials to condition based on the
% visual stimulus presented on the screen (group ="orientation"), while another grouping 
% could assign based on some intervention (before drug /after drug;
% group="drug"). That would allow one analysis to pool over drug state, but
% distinguish among orientations, while another could analysis drug states
% separately.
%
% For an example how to use this table, see sbx.Tuning, which determines 
% tuning for a set of conditions. 
% BK - March 2023.
classdef Condition < dj.Manual
    methods (Access=public)
        function tbl = defineConditions(tbl,expt,plg,prm,pv)
            arguments
                tbl (1,1) ns.Condition
                expt (1,1) ns.Experiment {mustHaveRows}
                plg (1,:) {mustBeNonzeroLengthText}
                prm (1,:) {mustBeNonzeroLengthText}
                pv.left (1,1) double = NaN  % Reduce the names to this number of chars from the left
                pv.replace (1,1) logical = false; % Set to true to replace (all) existing conditions
                pv.group (1,1) string = ""
                pv.nameValueOnly (1,1) = false; % Set to true to define condition names based o the prm values alone (and not their name).
            end
            if ischar(plg);plg={plg};end
            if ischar(prm);prm={prm};end
            pvSEPARATOR = ":"; % Between parm and value
            ppSEPARATOR = "_"; % Between one parm and the next. 

            nrPlg = numel(plg);
            nrPrm = numel(prm);
            assert(nrPlg==nrPrm,"Please specify one parameter per plugin")

            existingConditions = ns.Condition &expt;
            if pv.replace
                del(existingConditions)
            end

            exptTpl = fetch(expt);
            nrExpt = numel(exptTpl);
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
                [uVal,~,ix] = unique(val,"rows"); % Sort ascending.
                nrConditions = size(uVal,1);

                %% Check which conditions already exist and insert new ones 
                
                conditionsThisExpt = ns.Condition & exptTpl(e);
                for c=1:nrConditions
                    tpl = exptTpl(e);                            
                    tpl.name = strjoin(uVal(c,:),ppSEPARATOR);                    
                    if count(conditionsThisExpt & tpl) >0
                        % ALready exists, skip                        
                    else
                        tpl.group = pv.group;      
                        conditionNr = count(conditionsThisExpt) +1;
                        tpl.conditionnr= conditionNr;                    
                        insert(tbl,tpl);
                        %% Create and insert ConditionTrial tuples
                        trials = find(ix==c);
                        nrTrials = numel(trials);
                        tpl =ns.stripToPrimary(ns.Condition,tpl);
                        trialTpl = repmat(tpl,[nrTrials 1]);
                        trials = num2cell(trials);
                        [trialTpl.trial] = deal(trials{:});
                        insert(ns.ConditionTrial,trialTpl)    
                    end
                end                
                   
            end
        end
    end
end