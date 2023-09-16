%{
# A condition refers to a set of trials in an Experiment with matching parameters.
-> ns.Experiment 
conditionnr : smallint unsigned  # Condition number 1:nrConditions
name        : varchar(255)      # Condition name (plg_prm_value)
%}
% As the definition of what constitutes a condition (i.e. the set of
% stimulus parameters) varies per paradigm, this table has to be populated 
% by the user (i.e., it is dj.Manual). The populate function provides the
% tool to do this. For instance, in an experiment where the 'frequency'
% parameter of the 'gabor' plugin was the only parameter that varied across
% trials (i.e., the frequency defines the condition), call:
% populate(ns.Condition, expt, 'gabor','frequency') 
% to create a table of conditions and matching trials (in the
% ns.ConditionTrial part table). 
% The .name field creates names for the conditions. For instance, if one of
% the frequencies was 10, its name would become gabor_frequency_10
% 
% With another table that contains trial-based information (for instance
% sbx.Activity)  you can then analyze only the trials in condition 1 with:
%  sbx.Activity & (ns.ConditionTrial & struct('condition',1))
% Alternatively you can use the 'name' of the condition in the restriction:
% If one of the conditions had 10 as the frequency, 
%  sbx.Activity & (ns.ConditionTrial & struct('name','gabor_frequency_10'))
%
% BK - March 2023.

classdef Condition < dj.Manual
    methods (Access=public)
        function populate(tbl,expt,plg,prm,pv)
            arguments
                tbl (1,1) ns.Condition
                expt (1,1) ns.Experiment {mustHaveRows}
                plg (1,:) {mustBeNonzeroLengthText}
                prm (1,:) {mustBeNonzeroLengthText}
                pv.left (1,1) double = NaN  % Reduce the names to this number of chars from the left
            end
            if ischar(plg);plg={plg};end
            if ischar(prm);prm={prm};end
            SEPARATOR1 = "-"; % Between parm and value
            SEPARATOR2 = "+"; % Between one parm and the next. 

            nrPlg = numel(plg);
            nrPrm = numel(prm);
            assert(nrPlg==nrPrm,"Please specify one parameter per plugin")

            exptTpl = fetch(expt);
            nrExpt = numel(exptTpl);
            for e =1:nrExpt
                val = [];
                for i=1:nrPlg
                    if isnan(pv.left)
                        prefix = string(plg{i})+SEPARATOR1 +string(prm{i})+SEPARATOR1;
                    else
                        % Reduce prefix
                        prefix = string(plg{i}(1:pv.left)) +SEPARATOR1 +string(prm{i}(1:pv.left)) +SEPARATOR1;
                    end
                    prmValues = get(ns.Experiment & exptTpl(e),plg{i},'prm',prm{i},'atTrialTime',0)';
                    val= [val prefix+string(prmValues)]; %#ok<AGROW> 
                end
                val = fillmissing(val,"constant","unknown");
                [uVal,~,ix] = unique(val,"rows"); % Sort ascending.
                nrConditions = size(uVal,1);

                %% Create and insert Condition tuples
                tpl = repmat(exptTpl(e),[nrConditions 1]);
                for c=1:nrConditions
                    tpl(c).name = strjoin(uVal(c,:),SEPARATOR2);
                    tpl(c).conditionnr= c;                    
                    if count(tbl & tpl(c))==0
                        insert(tbl,tpl(c));
                    end
                end
                
                
                %% Create and insert ConditionTrial tuples
                trialTpl = [];
                for c=1:nrConditions
                    trials = find(ix==c);
                    nrTrials = numel(trials);
                    thisTpl = repmat(tpl(c),[nrTrials 1]);
                    trials = num2cell(trials);
                    [thisTpl.trial] = deal(trials{:});
                    trialTpl = [trialTpl ;thisTpl]; %#ok<AGROW> 
                end
                insert(ns.ConditionTrial,trialTpl)
            end
        end
    end
end