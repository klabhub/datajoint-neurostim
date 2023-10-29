function [v,name] = trialsPerCondition(condition)
% Function to extract trials per condition based on a few ways to specify a
% condition : an ns.Condition table, a list of trials, or a cell array with lists of trials. 
% 
% INPUT
% condition
% OUTPUT
% v -  Cell array of trials in each condition
% name - Names of the conditions (retrieved from ns.Condition table, empty
% otherwise).
% 
name = cell(1,numel(condition));
[name{:}]  =deal('');
     %% Determine trial->condition mapping
            if isempty(condition)
                v = {[]};% All trials
            elseif isa(condition,'ns.Condition')
                % Collect trials from ns.Condition table
                v = cell(1,count(condition));
                cCntr = 0;
                for c=fetch(condition)'
                    cCntr = cCntr+1;
                    v{cCntr} = [fetch(ns.ConditionTrial & c,'trial').trial];
                    name{cCntr} = c.name;
                end
            elseif isnumeric(condition)
                v = {condition};
            elseif iscell(condition)
                v = condition;
            else
                error('Condition must be empty, an ns.Condition table, a list of trials, or a cell array with lists of trials. Not a %s',class(condition))
            end

end