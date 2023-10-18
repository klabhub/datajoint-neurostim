function [inSet,notInSet] = resampleTrials(condition,withReplacement,frac)
% Resample trials to use in bootstrapping. This resampling
% makes sure to resample trials from each of the unique
% stimulus conditions so that the resampled trials have all of
% the conditions.
%
% Input
% conditions - A vector of numbers identifying the conditions for each
%               trial/
% withReplacement - set to true to sample with replacement
% (used by bootstrapping).
% frac  - The fraction of trials to resample. 1 means resample
%           all, 0.5 with replacement=false means split halves.
% OUTPUT
% inSet - The set of selected trials  (index into condition)
% outSet - The trials not in the set.
arguments    
    condition (1,:) double
    withReplacement (1,1) logical = false
    frac (1,1) double {mustBeInRange(frac,0,1)} =  1
end
cCntr= 1;
uCondition = unique(condition);
trialPerCond = cell(1,numel(uCondition));
for u= uCondition
    trialPerCond{cCntr} = find(condition==u);
    cCntr = cCntr+1;
end
if withReplacement
    % Sampling with replacement
    inSet = cellfun(@(x) (x(randi(numel(x),[1 ceil(frac*numel(x))]))),trialPerCond,'uni',false);
else
    % Sampling without replacement (e.g. to split 80/20 or
    % split halves)
    inSet = cellfun(@(x) (x(randperm(numel(x),ceil(frac*numel(x))))),trialPerCond,'uni',false);
end
notInSet = cellfun(@(x,y) setxor(x,y),trialPerCond,inSet,'uni',false);
inSet = cat(2,inSet{:});
notInSet =cat(2,notInSet{:});
end
