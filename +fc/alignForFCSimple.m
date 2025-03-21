function [D,channelInMatrix] = alignForFCSimple(parms,channels)
% Prepare data to input to any FC method
% There are different ways to create a dataset depending on the analysis
% For ITI in Dirmap, use parms.align.kind = 'trial_conc'
% Can use the whole experiment for Resting, etc.
% TODO: complete the description

    arguments
        parms (1,1) struct
        channels (1,1)  ns.CChannel
    end
    
    % data preparation: make a function for this
    warning('off', 'DataJoint:longCondition')
    
    switch parms.align.kind
    
        case 'trial_conc'
            [t,~,channelInMatrix] = align(ns.C & channels ,'channel',proj(channels), start = parms.align.start, stop = parms.align.stop);
       
            % produce a timepoints x trials x rois 3D matrix
            D = timetableToDouble(t);
            % remove first and last trial to avoid NaNs problems
            D = D(:,2:end-1,:);
            % standardize each trial before concatenating to avoid possible
            % Simpson's paradox when mixing data, and also to guarantee
            % same scale of the data.
            D = zscore(D,0,1);    
            % concatenate all trials into one 2D (datapoints x rois) 
            [timepoints, trials, rois_count] = size(D);
            D = reshape(D, timepoints * trials, rois_count);
            % sort the channels and use that order in the dataset
            [channelInMatrix,sorted] = sort(channelInMatrix);
            D = D(:,sorted); 
            
    
        % TODO: data creation for resting-state data
        % check how the data is being recorded since there are no trials.
        
        otherwise
            error('Unknown data creation method. Stop function.')
    end