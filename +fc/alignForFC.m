function [D,channelInMatrix,bad_rois] = alignForFC(parms,channels)
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
            allC = ns.C & channels;
            D = [];
            bad_rois = [];
            % loop through all experiments of the paradigm selected
            for thisC = fetch(allC)'
                % threshold based on parms.skeleton.pcell
                if ~isempty(parms.skeleton.pcell)
                    rois = fetchn(sbx.PreprocessedRoi & thisC & ['pcell>',parms.skeleton.pcell],'roi');
                    disp(numel(rois))
                    [t,~,channelInMatrix] = align(allC & thisC ,channel = rois, start = parms.align.start, stop = parms.align.stop);
                else
                    [t,~,channelInMatrix] = align(allC & thisC ,'channel',proj(channels & thisC), start = parms.align.start, stop = parms.align.stop);
                end
           
                % produce a timepoints x trials x rois 3D matrix
                d = timetableToDouble(t);
                % remove first and last trial to avoid NaNs problems
                d = d(:,2:end-1,:);
                % standardize each trial before concatenating to avoid possible
                % Simpson's paradox when mixing data, and also to guarantee
                % same scale of the data.
                d = zscore(d,0,1);    
                % concatenate all trials into one 2D (datapoints x rois) 
                [timepoints, trials, rois_count] = size(d);
                d = reshape(d, timepoints * trials, rois_count);
                % sort the channels and use that order in the dataset
                [channelInMatrix,sorted] = sort(channelInMatrix);
                d = d(:,sorted); 
                % bad rois checking
                % where should this be saved?
                % check for rois with variance = 0 (all zeros probably)
                zero_var_idx = var(d, 1) == 0;    
                if any(zero_var_idx)
                    % Remove bad ROIs from d
                    d(:, zero_var_idx) = [];
                    % Return actual ROI IDs
                    br = channelInMatrix(zero_var_idx);
                    % remove also from channelInMatrix
                    channelInMatrix(zero_var_idx) = [];
                else
                    br = [];
                end
                bad_rois  = cat(1,bad_rois,br);
                % concatenate experiments datasets across datapoints
                % D is 2D matrix
                D = cat(1,D,d); 
            end
    
    
        % TODO: data creation for resting-state data
        % check how the data is being recorded since there are no trials.
        
        otherwise
            error('Unknown data creation method. Stop function.')
    end

    % turn back on the warning
    warning('on', 'DataJoint:longCondition');
