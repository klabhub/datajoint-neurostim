function final_summary = get_summary_se(y, sub_id, pv)
% GET_SUMMARY_SE Calculates summary statistics with optional Cousineau-Morey correction.
%
%   This function computes descriptive statistics (mean, median, n, se, ci)
%   for various experimental designs. It can handle between-subject,
%   within-subject, and mixed designs. For designs with within-subject
%   factors, it applies the Cousineau-Morey correction to the standard
%   error to better represent within-subject variability.
%
%   DOES NOT WORK WITH CONTINUOUS FACTORS!
%
%   Args:
%       y (vector): The dependent variable scores (numeric).
%       sub_id (vector): The subject/participant identifiers for each score.
%       pv (struct): A struct with optional name-value pairs:
%           'within' - A struct where each field is a within-subject factor vector.
%                      The field names are used as factor names.
%                      Example: struct('Time', time_vec, 'Condition', cond_vec)
%           'between' - A struct for between-subject factors, similar to 'within'.
%
%   Returns:
%       table: A table with summary statistics for each condition.

    arguments
        y (:,1) {mustBeNumeric}
        sub_id (:,1)
        pv.within struct = struct()
        pv.between struct = struct()
        pv.filter_extreme_outliers = false
        pv.outlier_threshold_robust_z = 4;
    end

    if pv.filter_extreme_outliers

        z_scores = do.robust_z(y);
        isValid = abs(z_scores) < pv.outlier_threshold_robust_z;

    else

        isValid = true(size(y));

    end
    % 1. Assemble the main data table from inputs
    data = table(y, sub_id, 'VariableNames', {'score', 'subject'});
    
    
    [data, within_vars] = add_factors_to_table(data, pv.within);
    [data, between_vars] = add_factors_to_table(data, pv.between);
    
    % Handle potential missing data by removing rows with any NaN scores.
    data = data(~isnan(data.score) & isValid, :);

    % 2. Split data by between-subject factors, if they exist
    if isempty(between_vars)
        % If no between-subject factors, process the entire dataset as a single group
        final_summary = process_group(data, within_vars);
    else
        % If between-subject factors exist, process each group separately
        [group_idx, group_info] = findgroups(data(:, between_vars));
        num_groups = height(group_info);
        summary_list = cell(num_groups, 1);

        for i = 1:num_groups
            subset_data = data(group_idx == i, :);
            group_summary = process_group(subset_data, within_vars);
            
            % Prepend the between-subject factor information to the result
            between_info_rep = repmat(group_info(i,:), height(group_summary), 1);
            summary_list{i} = [between_info_rep, group_summary];
        end
        final_summary = vertcat(summary_list{:});
    end
end

% --- Local Helper Functions ---

function summary = process_group(data, within_vars)
% PROCESS_GROUP Handles normalization and stats computation for a single group.
    arguments
        data table
        within_vars (1,:) string
    end
    
    % 1. Normalize data if within-subject factors are present
    if ~isempty(within_vars)
        [data, k] = normalize_within_(data, within_vars);
        % Apply Morey's correction factor for bias
        morey_correction = sqrt(k / (k - 1));
    else
        % No within-factors, so no normalization or correction is needed
        data.normalized_score = data.score;
        morey_correction = 1;
    end
    
    % 2. Compute the final statistics for this group
    grouping_vars = within_vars;
    n_subjects = numel(unique(data.subject));
    summary = compute_stats_(data, grouping_vars, n_subjects, morey_correction);
end


function [data, k] = normalize_within_(data, within_vars)
% NORMALIZE_WITHIN_ Applies Cousineau's normalization to the data.
    arguments
        data table
        within_vars (1,:) string
    end

    % Determine k, the total number of within-subject conditions
    k = height(unique(data(:, within_vars)));
    if k < 2
        warning('Cousineau-Morey correction is intended for 2 or more within-subject conditions. SE will not be corrected.');
        k=1; % to avoid division by zero
    end
    
    % Calculate grand mean and participant means *within the current group*
    grand_mean = mean(data.score);
    participant_means = groupsummary(data, 'subject', 'mean', 'score');
    participant_means.Properties.VariableNames{'mean_score'} = 'participant_mean';
    
    data = outerjoin(data, participant_means(:, {'subject', 'participant_mean'}), 'Keys', 'subject', 'MergeKeys', true);

    % Apply normalization formula
    data.normalized_score = data.score - data.participant_mean + grand_mean;
end


function stats = compute_stats_(data, grouping_vars, n, correction_factor)
% COMPUTE_STATS_ Calculates final summary statistics table.
    arguments
        data table
        grouping_vars (1,:) string
        n double
        correction_factor double
    end

    % Use groupsummary to calculate basic stats
    if isempty(grouping_vars)
        % Handle case with no grouping variables (e.g., between-subjects only)
        stats_raw = groupsummary(data, {}, {'mean', 'median', 'std'}, {'score', 'score', 'normalized_score'});
        stats_raw = removevars(stats_raw, 'GroupCount');
    else
        stats_raw = groupsummary(data, grouping_vars, {'mean', 'median', 'std'}, {'score', 'score', 'normalized_score'});
    end

    % Rename for clarity
    stats_raw.Properties.VariableNames{'mean_score'} = 'mean';
    stats_raw.Properties.VariableNames{'median_score'} = 'median';
    stats_raw.Properties.VariableNames{'std_normalized_score'} = 'sd_norm';

    % Calculate corrected standard error
    se_norm = stats_raw.sd_norm / sqrt(n);
    stats.se = se_norm * correction_factor;
    
    % Calculate 95% Confidence Interval
    df = n - 1;
    if df > 0
        t_crit = tinv(0.975, df);
        ci_margin = t_crit * stats.se;
        stats.ci_lower = stats_raw.mean - ci_margin;
        stats.ci_upper = stats_raw.mean + ci_margin;
    else
        stats.ci_lower = nan;
        stats.ci_upper = nan;
    end
    
    % Combine results into the final table
    stats = [stats_raw(:, setdiff(stats_raw.Properties.VariableNames, {'sd_norm', 'GroupCount'}, 'stable')), ...
             table(repmat(n, height(stats_raw), 1), 'VariableNames', {'n'}), ...
             struct2table(stats)];
end


function [data, var_names] = add_factors_to_table(data, factors_struct)
% ADD_FACTORS_TO_TABLE Helper to add factor columns from a struct to the main data table.
    arguments
        data table
        factors_struct struct
    end
    
    var_names = string(fieldnames(factors_struct))'; % Get field names as a row string array

    if ~isempty(var_names)
        temp_table = table();
        for i = 1:numel(var_names)
            factor_name = var_names(i);
            factor_data = cat(1,factors_struct.(factor_name));

            if size(factor_data, 1) ~= height(data)
                error('Factor "%s" has %d rows, but y has %d rows.', ...
                      factor_name, size(factor_data, 1), height(data));
            end
            temp_table.(factor_name) = factor_data;
        end
        data = [data, temp_table];
    end
end

