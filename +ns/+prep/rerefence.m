function signal = rereference(signal, params)

% Rereferences signal

arguments

    signal double
    params.channel_locations = []
    params.method {ismember(params.method, {'average', 'laplacian','channel'})} = 'average'
    params.nan_handle {ismember(params.nan_handle, ["omit", "interpolate"])} = 'omit'
    params.interpolation_method {ismember(params.interpolation_method, ["inverse_distance", "spherical_spline"])} = 'inverse_distance'
    params.interpolation_params = {}
    params.ref_channel

end

ch_dim = ndims(signal)-1;
n_ch = size(signal, ch_dim);
org_signal = signal;
if strcmp(params.nan_handle, "interpolate")

    fprintf('Interpolating nan channels.\n')
    isInterp = all(isnan(signal), setdiff(1:ndims(signal),ch_dim));

    goodChanIdx = find(~isInterp);
    badChanIdx = find(isInterp);

    [ref_signal, interp_opts] = do_interp(goodChanIdx, badChanIdx);
    params.interpolation_params = [params.interpolation_params, 'idwPower', interp_opts.idwPower];

    fprintf('\t Interpolation complete.\n')
else
    ref_signal = signal;
end

fprintf('Re-referencing through %s\n', upper(params.method));

switch params.method
    

    case 'average'

        signal = org_signal - mean(ref_signal, ch_dim, 'omitnan');

    case 'laplacian'

        signal = ref_signal;
        
        for iCh = 1:n_ch

            predChIdx = setdiff(1:n_ch, iCh);

            interp_signalN = do_interp(predChIdx, iCh);

            if ismatrix(interp_signalN) 

                ref_signal(iCh,:) = interp_signalN(iCh,:);

            else

                ref_signal(:,iCh,:) = interp_signalN(:,iCh,:);

            end            

        end
        signal = org_signal - ref_signal;

    case 'channel'

        signal = signal - signal(params.ref_channel,:);

end

fprintf("\t Rereferencing complete.\n");


    function [interp_signal, interp_options] = do_interp(pred_ch_idx, tar_ch_idx)

        switch params.interpolation_method

            case 'inverse_distance'

                [interp_signal, interp_options] = ns.interpolate_by_inverse_distance(signal, params.channel_locations, pred_ch_idx, tar_ch_idx, params.interpolation_params{:});

            case 'spherical_spline'

                % to be developed
        end

    end

end