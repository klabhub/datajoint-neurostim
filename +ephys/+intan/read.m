function [header,data] = read(filename,pv)
% Reads Intan Technologies data files
%
% RHS Files uses read_Intan_RHS2000.m Version 3.0, 8 February 2021
%
% Specify channels to read, separately for :
% digIn
% digOut
% amplifier
% dac
% adc
% stim
% By default none of the channels are returned.To return all, specifiy NaN.
arguments
    filename (1,1) string
    pv.headerOnly (1,1) logical =false
    pv.digIn (1,:) double {mustBeNonnegative,mustBeInteger}= [];  
    pv.digOut (1,:) double {mustBeNonnegative,mustBeInteger}= [];  
    pv.amplifier (1,:) double {mustBeNonnegative,mustBeInteger}= [];  
    pv.dac (1,:) double {mustBeNonnegative,mustBeInteger}= [];  
    pv.adc (1,:) double {mustBeNonnegative,mustBeInteger}= [];  
    pv.stim (1,:) double {mustBeNonnegative,mustBeInteger}= [];  

end


if ~exist(filename,"file")
    error('Intan RHS file %s does not exist',filename)
end

[~,~,ext] =fileparts(filename);

switch upper(ext)
    case '.RHS'
        %% Call the Intan Function        
        %#ok<*USENS> % Many variables are created in the script
        ephys.intan.read_Intan_RHS2000_file;
        %% Rearrange to put data in a more usable format.        
        header.frequency = frequency_parameters;
        header.stim  = stim_parameters;
        header.trigger = spike_triggers;
        header.amplifier = amplifier_channels;
        header.digIn = board_dig_in_channels;
        header.digOut = board_dig_out_channels;
        header.adc  = board_adc_channels;
        header.dac = board_dac_channels;
        if ~pv.headerOnly
            if isnan(pv.digIn)
                pv.digIn = 1:numel(header.digIn);
            end
            if isnan(pv.digOut)
                pv.digOut = 1:numel(header.digOut);
            end
            if isnan(pv.amplifier)
                pv.amplifier= 1:numel(header.amplifier);
            end
            if isnan(pv.dac)
                pv.dac = 1:numel(header.dac);
            end
            if isnan(pv.adc)
                pv.adc= 1:numel(header.adc);
            end
            if isnan(pv.stim)
                pv.stim = 1:numel(header.stim);
            end
            % Keep only the requested channels            
            data.digIn = board_dig_in_data(pv.digIn,:)'; % Channels to columns
            data.digOut = board_dig_out_data(pv.digOut,:)'; % Channels to columns
            data.amplifier = amplifier_data(pv.amplifier,:)'; % microvolts
            data.dac = board_dac_data(pv.dac,:)'; % volts
            data.adc = board_adc_data(pv.adc,:)'; % volts
            data.stim = stim_data(pv.stim,:)';  % microamps
            data.time  = t; %time in seconds
        else
            data = struct('digIn',[],'digoOut',[],'amplifier',[],'dac',[],'adc',[],'stim',[],'time',[]);
        end
    case '.RHD'
        error('RHD Not implemented yet.')
    otherwise
        error('Unknown Intan data file extension %s .',ext);
end


end