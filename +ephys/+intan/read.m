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
    %% Specify the channels that should be read. NaN means all.
    pv.digIn (1,:) double = [];  
    pv.digOut (1,:) double = [];  
    pv.amplifier (1,:) double = [];  
    pv.dac (1,:) double = [];  
    pv.adc (1,:) double = [];  
    pv.stim (1,:) double = [];  
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
    case '.RHD'
        error('RHD Not implemented yet.')
    otherwise
        error('Unknown Intan data file extension %s .',ext);
end


end