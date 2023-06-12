%{
#  Preprocessing parameter sets used for ephys data
-> ephys.Array    
prep:  varchar(128)
---
fun : varchar(255)    # The function that reads the data and deso the preprocessing. See ephys.Preprocessed
description : varchar(1024)  # Short description
parms : longblob  # struct containing all parameters
%}
%
% EXAMPLE
% The user fills this table with named parameter sets of their chosing that
% they then use to populate the Preprocessed table.
% For instance, if ephys.ripple.prep is a function that takes an ns.Experiment row
% and a struct with parameters as its input, and returns the MUAE [nrSamples nrChannels]
% and the channel numbers [nrChannels 1]
% muaeParms =struc(....) % Whatever is needed by ephys.ripple.prep
% Then add this PrepParm to the table:
% insert(ephys.PreParm, struct('prep','muae', ...
%                                 'fun','ephys.ripple.prep', ...
%                                 'description', 'Standard Super and Roelfsema MUAE', ...
%                                 'parms',muae))
%
% Then populate(ephys.Preprocessd,'prep=''muae'') will populate the
% Preprocessed table with MUAE using these settings
%
% BK - June 2023
classdef PrepParm < dj.Lookup

end