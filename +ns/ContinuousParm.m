%{
#  Preprocessing parameter sets used for continuous data
tag:  varchar(128)
---
fun : varchar(255)    # The function that reads the data and does the preprocessing. 
extension: varchar(10) # Extension of the data files to which this applies. 
description : varchar(1024)  # Short description
parms : longblob  # struct containing all parameters
%}
%
% EXAMPLE
% The user fills this table with named parameter sets of their chosing that
% they then use to populate the Continuous table.
%
% BK - June 2023
classdef ContinuousParm < dj.Lookup

end