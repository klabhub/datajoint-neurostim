%{ 
# Functional Connectivity method parameters
fctag : varchar(50)
---
description : varchar (1024) # Description of the method
parms : blob                 # Struct with the parameters. 
%}
% The parms struct should have a memnber called method
% this can be a function_handle, or a string ("PEARSON")
%
% EXAMPLE
% tpl  =struct('parms',struct('method','PEARSON'),
%               'fctag','pearson',
%               'description','Pearson correlation based on the entire timecourse');
% insert(ns.FcParm,tpl)
% See ns.Fc/makeTuples how this is used.
%
classdef FcParm <dj.Lookup

end