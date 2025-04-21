%{ 
# Functional Connectivity method parameters
fctag : varchar(32)
ctag : varchar(32)  
---
description : varchar (1024) # Description of the method
parms : blob                 # Struct with the parameters. 
%}
% The parms struct can have the following members:
% .method =  the name of a function  that computes the functional
% connectivity. This function takes the parms struct and an ns.CChannel
% object as its input and returns [fc,p,err,src,trg]. For example, see
% fc.pearson. Note that this method will only be called with the CChannels
% in its skeleton.
%
% See fc.Skeleton for the definition of these parameters:
% .skeleton.method  - Which function to call to define a skeleton
% .skeleton.paradigm - Which paradigms to include when defining a skeleton
%
% EXAMPLE
% tpl  =struct('parms',struct('method','PEARSON'),
%               'fctag','pearson',
%               'description','Pearson correlation based on the entire timecourse');
% insert(fc.Parm,tpl)
% See fc.Fc/makeTuples how this is used.
%
classdef Parm <dj.Lookup

end