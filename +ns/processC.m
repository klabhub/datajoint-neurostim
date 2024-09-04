function [signal,time,channelInfo,recordingInfo] = processC(key,parms)
% Normally a C-table is filled from data in afile, but sometimes you want to 
% further process data that are already in ns.C. This function handles that
% option. The user addas a CParm entry with ns.processC as the 'fun' so
% that this file gets called, but the parms structure must also contain a
% .fun entry to specify which (user-defined) function should perform the
% actual processing on a row of ns.C.   This 'fun' must meet the same
% requirements as a regular ns.C preprocessing function (see ns.C and an
% example in sbx.read).  The ns.CParm parms structure also
% contains further (user-defined) instructions for processing. 

assert(isfield(parms,'fun'), 'To process an entry in the ns.C table, the ns.CParm parms struct should specify which function (''fun'') to use');

[signal,time,channelInfo,recordingInfo] = feval(parms.fun,key,parms);