%{
# ns.EvokedChannel is part table of ns.Evoked
-> ns.EvokedChannel
channel: int
---
signal: blob
%}

classdef EvokedChannel < dj.Part & dj.DJInstance
    properties (SetAccess = protected)
            master = ns.Evoked
    end
end