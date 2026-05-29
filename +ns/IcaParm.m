%{
# Parameter sets for storing EEGLAB ICA results
itag         :  varchar(32)      # A unique name for this ICA result definition
ctag         :  varchar(32)      # The ctag of the ns.C data to which this applies
---
description  :  varchar(1024)    # Short description
parms        :  longblob         # struct with ICA-specific options/metadata
%}
%
% EXAMPLE
% struct('itag','runica_default', ...
%        'ctag','eeg', ...
%        'description','Store EEGLAB pop_runica decomposition from EEG preprocessing')
%
% See also ns.Ica
%
% BK - May 2026
classdef IcaParm < dj.Lookup

end
