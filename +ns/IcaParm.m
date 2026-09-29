%{
# Parameter sets for storing EEGLAB ICA results
itag         :  varchar(32)      # A unique name for this ICA result definition
ctag         :  varchar(32)      # The ctag of the ns.C data to which this applies
session      : tinyint           # If 0 (default): one ICA per experiment. If 1: one ICA per session (see ns.IcaSession)
---
description  :  varchar(1024)    # Short description
parms        :  longblob         # struct with ICA-specific options/metadata
%}
%
% When session=false (default), ns.Ica is populated independently per
% experiment. When session=true, ns.IcaSession first concatenates all
% experiments in the session and runs ICA once; ns.Ica then copies the
% resulting weights into one row per experiment so all downstream tables
% (ns.Label, ns.Epoch) remain unchanged.
%
% The parms struct is passed to ephys.eeglab.preprocess via parms.eeglab.ica.
% An optional parms.filt struct (with fields locutoff/hicutoff) applies a
% high-pass filter before ICA.
%
% EXAMPLE: per-experiment ICA (default)
%   tpl = struct('itag', 'runica_default', ...
%                'ctag', 'eeg', ...
%                'session', false, ...
%                'description', 'runica with default settings', ...
%                'parms', struct('extended', 1));
%   insert(ns.IcaParm, tpl)
%
% EXAMPLE: session-level ICA with pre-filter
%   tpl = struct('itag', 'runica_session', ...
%                'ctag', 'eeg', ...
%                'session', true, ...
%                'description', 'runica across full session, high-passed at 1 Hz', ...
%                'parms', struct('extended', 1, ...
%                                'filt', struct('locutoff', 1, 'hicutoff', 0)));
%   insert(ns.IcaParm, tpl)
%
% See also ns.Ica, ns.IcaSession
%
% BK - May 2026
classdef IcaParm < dj.Lookup

    methods
        function insert(self, tuples, varargin)
            
            % Validate then insert into IcaParm.
            for i = 1:numel(tuples)
                pv = namedargs2cell(tuples(i).parms);
                tuples(i).parms = ns.IcaParm.validate(pv{:});
            end
            insert@dj.Lookup(self, makeMymSafe(tuples), varargin{:});
        end
    end

    methods (Static, Access = protected)
        function tpl = validate(tpl)
            % Validate and supply defaults for a single IcaParm tuple before insertion.
            arguments
                tpl.filt (1,1) struct = struct([])  % Pre-ICA filtering (high-pass recommended)
                % inputs of pop_runica /runica see help runica for details
                tpl.pca (1,1) double = inf          % Number of components (Inf = all)                
                tpl.icatype (1,1) string  {mustBeMember(tpl.icatype,["runica","binica","jader", "fastica"])}= "runica"
                tpl.extended (1,1) double = 1
                tpl.sphering (1,1) string {mustBeMember(tpl.sphering,["on" "off"])}= "on"
                tpl.maxsteps (1,1) double = 512
                tpl.stop (1,1) double = 1e-7
            end
                        
            % Optional filt sub-struct
            if isfield(tpl, 'filt')
                f = tpl.filt;
                assert(isstruct(f), '"parms.filt" must be a struct.');
                assert(isfield(f, 'locutoff') && isnumeric(f.locutoff) && isscalar(f.locutoff), ...
                    '"parms.filt.locutoff" must be a numeric scalar (Hz; use 0 to skip).');
                assert(isfield(f, 'hicutoff') && isnumeric(f.hicutoff) && isscalar(f.hicutoff), ...
                    '"parms.filt.hicutoff" must be a numeric scalar (Hz; use 0 to skip).');
                assert(f.locutoff >= 0 && f.hicutoff >= 0, ...
                    '"parms.filt" cutoff frequencies must be non-negative.');
                assert(f.locutoff > 0 || f.hicutoff > 0, ...
                    '"parms.filt" must specify at least one non-zero cutoff frequency.');
            end

            % Warn if ctag not yet in ns.C (may not exist at definition time, so just warn)
            if isfield(tpl, 'ctag') && count(ns.CParm & struct('ctag', tpl.ctag)) == 0
                warning('ns:IcaParm:unknownCtag', ...
                    'ctag "%s" has no matching rows in ns.CParm yet.', tpl.ctag);
            end
        end
    end

end
