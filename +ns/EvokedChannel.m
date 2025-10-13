%{
# ns.EvokedChannel is part table of ns.Evoked
-> ns.Evoked
channel: int
---
signal: blob
%}

classdef EvokedChannel < dj.Part & dj.DJInstance
    properties (SetAccess = protected)
            master = ns.Evoked
    end

     methods

        function plot(self, pv)

            % arguments
            % end

        end

        function summarize(self)

            arguments
                
                self      
                                
            end

            % fetch data
            dat = fetchtable(self, '*');
            % collapse across channels
            [gru_idx, gru_id] = findgroups(dat(:,{'subject', 'starttime', 'session_date', 'name', 'group'}));
            % Downsample to smallest sampling rate
            dat = [dat,rowfun(@(x) numel(x{:}), dat, InputVariables='signal', OutputVariableNames='length')];
            % Find max and min sample length
            min_len = min(dat.length);            
            % if different sampling rates, downsample to min_len
            if ~isscalar(unique(dat.length))

                dat.signal = arrayfun(@(signal, len) decimate(signal{:}, len/min_len), dat.signal, dat.length, 'UniformOutput', false);

            end

            gru_id.signal = splitapply(@(s) mean(cat(1, s{:}),1), dat.signal, gru_idx);
            % for now, it calculates within-subject error at every
            % timepoint individually
            ep_tbl = fetch(ns.EpochParm & self,'epoch_win');
            t = linspace(ep_tbl.epoch_win(1), ep_tbl.epoch_win(2), min_len);
            within_factors = struct(session=join([gru_id.session_date,gru_id.starttime],'/'), name=gru_id.name, group=gru_id.group);
            summ_tbl = [];
            for ii = 1:min_len

                datN = gru_id.signal(:,ii);
                summN = ns.get_summary_se(datN, gru_id.subject, within = within_factors);
                summN = addvars(summN,t(ii)*ones(height(summN),1), Before='mean', NewVariableNames='time');
                summ_tbl = [summ_tbl; summN];
            end

        end

    end
end