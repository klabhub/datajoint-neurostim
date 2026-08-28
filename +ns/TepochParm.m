%{
# TepochParm: parameters to transform epoch data.
ttag : varchar(32)  # tag for this transformed epoch
etag : varchar(32) # which epochs to transform 
---
fun : longblob # Function to do the transform (see ns.cache/compute)
window  = NULL      : tinyblob  # Start and stop time of the epoch. Defaults to
                            #   entire epoch
channels            : blob  # Channels to include/exclude. Defaults to all in etag
                            #   list of channels or a function to select
                            #   channels for a given key, or one of the
                            #   following
                            #       1. (default) {'good'}  % includes channels without
                            #           any flags
                            #       2. {flag1, flag2,...} cell (1,:) 
                            #           includes only the channels with
                            #           the specified flags which must 
                            #           correspond to the flags 
                            #           in ns.C/info.etc.noiseDetection.reference.badChannels                            
                            #       3. {'good', {flag1, flag2,...}} % combines 1&2
                            #       4. {'exclude', {flag1, flag2,...}} excludes
                            #           only the specified flag names                            
                            #       5. {..., 'interpolate', interp_args} interpolates all the
                            #           bad channels where interp_args is a
                            #           a cell array containing the inputs
                            #           to `ns.interp_idw()`
                            #       6. {{'interpolate', {flag1, flag2,...}, interp_args}
                            #           interpolates only specific channels
                            #       7. [channel1,...] where channel1, ...2,
                            #            etc. are the list of channels to include
                            #       8. {func_str, {arg1, arg2,...}}
                            #           func_str is a function to select
                            #           channels. Must be entered by itself
                            #           or as the final entry after 1-7 above.
                            #           USAGE:
                            #               func = str2func(func_str);
                            #               ecTbl_select = func(...
                            #                   ns.EpochChannel & key, ...
                            #                   arg1, arg2,...)
                            #       9. {'good', 'exclude', 'interpolated', {flag1,...}}
                            #           excludes channels interpolated
                            #           during PrepPipeline. If
                            #           followed by {flag1, ...}, excludes only the
                            #           interpolated channels with the specified flags.
trials              : blob  # Defaults to all in etag
                            #   a function string to select trials for a 
                            #   given key, or any of the entries to
                            #   `channels` specified above
                            #   [1:4, 8] above
                            #   Caveats:
                            #       2. flags must correspond to the flags 
                            #           in ns.Epoch/plg.categories or
                            #           ns.Epoch/art.categories
conditions = NULL   : blob  # Defaults to all, conditions to
                            #   include/exclude
average = NULL      : blob  # factors to average over. string or cell array
                            #    with one ore more:
                            #    'trials', 'channels', 'conditions', func_str
                            # If func_str given, the function must
                            #    take the same inputs (self,key) as
                            #    ns.Tepoch.makeTuples(). rest of the
                            #    args must be given in the trailing element in the cell array.
                            #       {func_str, {arg1,arg2,...}}
                            #    Outputs:
                            #        group_assignment: int group ids 
                            #            assigned to trials
                            #        group_ids (optional): list of
                            #            unique groups
%}

classdef TepochParm < dj.Lookup & dj.DJInstance
    
    methods
        function insert(tbl,tpl)
            arguments
                tbl (1,1)
                tpl (:,1) struct
            end

            tpl = namedargs2cell(tpl);
            tpl = set_defaults(tpl{:});
            
            insert@dj.Lookup(tbl,tpl);
        end
    end

end

function tpl = set_defaults(tpl)

arguments
    
    tpl.ttag char
    tpl.etag char
    tpl.fun struct
    tpl.window (1,2) double
    tpl.channels = {'good'}
    tpl.trials cell = {'good'}
    tpl.conditions cell
    tpl.average cell

end

tpl = makeMymSafe(tpl);
end