%{
# TepochParm: parameters to transform epoch data.
ttag : varchar(32)  # tag for this transformed epoch
etag : varchar(32) # which epochs to transform 
---
fun : longblob # Function to do the transform (see ns.cache/compute)
window  = NULL  : tinyblob  # Start and stop time of the epoch. Defaults to
                            #   entire epoch
channels = NULL     : blob  # Channels to include/exclude. Defaults to all in etag
                            #   list of channels or a function to select
                            #   channels for a given key, or one of the
                            #   following
                            #       1. {'good'}  % includes channels without
                            #           any flags
                            #       2. flag_names cell (1,:) % includes
                            #           chnnels with certain flag_names
                            #       3. ['good', flag_names] % combines 1&2
                            #       4. ['exclude', flag_names] 
trials = NULL       : blob  # Defaults to all in etag
                            #   a function string to select trials for a 
                            #   given key, or one of the
                            #   following
                            #       1. {'good'}  % includes trials without
                            #           any flags
                            #       2. flag_names cell (1,:) % includes
                            #           trials with certain flag_names
                            #       3. ['good', flag_names] % combines 1&2
                            #       4. ['exclude', flag_names] 
conditions = NULL   : blob  # Defaults to all, conditions to
                            #   include/exclude
average = NULL      : blob  # factors to average over. string or cell array
                            #    with one ore more:
                            #    'trials', 'channels', 'conditions', func_str
                            # to restrict averaging to certain
                            #    conditions, trials, and channels
                            #    provide the values or function string 
                            #    in the cell array
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

            % TODO define insert to check valid values
            %               assert(isempty(parms.average) || all(ismember(parms.average,["trial" "channel"])),"Tepoch averaging must ...")

            tpl = makeMymSafe(tpl);
            insert@dj.Lookup(tbl,tpl);
        end
    end

end