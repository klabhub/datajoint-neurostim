%{
# Epoching parameters to segment and preprocess trial data. 
etag : varchar (32) # unique tag
---
ctag            : varchar(32)       # C data that can be epoched with these parms
dimension       : varchar(32)       # Condition from the dimension table
window          : tinyblob          # Start and stop time of the epoch.
channels =NULL  : blob              # Channels to include. Defaults to all in the ctag
align           : blob              # struct defining the align event
prep            : blob              # struct containing the preprocessing parameters
art             : blob              # struct array containing the artifact removal parameters
%}

% MOz Feb, 2025
classdef EpochParm < dj.Lookup & dj.DJInstance
    methods
        function insert(self, tuples, varargin)
           % Overload the insert method to do argument validation and setup
           % defaults
            pv = namedargs2cell(tuples);
            tuples = self.validate(pv{:});
            % Call the superclass function after validation
            insert@dj.Lookup(self, tuples, varargin{:})
        end
    end

    methods (Static, Access = protected)
        function pv = validate(pv)
            arguments
                pv.etag (1,1) string
                pv.ctag  (1,1) string
                pv.dimension (1,1) string
                pv.window (1,2) 
                pv.channels (1,:) {mustBeNumeric} = []
                pv.prep (1,1) struct {prep.mustBePrepParm}  = struct('dummy',true);
                pv.art  (1,1) {prep.mustBeArtParm} = struct('dummy',true);
                pv.align (1,1) struct =struct('dummy',true);
            end  

             % validate 'dimension' and 'plugin' exist in ns.Dimension
            dimTbl = ns.Dimension & struct('dimension',pv.dimension);
            assert(count(dimTbl), ...
                'Dimension table does not contain dimension value of "%s"', pv.dimension);
            cTbl = ns.C & struct('ctag',pv.ctag);
            assert(count(cTbl), ...
                'C table does not contain ctag value of "%s"', pv.ctag);

            if isfield(pv.align,'dummy')
                % Default to the startTime of the plugin that defined the
                % dimension. 
                %  Check that there is only one plugin for this dimension
                G = proj(dimTbl, 'dimension');                         
                multipleOptions = aggr(G, dimTbl, 'count(distinct plugin)->n') & 'n>1';
                assert(count(multipleOptions)==0,"The %s dimension links to multiple plugins. Cannot pick a default align.",pv.dimension);
                % Setup the align struct
                plg = fetch1(dimTbl,'plugin','LIMIT 1');
                pv.align = struct('plugin',plg{1},'event','startTime');
            end
                          
        end
    end

end



