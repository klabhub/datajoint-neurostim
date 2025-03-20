%{
# Functional Connectivity 
-> ns.C         # Continuous data used to compute FC.
-> fc.Skeleton  # Parameters that define the FC computation 
---
quality = NULL : float # Some measure of quality of the FC (optional)
%}
%
% Functional connectivity table. 
%  The ns.C foreign key identifies the continuous data for which functional
%  connectivity is computed.
% The fc.Parm is a lookup table that defines how FC is computed. The user
% has full control by specifying a method in the fc.Parm parms. 
%
% The dj.Part table fc.FcPair stores the actual FC values per pair
%
% See Also fc.Parm
classdef Fc < dj.Computed
     methods
        
         function plot(tbl,pv)
            % Rudimentary plot function. Needs to be extended with options
            arguments
                tbl (1,1) {mustHaveRows(tbl)}               
                pv.bothSides (1,1) logical = true; % Show symmetric matrix
                pv.alpha (1,1) double = Inf; % so far for testing                
            end
            tiledlayout('flow')
            glob_max = [];
            glob_min = [];
            for tpl = fetch(tbl)'
                nexttile
                FC = fcMatrix(tbl & tpl,type =["fc" "p"]);
                FC.fc(FC.p > pv.alpha) = NaN;
                if pv.bothSides
                    % get NaN entries in the lower triangle
                    mask = isnan(tril(FC.fc, -1));
                    % copy values from upper to lower
                    aux = FC.fc.';
                    FC.fc(mask) = aux(mask);
                    
                end
                glob_max = cat(1,glob_max,max(max(FC.fc)));
                glob_min = cat(1,glob_min,min(min(FC.fc)));


                imagesc(FC.fc);
                % to make the FC matrix square
                axis equal;        
                pbaspect([1 1 1]); 
                xlabel 'Channel'
                ylabel 'Channel'                
                title(sprintf('%s\n for %s-%s@%s',tpl.fctag,tpl.subject,tpl.session_date,tpl.starttime),'Interpreter', 'none')  
            end
            % use one unique color range for all plots to compare them
            rwb = centeredColormap(256,min(glob_min),max(glob_max));
            axHandles = findall(gcf, 'Type', 'axes');
            for i = 1:numel(axHandles)
                colormap(axHandles(i), rwb);
                colorbar(axHandles(i));
                set(axHandles(i), 'CLim', [min(glob_min) max(glob_max)]);
            end
        end

        function v = fcMatrix(tbl,pv)
            % Convenience function that returns a struct with FC as a
            % matrix, together with the channels that the FC was
            % calculated for and, optionally, the associated p and err
            % values form the FcPair table.
            arguments
                tbl (1,1) fc.Fc
                pv.type (1,:) string = "fc"  % Set to ["fc" "p" "err" ]  to include p and err values.
            end
                columns = cellstr([pv.type "source" "target"]);
                FC = fetch(fc.FcPair & tbl,columns{:});                
                [v.channel] = sort(unique([[FC.source];[FC.target]]));
                [~,srcIx] = ismember([FC.source],v.channel);
                [~,trgIx] = ismember([FC.target],v.channel);                
                nrChannel = numel(v.channel);                
                sz = [nrChannel nrChannel];
                for tp = pv.type
                    v.(tp) = nan(sz);
                    v.(tp)(sub2ind(sz,srcIx,trgIx)) = [FC.(tp)];                        
                end
        end

        % copy to datajoint and remove from here.
        function cmap = centeredColormap(n, dataMin, dataMax)
            % n is the number of colors (optional, default = 256)
            % dataMin and dataMax are the minimum and maximum of the data range
            if nargin < 1
                n = 256; % Default number of colors
            end
            % Define color for negative (blue), center (white), and positive (red)
            colors = [0 0 255/256; 1 1 1; 1 0 0]; % Blue, White, Red
            % Create interpolation points
            x = [dataMin, 0, dataMax]; % Map dataMin to blue, 0 to white, dataMax to red
            % Create the colormap using interpolation for each channel (R, G, B)
            cmap = interp1(x, colors, linspace(dataMin, dataMax, n), 'linear');
        end

    end

    methods (Access=protected)

        function makeTuples(tbl,key)
            % When calling (par)populate, this function is called with the
            % key containing the information on a set of continuous data
            % (i.e. limited to an experiment), and a skeleton (which is 
            % computed per session, so it can be based on multiple
            % experiments). 
            % In principle FC can be computed for any C data that has
            % multiple channels.


            %% Fetch the parms from the FcParms table
            parms = fetch1(fc.Parm & key,'parms');
            assert(~isempty(which(parms.method)),sprintf('Unknown FC method %s',parms.method))
           
           %% Call the user specified method function 
            % The method takes the parms and channels as input and returns 
            % fc, p, and err, for each src and trg. 
            %channels = ns.CChannel & fc.SkeletonChannel & key ; % Restrict channels to the skeleton                        
            % Use the fc.SkeletonFcPair to compute the experiment FC
            % Using channels is not enough if computing multiple
            % regression. Need all the sources for  a target for example.
            
            % now restrict to the paradigm in parms.skeleton.paradigm
            % is this correct, or do we want the fc for everything?
            

            channels = ns.CChannel & fc.SkeletonChannel & key;
            % worth it having a new function or adjust the original
            % alignForFC?
            [D,channelInMatrix] = fc.alignForFCSimple(parms,channels); 
            sk =  fc.SkeletonFcPair & key;
            % for now multiple regression for FC
            [FC,p,err,src,trg] = feval(parms.method,parms,D,channelInMatrix,sk);      

            % Insert the key into the Fc table
            % TODO: key.quality = add some generic quality measure...            
            insert(tbl,key);

            % Prepare tuples for the FcPair table
            tpl = struct('fc',num2cell(FC(:)),'p',num2cell(p(:)),'err',num2cell(err(:)),'source',num2cell(src(:)),'target',num2cell(trg(:)));
            tpl = mergestruct(ns.stripToPrimary(fc.Fc,key),tpl);
            % Insert.
            chunkedInsert(fc.FcPair,tpl);
        end
    end
end
