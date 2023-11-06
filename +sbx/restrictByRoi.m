function r = restrictByRoi(tbl,varargin)
% Function to restrict a ns.CChannel by properties of an ROI 
% This is needed because the ROI in sbx.PreprocessedRoi are not directly
% linked (through a foreign key in the schema) to the channels of the ns.CChannel 
% even though the code ensures that channel==roi. 
% INPUT
% tbl  - a ns.CChannel table
% varargin - any restrictions that can be applied to sbx.PreprocessedRoi
%           e.g. 'pcell>0.75', or 'compact<1'
% OUTPUT
% tbl - a ns.CChannel table restricted to the channels that contain
% information on ROIs that match the restrictions.
% 
% BK - Nov 2023
arguments
    tbl (1,1) ns.CChannel   
end
arguments (Repeating)
    varargin  % Anything that can be a restriction (struct, string, char, other dj.GeneralRelvar
end
sessions  = ns.Session & tbl;
allTpls=  [];
% Because roi are defined per session, loop over sessions
for sess = fetch(sessions)'
    % Extract the relevant ROI
    thisTbl = sbx.PreprocessedRoi &sess;
    thisTbl.restrict(varargin{:}); % Keep only rois that match the restrictions.
    keep= fetch(thisTbl,'roi');  % Fetch roi
    thisRois = num2cell([keep.roi]);
    thisTpls = fetch((tbl & sess) & struct('channel',thisRois'));  % Restrict by assuming/knowing that roi=channel.
    if isempty(allTpls)
        allTpls =thisTpls;
    else
        allTpls = cat(1,allTpls,thisTpls);
    end
end

r =proj(tbl & allTpls); % Restrict the original table with the list of tpls collected above.