function T = ripplePrep(key)
% Read data from Ripple data files, preprocess them according to the parameters 
% passed in the parms struct.  This handles MUAE, LFP, SPIKES, and EEG
% data and relies on the Ripple neuroshare tools to read the nve, nsx etc
% files.
arguments 
    key % The keysource of the Preprocessed table (Experiment and PrepParm tuple)   
end 
T = [];

parms = fetch(ephys.PrepParm & key,'*');
files  =  fetch(ns.File & key);

% Fetch which files to read

% Open them with neuroshare. 
[entityInfo,fileInfo,hFile]= ephys.ripple.open(filename);

if isempty(e)
    % Get all electrocdes
end

% Check that the requested signals exist

% Read the signals 

% Preprocess

switch upper(parms.type)
    case 'MUAE'

    case 'EEG'
    case 'LFP'
    case 'SPIKES'
end


end


function [entityInfo,fileInfo,hFile]= open(filename)

            if ~exist(filename,'file')
                error('ripple:open','File %s does not exist',filename)
            else
                [errCode, hFile] =ns_OpenFile(filename);
            end
            
            if ~strcmpi(errCode,'ns_OK');siberr(['ns_OpenFile failed with ' errCode ]);end            
            [errCode, fileInfo] = ns_GetFileInfo(hFile);
            if ~strcmpi(errCode,'ns_OK');siberr(['ns_GetFileInfo failed with ' errCode ]);end            
            % Initialize the entity info and then read it in.
            entityInfo = repmat(struct('EntityLabel', '','EntityType',0,'ItemCount',0),[fileInfo.EntityCount 1]);
            
            for i = 1:sFileInfo.EntityCount
                [~, entityInfo(i,1)] = ns_GetEntityInfo(hFile, i);               
            end
            % Add the ChannelID information
            for i = 1:fileInfo.EntityCount
                entityInfo(i,1).ChannelID = hFile.Entity(i).ElectrodeID;
                entityInfo(i,1).Scale = hFile.Entity(i).Scale;
            end
end
