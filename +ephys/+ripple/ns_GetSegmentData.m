function [ns_RESULT, TimeStamp, Data, SampleCount, UnitID] = ...
    ns_GetSegmentData(hFile, EntityID, Index,NoWaveForms)
% ns_GetSegmentData - Retrieves segment data by index
% Usage:
% 
% [ns_RESULT, TimeStamp, Data, SampleCount, UnitID] = ...
%             ns_GetSegmentData(hFile, EntityID, Index, NoWaveForms)
% Description:
% 
% Returns the Segment data values in entry Index of the entity EntityID 
% from the file referenced by hFile. The data values are returned in Data. 
% The timestamp of the entry id returned in TimeStamp.  The number of 
% samples written to Data is returned in SampleCount.  The data variable 
% should be accessed as a 2-dimensional array for samples and sources.
% 
% Parameters:
% 
% hFile               Handle/Identification number to an open file.
% 
% EntityID            Identification number of the entity in the data file.
% 
% Index               Index number of the requested Segment data item.
%       
% NoWaveForms        Boolean to set whether to read actual waveforms or
%                    time /id only.
% Remarks:
% 
% A zero unit ID is unclassified, then follow unit 1, 2, 3, etc.
%
% Return Values:
% 
% TimeStamp           Time stamp of the requested Segment data item.
% 
% Data                Variable to receive the requested data.
% 
% SampleCount         Number of samples returned in the data variable.
% 
% UnitID              Unit classification code for the Segment Entity.
% 
% ns_RESULT           This function returns ns_OK if the file is 
%                     successfully opened. Otherwise one of the following 
%                     error codes is generated:
% 
% ns_BADFILE          Invalid file handle passed to function
% 
% ns_BADENTITY        Invalid or inappropriate entity identifier specified
% 
% ns_BADINDEX         Invalid entity index specified
% 
% ns_FILEERROR        File access or read error

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     The Wisteria Neuroshare Importer is free software: you can 
%     redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     The Wisteria Neuroshare Importer is distributed in the hope that it 
%     will be useful, but WITHOUT ANY WARRANTY; without even the implied 
%     warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
%     See the GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with the Wisteria Neuroshare Importer.  If not, see 
%     <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin <4
    NoWaveForms = false;
end
TimeStamp    = [];
Data         = [];
SampleCount  = [];
UnitID       = [];

%check hFile
if ~isstruct(hFile)
    ns_RESULT = 'ns_BADFILE';
    return
end

%check Entity
if ~(isnumeric(EntityID))||...
    (uint16(EntityID)~=EntityID)||...
    (EntityID > length(hFile.Entity))||...
    (EntityID < 1)||...
    ~strcmp(hFile.Entity(EntityID).EntityType, 'Segment')
    ns_RESULT = 'ns_BADENTITY';
    return
end

%check Index
if any(Index<1| Index>hFile.Entity(EntityID).Count)
    ns_RESULT = 'ns_BADINDEX';
    return
end

ns_RESULT = 'ns_OK';
nrIndex = numel(Index);

%create fileInfo Structure
fileInfo = hFile.FileInfo(hFile.Entity(EntityID).FileType);
SampleCount = (fileInfo.BytesDataPacket - 8)/2;
PacketIndex = find(fileInfo.MemoryMap.Data.PacketID==...
    hFile.Entity(EntityID).ElectrodeID);

PacketIndex =PacketIndex(Index);
TimeStamp =...
    double(fileInfo.MemoryMap.Data.TimeStamp(PacketIndex))/30000;
UnitID = fileInfo.MemoryMap.Data.Class(PacketIndex);

% skip to event wave data. But only read if requested
 if NoWaveForms
    Data = [];
 else     
     Data=nan(nrIndex,SampleCount,'int16');

     % Read wave form data 
    % These consecutive reads with a fixed skipping of 8 bytes shoulwd work
    % but dont... left for later BK 
%     if nrIndex==1 || all(diff(Index)==1)
%         % Consecutive file reads. Read them at once
%     offset = double(fileInfo.BytesHeaders) + 8 +...
%             fileInfo.BytesDataPacket*(PacketIndex(1)-1);    
%          fseek(fileInfo.FileID, offset, -1);    
%     Data = fread(fileInfo.FileID, nrIndex*SampleCount, '52*int16=>int16',8); 
%     Data = reshape(Data,SampleCount,nrIndex)';
%     else % Non consecutive reads. Have to loop (could chunk over consecutive bits...)          
       
        for i=1:nrIndex
            offset = double(fileInfo.BytesHeaders) + 8 + ...
                double(fileInfo.BytesDataPacket) * double(PacketIndex(i)-1);    
            fseek(fileInfo.FileID, offset, -1);    
            Data(i,:)= fread(fileInfo.FileID, SampleCount*1, 'int16=>double');  
        end
%     end        
    % scale the data
    Data = Data*hFile.Entity(EntityID).Scale;
end