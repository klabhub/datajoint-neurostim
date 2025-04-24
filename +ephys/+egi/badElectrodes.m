function [perExpt,perChannel] = badElectrodes(C,parms)
% Read file with a list of bad electrodes
% parms.filename = string that must be in the filename (typically  badElectrodes)
% parms.extension = extension of the file (typically .xlsx)
%
% For instance, for subject 123 there will be a file called
% 123.badElectrodes.xlsx in the session folder. It contains a list of
% electrodes considered "bad".
% See Also ns.Artifact
arguments
    C (1,1) {mustBeCOrExperiment(C), mustHaveRows(C,1)}
    parms (1,1) struct
end

%% Find the relevant file
if isa(C, "ns.C")
    expt =(ns.Experiment & C);
else
    expt = C;
end
file = (ns.File & expt) & struct('extension',parms.extension) & ['filename LIKE "%' char(parms.filename) '%"'];
fldr = folder(expt);
T=table;
if exists(file)
    fn = fetch1(file,'filename');
    ff= fullfile(fldr,fn);
    if exist(ff,'file')
        T = readtable(ff);
        if isempty(T)
            fprintf(2,'File (%s) found but table is empty in  %s\n',fn,fldr)
        end
    else
        fprintf(2,'File (%s) not found in folder %s\n',fn,fldr)
    end
    badElectrodes = T{:,1};
    perChannel =struct('start',[],'stop',[],'trial',[],'channel',num2cell(badElectrodes));
else
    fprintf('No %s for %s.\n',parms.filename,fldr)
    perChannel =struct('start',[],'stop',[],'trial',[],'channel',[]);
end

perExpt = struct('start',[],'stop',[],'trial',[]); % Nothing applicable to all channels

end

%% checker function for arguments block
function mustBeCOrExperiment(inputArg)
    % Checks if the input is either ns.C or ns.Experiment
    if ~(isa(inputArg, 'ns.C') || isa(inputArg, 'ns.Experiment'))
        eid = 'Validation:InvalidClass';
        msg = sprintf('Input must be of type ns.C or ns.Experiment. Got %s.', class(inputArg));
        error(eid, msg);
    end
end