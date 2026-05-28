function EEG = addEvents(EEG,plg,prm)
% Add events from PluginParameter to an EEG dataset. This can be used for
% visual inspection or comparison to ICA components.
% The events will be added to the EEG.event structure with the type field
% identifying the prm name.
%
% See Also ephys.eeglab.dataset
arguments
    EEG (1,1) struct % EEGLAB dataset  created by ephys.eeglab.dataset
    plg (1,1)  string % One plugin's name   (e.g. edf)
    prm (1,:) string % Event names. A strgin array of event names (e.g. ["startblink" "startsacc"])
end

assert(isfield(EEG.etc,'neurostim'),'This EEG struct does not contain the necessary neurostim information.')
plgRel  =  ns.Plugin & EEG.etc.neurostim.expt  & struct('plugin_name',plg);
assert(exists(plgRel),"The %s plugin could not be found in this experiment (%s@%sT%s)",plg,EEG.etc.neurostim.expt.subject,EEG.etc.neurostim.expt.session_date,EEG.etc.neurostim.expt.starttime);
prmRel = ns.PluginParameter & plgRel & in('property_name',prm);
if count(prmRel)==0
    fprintf(2,'No %s:%s events found.\n',plg,strjoin(prm,'/'));
else
    % Map neurostim time to egi time.
    inverseClockParms= [1 -EEG.etc.neurostim.clockParms(2)]./EEG.etc.neurostim.clockParms(1);
    
    type  = fetchn(prmRel,'property_name');
    nrEventTypes = numel(type);
    egiLatency = cell(nrEventTypes,1);
    trial = cell(nrEventTypes,1);
    for typeCntr= 1:nrEventTypes
        p = fetch(prmRel & struct('property_name',type{typeCntr}),'property_nstime','property_trial');
        egiTimeSeconds  = polyval(inverseClockParms,p.property_nstime);
        egiLatency{typeCntr} = egiTimeSeconds(:)'*EEG.srate; % Must be in samples
        trial{typeCntr} = p.property_trial(:)';
        type{typeCntr} = p.property_name;
    end
    % Add the events to the EEG dataset
    fprintf('Adding %d %s events from %s plugin to the EEG.event struct.\n',numel([egiLatency{:}]),strjoin(prm,'/'),plg)
    EEG = eeg_addnewevents(EEG, egiLatency, type);
end
