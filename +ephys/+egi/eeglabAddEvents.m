function EEG = eeglabAddEvents(EEG,plg,prm)
arguments
    EEG (1,1) struct % EEGLAB dataset 
    plg (1,1)  string % Plugin name 
    prm (1,:) string % Event names
end

assert(isfield(EEG.etc,'neurostim'),'This EEG struct does not contain the necessary neurostim information.')

plgRel  =  ns.Plugin & EEG.etc.neurostim.expt  & struct('plugin_name',plg);
assert(exists(plgRel),"The %s plugin could not be found in this experiment (%s)",plg,EEG.etc.neurostim.expt)
prmRel = ns.PluginParameter & plgRel & in('property_name',prm);
if count(prmRel)==0
    fprintf(2,'No %s:%s events found.\n',plg,strjoin(prm,'/'));
else
    inverseClockParms= [1 -EEG.etc.neurostim.clockParms(2)]./EEG.etc.neurostim.clockParms(1);
    
    type  = fetchn(prmRel,'property_name');
    nrEventTypes = numel(type);
    egiLatency = cell(nrEventTypes,1);
    trial = cell(nrEventTypes,1);   
    for typeCntr= 1:nrEventTypes 
        p = fetch(prmRel & struct('property_name',type{typeCntr}),'property_nstime','property_trial');        
        egiTimeSeconds  = polyval(inverseClockParms,p.property_nstime);
        egiLatency{typeCntr} = egiTimeSeconds(:)'*EEG.srate;
        trial{typeCntr} = p.property_trial(:)';
        type{typeCntr} = p.property_name;
    end
EEG = eeg_addnewevents(EEG, egiLatency, type);

end
