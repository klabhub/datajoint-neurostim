function [keepFrameIx,frameNsTime] = framesForExperiment(key)
% Use the metadata added to the experiment by sbx.Preprocessed to match
% frames to experiments
% Restrict to experiments that sbx.Preprocessed uses.
analyzeExptThisSession = analyze(ns.Experiment & (ns.Session & key),strict=false);
exptT = ns.getMeta( analyzeExptThisSession,["nrframes" "nrplanes"]);
exptT = sortrows(exptT,"starttime");
exptT = convertvars(exptT,["nrframes" "nrplanes"],"double");
hasNoFrames = isnan(exptT.nrframes) |  exptT.nrframes==0;
exptT = exptT(~hasNoFrames,:);

row = find(key.starttime==exptT.starttime);
cumFrames = cumsum(exptT.nrframes);
if row>1
    start = cumFrames(row-1);
else
    start = 0;
end
keepFrameIx =start + (1:exptT.nrframes(row));
frameNsTime = get(ns.Experiment & key,'mdaq','prm','laserondighigh','what',"clocktime");
nrTTL = numel(frameNsTime);
% There always appears to be 1 extraneous TTL; check the match with
% this assumption
assert(abs(exptT.nrframes(row)-(nrTTL/exptT.nrplanes(row)))<2,'Cannot map SBX frames to trials; TTL-Frame mismatch (%d TTL in mdaq %d frames for %d planes in sbx).\n',nrTTL,exptT.nrframes(row),exptT.nrplanes(row));
frameNsTime(1) =[];
end