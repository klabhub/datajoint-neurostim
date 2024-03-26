function pick = pickOne(mffInfo,key,which)
% Pick one relevant file from a struct array of candidates.
arguments 
    mffInfo (:,1) struct
    key (1,1) struct
    which (1,1) string {mustBeMember(which,["BEFORE" "AFTER" "NEAREST"])} = "NEAREST"    
end

if numel(mffInfo)==1
    % A single file found - we'll use that
    pick = mffInfo;
else
        % More than one candidate. This assumes the clocks are not wildly
        % out of sync on the Neurostim and EGI machine...        
        mffTime  = cellfun(@(x) datetime(x,"InputFormat","HHmmss","Format","HHmmss"),{mffInfo.time});            
        neurostimTime = datetime(key.starttime,"InputFormat","HH:mm:ss","Format","HHmmss");
        switch upper(which)
            case "BEFORE"
                ix =  find(mffTime<neurostimTime & (neurostimTime-mffTime)>0,1);    %find the most recent one just before the neurostim file
            case "AFTER"
                ix =  find(mffTime>neurostimTime & (mffTime-neurostimTime)>0,1);    %find the most recent one just after the neurostim file
            case "NEAREST"
                [~,ix] =  min(abs(mffTime-neurostimTime));            
        end      
        pick =mffInfo(ix);
end
end