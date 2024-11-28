%% Various QC queries

%% Plugin parameter type/storage
% Check consistency of property_type across experiments (e.g., a plugin
% parameter that is sometimes stored as a bytestream and in other
% experiments as a parameter. Such differences cause problems in the user
% code. 
    T = fetchtable(ns.PluginParameter,'property_name','plugin_name','property_type');
    % Check same type across experiments
    G = groupsummary(T,["plugin_name" "property_name" ],"numunique","property_type");
    mltpl = G.numunique_property_type>1;
    if any(mltpl)
        M = innerjoin(T,G(mltpl,:),"Keys",["plugin_name" "property_name"]);        
        [~,ix] = unique(M(:,["plugin_name" "property_name" "property_type"]));
        M = M(ix,:)
    end
    %% Load cic for the affected experiments for inspection
    expts = unique(M(:,["starttime" "session_date" "subject"]));
    clear c
    for i=1:height(expts)
        c(i)= open(ns.Experiment & struct(table2struct(expts(i,:))));
    end

    %% Show the difference    
    clc
        for rw=1:height(M) 
            thisC = string({c.subject}) == M.subject(rw) & string({c.startTimeStr})==M.starttime(rw);
            fprintf(2,"%s:%s -> %s\n",c(thisC).file,M.property_name(rw),M.property_type(rw));
            prm = c(thisC).(M.plugin_name(rw)).prms.(M.property_name(rw));
            fprintf(1,"%s %d \n",class(prm.log{1}),numel(prm.log))                        
        end  
    
        %% Fix DJ
        % After fixing the code, can bring the pluginparameter storage in DJ up to date by running
        % This will re-load all experiments (without deleting them first).
        % If you want to ensure referential integrity, set pedantic=true,
        % but this will delete all derived tables (C, tuning,etc.)
        % Set safemode to 0 to avoid the confirmation of the deletion of
        % the pluginparameters for every file
        % dj.config('safemode',0)
        % updateWithFileContents(ns.Experiment,newOnly=false)
        % Set it back
        % dj.config('safemode',1)
    