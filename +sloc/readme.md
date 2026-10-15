# Source Localization

This package uses Fieldtrip to perform source localization on data
in a datajoint-neurostim database. 

The computations are all Fieldtrip based, install that first and run
```ft_defaults``` before starting. The interactive steps (e.g., fiducial markin)
are also Fieldtrip-based; consult the Fieldtrip tutorials for more info.

Although it may be used more generally, the workflow is based on our setup 
where we record EEG with a EGI hydrocel 256 channel net, acquire photogrammetry
in each recording session, and "solve" these models using EGI tools. 


## Example Workflow

1. Identify the nasion, left and right auricular fiducial points for a
selection of subjects. Here for subjects that participated in a specific experiment.


```sloc.Fiducial(qry = ns.Subject & (ns.Experiment & 'paradigm ="rvsBlocked"'))```

1. Define your head and source model by defining a headmodel parm. 
    The parms fields specify the cfg's that are passed to the fieldtrip functions.
    Fields that are not defined will get default values (see ```sloc.HeadmodelParm```).

    ```matlab 
    tpl = [];
    tpl.hmtag = 'brainnetome';  % A unique tag for this headmodel/sourcemodel combo
    tpl.description = 'MNI head model with 5 mm template and 2mm mesh resolution, mapped to the Brainnetome atlas';
    tpl.parms.sourcemodel.cfg.method = 'basedonmni';
    tpl.parms.sourcemodel.template  = 'standard_sourcemodel3d5mm';
    tpl.parms.sourcemodel.atlas   = 'brainnetome/BNA_MPM_thr25_1.25mm.nii';
    tpl.parms.elec = 'GSN-HydroCel-256.sfp';
    tpl.parms.mesh.cfg.downsample = 2;
    insertIfNew(sloc.HeadmodelParm, tpl);
    ```

    You can visualize the models using ```plot(sloc.Headmodel,what="segmentation")```
    or other options for ```what``` such as ```"mesh", "allmeshes", "sourcemodel"```.

1. Let Fieldrtip compute the models: ```populate(sloc.Headmodel)```

1. Define the Leadfield parameters. The deafults (defined in ```sloc.LeadfieldParm```) are 
   probably ok. The auto mode computes the leadfield most efficiently and can use a parpool if it
   has been defined with nsParPool (and the environment variable NS_PARPOOL).
    ```matlab
        tpl = [];
        tpl.lftag= 'default';
        tpl.description ='Leadfield computation using all default settings';
        tpl.parms = struct('mode','auto');  % Validation will fill in FT defaults
    ```

    You will also have to list for which paradigms a leadfield should be computed.

    ```
    tpl.name  = 'rvsBlocked'; % Name of the paradigm
    tpl.lftag = 'default';  % Name of the leadfield computation parameters
    insertIfNew(sloc.LeadfieldParmParadigm, tpl); % Link the two.
    ```

    Then ``` populate(sloc.Leadfield) ``` will start computing leadfields for each session with 
    an experiment that has the rvsBlocked paradigm. The code will check for a solved GPSC file in that
    folder (the final result of the photogrammetry process).  We call those files ```subject.solved.gpsr```.
    If your convention differs you will have to adjust the code in sloc.Leadfield/makeTuples a bit. 
  

    The sloc.Leadfield model has a plot function that visualizes the norm of the leadfield
    at every dipole location. 
    
1.  Localize neural activity to the sources.
    The sloc.SourceParm table stores the parameters for this. Here is an example that
    determines sources for epochs defined with etag="flicker". This is an SSVEP paradigm
    with a 500 ms baseline and 2000 ms SSVEP response. The ssvep mode determines a common spatial
    filter based on the base and the first response window and then applies this filter to all 
    response windows. The power stored in the sloc.Source table is the ratio of the response windows
    to the base window.


    ```matlab
        tpl = [];
        tpl.stag = 'ssvep';
        tpl.description = 'Source localized SSVEPs comparing post/pre';
        tpl.parms.epoch = "flicker";
        tpl.parms.baseWindow = [-500 0];
        tpl.parms.responseWindow = {[500 1000],[1000 1500],[1500 2000]};
        tpl.parms.mode = 'SSVEP';
        insertIfNew(sloc.SourceParm,tpl);
    ```

    For source SNR, use `mode = 'SNR'`. The response power in the signal band
    around each condition frequency is divided by the mean power in the two
    flanking noise bands. `signalHalfWidth` and `noiseHalfWidth` specify the
    half-widths of these bands in Hz; the noise half-width must be larger than
    the signal half-width. For example:

    ```matlab
        tpl.parms.mode = 'SNR';
        tpl.parms.signalHalfWidth = 0.5;
        tpl.parms.noiseHalfWidth = 2;
        insertIfNew(sloc.SourceParm,tpl);
    ```
    Note that for a signalHalfWidth of 0.5 you'd want a frequency resolution of at least 1 Hz. So at least
    1 s of data. 
    
## ToDo
. Source localization without a subject MRI (template brain)
. Source localization withjout photogrammetry. 

