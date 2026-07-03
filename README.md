# DataJoint schema and tools for Neurostim

This folder contains scripts to start a Datajoint schema and database analysis pipeline for
data collected with Neurostim (<https://github.com/klabhub/neurostim>).

To read more about DataJoint, how to setup a DataJoint SQL server, see <https://datajoint.com>.

## Installation

After completing the DataJoint installation, make sure that you can connect to your server by typing
`dj.conn` on the Matlab command prompt. This should return some information on the connection to your server(e.g., the host computer running the MySQL server and the user name). You should only proceed to the following steps once `dj.conn` works.

1. Clone this reposiory to a folder on your machine (`'c:\github\datajoint-neurostim'`). 
    This includes a submodule that is a fork of the datajoint-matlab toolbox that has additional features. 
1. Add it to your matlab path:
```matlab 
    fldr = 'c:\github\'
    subFldrs = fullfile(fldr,{'datajoint-neurostim' 'datajoint-neurostim/datajoint-matlab' , ['datajoint-neurostim/tools/mym/' char(mexext)]});
    addpath(subFldrs{:})
    ```
1. To try if it worked, go to a temp folder and initialize a database:

```matlab
dj.conn
cd c:\temp;
root = 'the root folder containing your neurostim data files. This folder has years as subfolders.';
 nsInitializeDataJoint(pwd,'mytest','ns')
```
That will create a database called 'mytest' on the MySQL server and give you access to the basic tables for Neurostim files (Subject, Session, Experiment,Plugin,etc.) in the `+ns` package. 
### Scanbox
The `+sbx` package sets up additional tables for two-photon imaging data collected with the ScanBox software (Neurolabware.com). 
The preprocessing uses various Python packges (E.g. [suite2p] (https://github.com/MouseLand/suite2p#readme)and [Cascade](https://github.com/PTRRupprecht/CascadeTorch) package )
These Python packages are run  "OutOfProcess" to avoid conflicts with Matlab libraries. 

1. Follow the instructions from suite2p or cascade to install the Python packages in their own conda/mamba environment. 
1. Use an environment variable to specify 


## Usage

Once DataJoint is installed here are the functions that setup a pipeline:

1. `nsInitializeDataJoint`  (create a schema/database in the SQL database)
1. `nsScan`                 (scan a folder with Neurostim output files and add them to the database)  

`nsScan` can be called from the command line, or run interactively from the `nsMeta` Matlab app.

Here is a complete example:

```matlab
% Define the root folder 
root = 'Z:\data';
setenv('NS_ROOT',root)  % May need to set again when Matlab restarts
% Create a database on the Datajoint server
codeFolder =  'c:\temp\project\code';
nsInitializeDataJoint(codeFolder,'ProjectName','ns')
% Scan all files in the year 2023 (underneath the root folder)
nsScan('date','01-Jan-2023','schedule','y','root',root,...
                    'readFileContents',true,'addToDataJoint',true,'safeMode',false);
% The database tables have now been created and filled:
ns.Experiment % Shows the Experiment table in the ProjectName project
% The get() function of the ns.Experiment class gives access to all Plugin properties in the experiment.
% For instance in an experiment where the orientation of the grating stimulus changed in each trial, the following call retrieves the orientation for each trial. 
ori = get(ns.Experiment ,'gabor','prm','orientation','atTrialTime',0)
```


BK -  June 2026.
