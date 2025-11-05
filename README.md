# DataJoint schema and tools for Neurostim

This folder contains scripts to start a Datajoint schema and database analysis pipeline for
data collected with Neurostim (<https://github.com/klabhub/neurostim>).

To read more about DataJoint, how to setup a DataJoint SQL server, see <https://datajoint.com>.

## Installation

After completing the DataJoint installation, make sure that you can connect to your server by typing
`dj.conn` on the Matlab command prompt. This should return some information on the connection to your server(e.g., the host computer running the MySQL server and the user name). You should only proceed to the following steps once `dj.conn` works.

1. Clone this reposiory to a folder on your machine (`'c:\github\datajoint-neurostim'`)
1. Add it to your matlab path (`addpath('c:\github\datajoint-neurostim'`)
1. This should work with the regular datajoint toolbox from datajoint.com, but the repo also contains a fork of this toolbox that may have some additional features. You can use this fork instead by adding the datajoint-matlab submodule to your path. ((`addpath('c:\github\datajoint-neurostim\datajoint-matlab'`)
1. To try if it worked, go to a temp folder and initialize a database:

```matlab
dj.conn
cd c:\temp;
root = 'the root folder containing your neurostim data files. This folder has years as subfolders.';
 nsInitializeDataJoint(pwd,'mytest',root,'ns')
```
That will create a database called 'mytest' on the MySQL server and give you access to the basic tables for Neurostim files (Subject, Session, Experiment,Plugin,etc.) in the `+ns` package. 
### Scanbox
The `+sbx` package sets up additional tables for two-photon imaging data collected with the ScanBox software (Neurolabware.com). The preprocessing uses the [suite2p] (https://github.com/MouseLand/suite2p#readme), which requires Python 3.8 to run.

1. Install miniconda
1. Follow the instructions from the suite2p readme to install suite2p in its own conda environment named 'suite2p'

There are two ways in which Matlab can call the Python code. The first ('InProcess') is by letting Matlab know where to find the python you just installed.

For instance, I did this:

```matlab
pyenv('ExecutionMode','InProcess','Version','c:\Users\bartk\miniconda3\envs\suite2p\python.exe');
```

To check that this works. Load the default suite2p options. This should return a Python dict in Matlab. The command line will show, among other things, the suite2p version that is installed.

```matlab
py.suite2p.default_ops()
```

If your pyenv is setup correctly in Matlab, the sbx pipeline will call the suite2p pacakge automatically .

This 'InProcess' version seems to work fine on my Windows desktop, but I had trouble getting this to work on our HPC cluster. The problem appears to be that some of the libraries that Matlab uses conflict with those in the Python install. To get around this, the current code defaults to calling Python outside of Matlab (using the system() command.) This has the advantage that there are no library conflicts. For this version, you have to set the NS_CONDA environment variable. For instance, on a unix installation:

```matlab
setenv('NS_CONDA','/home/bart/miniconda3')
```

(By not setting this environment variable, or setting it to empty, the InProcess variant is used).

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
nsInitializeDataJoint(codeFolder,'ProjectName',root,'ns')
% Scan all files in the year 2023 (underneath the root folder)
nsScan('date','01-Jan-2023','schedule','y','root',root,...
                    'readFileContents',true,'addToDataJoint',true,'safeMode',false);
% The database tables have now been created and filled:
ns.Experiment % Shows the Experiment table in the ProjectName project
% The get() function of the ns.Experiment class gives access to all Plugin properties in the experiment.
% For instance in an experiment where the orientation of the grating stimulus changed in each trial, the following call retrieves the orientation for each trial. 
ori = get(ns.Experiment ,'gabor','prm','orientation','atTrialTime',0)
```


BK -  April 2023.
