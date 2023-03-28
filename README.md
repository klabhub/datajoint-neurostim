# DataJoint schema and tools for Neurostim data files

This folder contains scripts to start a Datajoint schema and database analysis pipeline for
data collected with Neurostim (<https://github.com/klabhub/neurostim>).

To read more about DataJoint, how to setup a DataJoint SQL server, see <https://datajoint.com>.

## Installation

1. Clone this reposiory to a folder on your machine (`'c:\github\datajoint-neurostim'`)
1. Add it to your matlab path (`addpath('c:\github\datajoint-neurostim'`)
1. This should work with the regular datajoint toolbox from datajoint.com, but the repo also contains a fork of this toolbox that may have some additional features. You can use this fork instead by adding the datajoint-matlab submodule to your path. ((`addpath('c:\github\datajoint-neurostim\datajoint-matlab'`)
1. To try if it worked, go to a temp folder and initialize a database:

```matlab
cd c:\temp;
root = 'the root folder containing your neurostim data files. This folder has years as subfolders.';
 nsInitializeDataJoint(pwd,'mytest',root,'ns')
```

### Scanbox
The +sbx package sets up a pipeline for two-photon imaging data collected with the ScanBox software (Neurolabware.com). The preprocessing uses the [suite2p] (https://github.com/MouseLand/suite2p#readme), which requires Python 3.8 to run.

1. Install miniconda
1. Follow the instructions from the suite2p readme to install suite2p in its own conda environment named 'suite2p'
1. Tell Matlab where to find Python. For instanc, I did this:
```matlab
pyenv('ExecutionMode','InProcess','Version','c:\Users\bartk\miniconda3\envs\suite2p\python.exe');
```
To check that this works. Load the default suite2p options. This should return a Python dict in Matlab. The command line will show, among other things, the suite2p version that is installed.
```matlab
py.suite2p.default_ops()
```

The Matlab pipeline will call the pacakge automatically if your pyenv is setup correctly in Matlab.
## Usage

Once DataJoint is installed here are the functions that setup a pipeline:

1. nsInitializeDataJoint  (create a schema/database in the SQL database)
1. nsScan                 (scan a folder with Neurostim output files and add them to the database)  

This can be scripted, or run interactively from the nsMeta Matlab app, which also allows you to add additional meta information (which will be written to json files in the data folder structure).

The elements of these additional meta data can be pre-specified for subjects, sessions, and experiments in a _definition.json file, or they can be defined ad-hoc in the nsMeta app. For instance, consider this element in the session_definition.json file (placed at the root of the data folder):

```
{
 "Fields":{
"comments":{"Description":"Free form comments"}},
 "Locked":1
} 
```

This states that for each Session meta data representing "comments" can be added but no other meta data are allowed (Locked:1) for the experiments below the root folder where the _definition file is stored.  See subject_definition.json and experiment_definition.json for additional examples.

BK -  January 2023.
