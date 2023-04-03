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
1. Tell Matlab where to find Python. For instanc, I did this:

```matlab
pyenv('ExecutionMode','InProcess','Version','c:\Users\bartk\miniconda3\envs\suite2p\python.exe');
```

To check that this works. Load the default suite2p options. This should return a Python dict in Matlab. The command line will show, among other things, the suite2p version that is installed.

```matlab
py.suite2p.default_ops()
```

If your pyenv is setup correctly in Matlab, the sbx pipeline will call the suite2p pacakge automatically .

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

## nsMeta
The `nsMeta` Matlab app provides a graphical interface to the metadata of a Neurostim data tree. This information is partially represented in the file and folder names (e.g., subject identifiers and paradigm names are in the Neurostim output file), but can be complemented with additional meta information stored in JSON files. 

### Meta data definitions
Meta data for an experiment are stored in a JSON file with the same name as the Neurostim output file, but with the .json extension.  
Meta data for a session (i.e., the unique combination of a date and a subject) is stored in a filed called `subject.json`, stored in the day of the session.
Meta data for all subjects are stored at the root of the data tree in a file called `subject.json`.

The elements of these additional meta data can be pre-specified for subjects, sessions, and experiments in a _definition.json file placed at the root of the data tree. This implies that the same default set of metad data fields is used for all Neurostim data in the tree. For instance, consider this element in the session_definition.json file:

```json
{
 "Fields":{
"comments":{"Description":"Free form comments"}},
 "Locked":1
} 
```

This states that for each Session meta data representing "comments" can be added but no other meta data are allowed (`Locked:1`).  See subject_definition.json and experiment_definition.json for additional examples.
You''l find example _definition.json files in the neurostim_datajoint folder. To use them for a specific data tree, copy them to the root folder of the tree, edit to your liking, then run `nsMeta`. 

Unless the _definition files state that the meta data are `Locked`, the nsMeta app also allows users to add new meta data fields.  

### Usage

Type `nsMeta` to start the app.
Go to File|Change Root Folder (Ctrl-R) and select the root of your Neurostim data tree (the folder with subfolders named after the years).
nsMeta will scan these folders for Neurostim files, and fill the Year/Month/Day drop downs.
Select a date using the drop downs (or use the previous (<) or next (>) buttons to move to the previous/next day; days for which there is no folder in the tree will be skipped). 
If the Refresh button is green, then you should press it to scan the folder. This will populate the Experiment, Session, and Subject tables with all the meta data based on file names and json files (if any).

Clicking on an (editable) meta data field will allow you to edit it.
A double-click on the ID field will open the json file associated with that table/row (if it exists).

Once you make changes to any of the meta data, the Update button will turn green. Press Update to save changes to the meta data (they are not saved automatically, but nsMeta will warn you if you try to move to a different day without updating).

To test what would be done without making any changes, you select dryrun and read was is shown on the command line. If that looks good, deselect dryrun and press Update again.

### Interaction with DataJoint

By default, nsMeta only writes the changes to the json files in the data tree. In principle you could then run nsScan to update the DataJoint database, but it may be easier to select the DataJoint check box to do this at the same time.

For this to work, the connection parameters under File | DataJoint should be set correctly. The Safemode parameter in that menu determines whether deletes from the DataJoint database require manual confirmation. Set it to 0 to sidestep this.

In addition, the target database on the server is determined by your working folder in Matlab. After you select the DataJoint checkbox, the currently selected database will be shown to the left of it. Before Updating with DataJoint selected, make sure you are writing to the correct database!

Also, you may want to use the Paradigm or Subject selectors at the top of the interface to only add specific experiments to your database (and ignore the other experiments that may have happened on the same day).



BK -  April 2023.
