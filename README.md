# DataJoint schema and tools for Neurostim data files
This folder contains scripts to start a Datajoint schema and database analysis pipeline for 
data collected with Neurostim (https://github.com/klabhub/neurostim).

To read more about DataJoint, how to setup a DataJoint SQL server, and to install the DataJoint toolbox for Matlab, see https://datajoint.com.

Once DataJoint is installed here are the functions that setup a pipeline:

1. nsInitializeDataJoint  (create a schema/database in the SQL database)
1. nsScan                 (scan a folder with Neurostim output files)  
1. nsAddToDataJoint       (Add the scanned files)

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