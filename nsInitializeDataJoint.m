function nsInitializeDataJoint(code,databaseName,packageName,pv)
% Setup a new datajoint pipeline for a Neurostim project.
% For instance, if Alice's data are all under the root folder x:\data\
% and she wants to start a new project called 'alice_memory' with Matlab code
% for the pipelin stored in u:\projects\memory\code, she uses:
%
%  nsInitializeDataJoint('u:\projects\memory\code','alice_memory','x:\data\')
%  And then she runs
%  nsScan('date','01-May-2022','schedule','m','readFileContents',true)
% to scan files collected in May '22 and add them to the pipeline.
% The code in the datajoint-neurostim repository handles file scanning and
% stores all Neurostim data (i.e. all neurostim parameters) in the
% database.
% To retrieve the data for a single Neurostim experiment (i.e. file), you
% use a query to define the experiment :
% For instance, this will return a struct array with the data for all
% experiments for subject #12:
% data = get(ns.Experiment & 'subject=12')
%
% Alice can now add files to the u:\projects\memory\code\+ns that define
% the analysis pipeline.
%
% If Alice also wants to analyze Calcium imaging data, she can add the Ca
% package by calling
% nsInitializeDataJoint('u:\projects\memory\code','alice_memory','','ca')
% This will create the lookup tables from CA Element in the alice_memory
% database/schema, and create a +ca folder in the /code folder where Alice
% can add her own additions to the analysis pipeline.
%
% INPUT
% code -  The root level folder where your project's code lives. The script
%           will create a +ns package subfolder with a getSchema.m file.
% databasename - The name of this project in your SQL database (e.g.
%                   'alice_memory')
% packageName- The package(s) to install in the database. Default is ns
% (package to handle Neurostim files). Other options are:
%               'sbx' - Calcium imaging with ScanBox
%               'ephys' -  EEG , LFP electrophysiology 
%               To use multiple packages, pass a cell array of strings with
%               package names ({'ns','sbx'}). Each package will create a
%               separate schema within the same database
% Parm/Value pairs:
% 'dataRoot' - The root folder that contains all data files. The
%               folders below this are the years. Defaults to
%               getenv('NS_ROOT')

%  BK - April 2022

arguments
    code {mustBeText} 
    databaseName {mustBeText}    
    packageName {mustBeText} = {'ns'}
    pv.dataRoot {mustBeText} = getenv('NS_ROOT');
end
% Make sure the database name begins with a lower case letter
if isempty(regexp(databaseName,'^[a-z]+', 'once'))
    error('The schemaName (%s) must start with a lower case character\n', databaseName);
end
if ~exist(code,'dir')
    mkdir(code);
end
if ~iscell(packageName)
    packageName  = {packageName};
end

%% Add the schema and utilities (dj*)
% Move to the code folder and add it to the path
cd(code)
%% Define the ROOT as an environment variable
if ~isempty(pv.dataRoot)
    setenv('NS_ROOT',pv.dataRoot);
end
%% Create the schema on the SQL server

query(dj.conn, sprintf('CREATE DATABASE IF NOT EXISTS `%s`',databaseName))

for i=1:numel(packageName)
    %% Create a package folder in the project to extend the schema
    packageDir =fullfile(code,['+' packageName{i}]);
    if ~exist(packageDir,"dir")
        mkdir(packageDir);
    end

    %% Create the getSchema function - use databaseName/packageName to setup different schemas for each package.
    gs= 'function obj = getSchema\n persistent OBJ \n if isempty(OBJ) \n     OBJ = dj.Schema(dj.conn,''%s'', ''%s/%s'');\n end\n obj = OBJ;\n end \n';
    fid = fopen(fullfile(code,['+' packageName{i}],'getSchema.m'),"w");
    fprintf(fid,gs,packageName{i},databaseName,packageName{i});
    fclose(fid);
    switch (packageName{i})
        case 'ns'          
            fprintf('The datajoint pipeline for %s has been setup. Run nsScan to add files.\n',databaseName);     
        case 'sbx'
            % Nothing to do
        otherwise
    end
end
%% Go to the project folder
cd(code)




end