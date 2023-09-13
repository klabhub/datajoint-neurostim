%{
#  Table with analysis parameter sets.  
    tag:  varchar(128)    # Short name to idenfity this set.   
    ---
    target: varchar(128) # Name of the table to which this parameter set applies   
    description  = NULL : varchar(1024)  # Description of the parameter (optional) 
    parms : longblob  # struct of all parameters
%}
% 
% Tables that depend on parameters to be filled use this table as a foreign
% key, and in the keySource specifiy that only a subset of rows (defined by
% the target key) applies. 
% 
% The advantage of using such a generic table with all kinds of parameters
% is that not every table that depends on parameters has to define its own
% lookup table. 
% 
% The disadvantage is that, because the actual parms are stored as blob,
% one cannot query on their contents. For this reason, using a good 'tag'
% and 'Description' is helpful. Also, a computed table that has a foreign key in
% this table must define a getKeySource function to restrict the rows from
% this table. 
% 
% For example, see sbx.Tuning 
% 
% BK - Sept 2023
classdef Parms < dj.Lookup
    
end