function writeJson(str,file)
%Utility to write human readable json files
str = jsonencode(str,"PrettyPrint",true);
fid  = fopen(file,'wt');
fprintf(fid,'%s',str);
fclose(fid);
end
