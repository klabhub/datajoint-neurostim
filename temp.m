function v =temp(td)
% A function that behaves just as the built-in tempname, except that it
% returns names in the folder specified in the NS_TEMPDIR environment
% variable (if it is defined).
%
arguments
    td  (1,1) string =""
end
if td==""
    td = string(getenv('NS_TEMPDIR'));
end

if td==""
    v = tempname;
else
    v = tempname(td);
end
