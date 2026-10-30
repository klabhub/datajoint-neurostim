function v =temp(td)
% A function that behaves just as the built-in tempname, except that it
% returns names in the folder specified in the NS_TEMPDIR environment
% variable (if it is defined).
%
arguments
    td  (1,1) string =string(getenv('NS_TEMPDIR'))
end

if td==""
    v = tempname;
elseif exist(td,"dir")
    v = tempname(td);
else
    fprintf('This NS_TEMPDIR does not exist: %s. Using %s instead.\n',td,tempdir)
    v = tempname;
end
