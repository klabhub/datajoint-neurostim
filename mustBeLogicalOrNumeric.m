function mustBeLogicalOrNumeric(x)

assert(isnumeric(x) && (x == 0 || x == 1) || islogical(x), "Must be logical or 0 or 1.");

end