function assertRequiredMembers(parms,required,fun)
% Convenience function to alert the user to required parameters for a
% ns.CParm. Call this from the function that is defined as the 'fun' in
% ns.CParm. 
%
% See Also : sbx.spikeML
arguments
    parms (1,1) struct
    required (1,:) string
    fun (1,1) string
end
specified=fieldnames(parms);
missing = setdiff(required,specified);
if ~isempty(missing)
    fprintf(2,"ns.C processing function %s has missing members in the ns.CParm parms struct:\n",fun);
    fprintf(2,"* %s\n",missing);
    error("ns.C - %s failed",fun);
end
