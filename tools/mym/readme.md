These are the compiled binaries available in the mym project on github. 
Copying these avoids the need to install GHTB ,compareVersions etc. 
The user only needs to add the relevant subfolder to their path.

Oct 2026:
Forked [the original mym ](https://github.com/datajoint/mym) to [klabhub](https://github.com/klabhub/mym).
Added cleanup on exit (closing connections).
Recompiled on RHEL and Windows. 
This fixed matlab crashes that only occurred on exiting Matlab on RHEL.

Fixed "Commands out of sync" when a query contains multiple statements (e.g. a
multi-statement `connectionInit_function`): mym now consumes the results of all
statements (`drainResults` in mym.cpp). Rebuilt on RHEL (mexa64) only; mexw64 and mexmaci64 still need to be recompiled.
