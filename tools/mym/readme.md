These are the compiled binaries available in the mym project on github. 
Copying these avoids the need to install GHTB ,compareVersions etc. 
The user only needs to add the relevant subfolder to their path.

Oct 2026:
Forked [the original mym ](https://github.com/datajoint/mym) to [klabhub](https://github.com/klabhub/mym).
Added cleanup on exit (closing connections).
Recompiled on RHEL and Windows. 
This fixed matlab crashes that only occurred on exiting Matlab on RHEL.
