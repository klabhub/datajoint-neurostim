

The 64 bit version has a slightly different header than the 32bit version, that's why the sources need to be kept in 
separate directories.

To generate the mex files:

switch (mexext)
case 'mexw32'
    cd win32
	mex -largeArrayDims win32/edfmex.cpp win32/edfapi.lib -Iwin32
case 'mexw64'
    cd win64
	mex -largeArrayDims edfmex.cpp libedfapi64.lib -Iwin64
case 'mexa64'

   mex -largeArrayDims edfmex.cpp libedfapi.so -Ia64
   
end 


Use the makeheader.m file in the source directories to generate the appropriate header for a machine.

Tested 32 bit and 64 bit windows; worked ok.

BK - Oct 2008.


To compile on linux, I used matlab r20178b and gcc 4.9.4. The gcc module must remain loaded
Also, to find the .so , the LD_LIRBRARY_PATH needs to include the +eyelink folder 
This can be set in .bashrc  with something like
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/datajoint-neurostim/+ns/+eyelink








