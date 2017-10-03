This is a file that contains commonly encountered problems when using MM.


# About user input files

## Can I change file names?
There are only serveral files we need to deal with:  
* `user.XXX.f`: user-specified file for providing PES
* `mm.x`: executable after compilation, takes two arguments, with usage `./mm.x fort.1 fort.2`
* `fort.1`: first argument of `mm.x`, input file
* `fort.2`: second argument of `mm.x`, output file.

Their names can all be user-specified. However, by convention we usually only change `user.XXX.f` to a system of intesest, such as `user.H2O.f`, and leave other filenames along. From now on, we use these defualt file names.


## What are the files that I need to provide?

There are 3 steps you have to do to run MM:
1. Provide PES in `user.XXX.f` to get executable `mm.x`
2. Provide parameters in input file `fort.1`
3. Execute `./mm.x fort.1 fort.2`





