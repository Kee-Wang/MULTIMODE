# MULTIMODE

## Indroduction

MULTIMODE is a general code that obtains the ro-vibrational eigenvalues and eigenfunctions of the “Watson Hamilonian” using the n-mode representation of the potential.  It has been applied to many polyatomic molecules from triatomics, e.g.,  H2O, tetratomics, e.g., H2CO, HOCO, H2O2, pentatomics, e.g., CH4, and larger molecules, e.g., C2H4, CH3OH, CH3NO2, CH3CHOO.  The approach is based on the VSCF/VCI method, which is in principle “exact”.

## Tested Compilation environment
* Our MM was tested with `gfotran` or `ifort` compilor in Linux system. 
* There is known compatibility problem compling using MacOS with `GNU Fortran (GCC) 6.1.0`, and we don't have a soluiton for that yet.

## A Quick Usage Guide -- Command Line version
Details of usage is elaborated in `/mm/src`, and here we only provide a glance of general procedure.

* Go to `/mm/src` to link user's potential in `user.XXX.f`
* Compile with provided `Makefile` to get executable `mm.x`
* User has to modify input file `fort.1` according to user's need
* Execute `./mm.x fort.1 fort.2` with `fort.2` the output file name
* Wait for program to end and get the result

**Common problems are discusseed in [`QnA.md`](https://github.com/Kee-Wang/MULTIMODE/blob/master/QnA.md) file.**

## A Quick Usage Guide -- try MM's Little Helper GUI!
Not all parameters in `fort.1` are essential. Some of the parameters can be hidden and some parameters has annoying correlations relationship with other parameters, which is hard to keep track for beginners. The difficulty has been noticed for a long time, which finally urges us to develp a user-frendly GUI version to help generate `fort.1`.

Meet `MM's little helper!` This is a GUI version of MM input file that hide non-essential parameters and have nicer user-interface to generate a comelete version of `fort.1` can be readily executed by `mm.x`.

Current'y only `generate fort.1` button is working. Other new features will be coming! Here is the screen shot:

![Image of MM's little helper](https://github.com/Kee-Wang/MULTIMODE/blob/master/mmhelper/mm_helper_screen_shot.png)





## TO-DO-LIST
### Essential

- [ ] A more detailed TO-DO-LIST
- [ ] A more detailed README.md

### Open Issue

- [ ] Should we only present 5.1.4 version


### Optional

- [ ] [Choose a liscence for MM](https://choosealicense.com)


# Version Updates

## Introduction

MM has been continuously developed since 1996. Initially MM was code for group research use only. Now the code has been made publicly available since version 5.1.4, thanks to the funding from NSF.

## Public version updates [Current project]
These are user-friendly versions, with current version first.

### ver. 5.1.4 Oct. 8, 2017
The first publicly-avaible version. The package was uploaded to Github with a myrid of supplementary materials to make the MM as user-friendly as possible.


## Important pre-public versions [Histroy]

### ver. 5.1.4
The package is updated with new features.

### ver. 4.9.0
One of most extensively used stable versions for years in the group. Lots of new features added.

### ver. 3.4
Lots of notes added.




# Ackowledgement

## Funding Agency
The original version was created by Stuart Carter, but the public user-friendly version was funded by NSF.

## Credit

### Theory
* Joel M. Bowman
* Stuart Carter
* TBD

### MM orginial package and documentation
* Stuart Carter

### User-friendly Realization
* TBD
