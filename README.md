Snaptool can read GADGET-2 snapshots (in formats 1 and 2, but not HDF-5), as well as ASCII files with particle X Y Z positions, and plot density maps.

Supported features:
* Read every N-th particle to reduce memory usage for huge snapshots
* Shift box center, using periodical boundary conditions
* Cut spherical region, resize and rotate the region to plot
* Change palette and image resolution
* Add a grid of lines on a plot
* Add halo positions from AHF halo catalog
* And many others

For example, Snaptool can be used to create an animation like this: [https://youtube.com/watch?v=NLlAhW7XYQQ](https://youtube.com/shorts/NLlAhW7XYQQ).

The documentation is still incomplete...

# Installation

You need gfortran. Type 'make' to compile the code.

# Usage

Prepare a file with commands, and run with:
'''
./snaptool file.plt
'''
The output will be a set of .gif files.


'''
      ---== SNAPTOOL HELP ==---

 In-file contains a sequence of commands.
 The line with the command name should not contain anything else:
 the arguments, if required, are specified on subsequent line(s).
 The in-file should finish with the end command.

 List of commands:

 # the rest of line is a comment

 set map
 isize jsize plane
   sets the map size in pixels and the projection plane,
   which is one of: xy, xz, yz.
   Default: 500 500 xy

 read pos
 path/to/snapshot_xxx.0
 ...
 path/to/shapshot_xxx.last
   reads positions of particles from either a single-file or
   multiple-file GADGET snapshot. The number of files to read
   is determined automatically from the header of the first
   file in the list. Clears the map.

 map
   fills the density map using particles currently read into memory.
   The previous contents are not cleared.

 plot
   saves the map as an image, clears the map.

 set read every
 n
   read only every n-th particle from snapshots. Default: 1.

 set center
 X Y Z
   sets center of the region to map and plot. Default: Box/2 Box/2 Box/2

 move center
 X Y Z
   the same as "set center", but affects only the positions of particles
   that are already read, and does not affects next "read pos".

 evolve center
 X1 Y1 Z1 zstart
 X2 Y2 Z2 zend
   defines the movement of the map center between two positions for each
   snapshot.

 set size
 S
   sets size of the working region in all three dimensions. Default: Box.

 set cut
   or
 unset cut
   chooses whether to cut the region of the size S and center X, Y, Z
   during reading the snapshot. This is useful to decrease memory usage.

 set N pos max
 N
   sets the maximum number of particles which can be allocated.
   Default: 2 000 000 000.

 set outfile
 filename
   sets the output filename base. Default: snapplot

 add suffix z
   adds _zX.XXX to the output filenames.
 add suffix axes
   adds _xy, _xz or _yz to the output filenames.
 add suffix num
   adds _NNNNN to the output filenames. The number is increased
   each time writing occures.

 include
 filename
   run the sequence of commands from the given file,
   than return to the current file.

 select sphere
 R
   marks all particles with distance >R from center
   as non-plottable.

 select type
 i
   selects only particles of type i=0..5

 print dist
   prints distance from previously set center to the
   closest particle from current selection.

 rotate
 angle x y z
   rotate everything on the given angle around the axis given
   by the three coordinates.

 AHF
 file
   display halos from AHF _halos catalog

 set palette
 N
   N=1 (default) - black-red-yellow. N>1 - white-black.

 set grid
 step
   plots a grid with given step size and origin at box origin.

 load ramses map
 filename
   load ascii file created by RAMSES part2map utility.

'''
