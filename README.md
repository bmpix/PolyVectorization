This is an implementation of "Vectorization of Line Drawings via PolyVector Fields" by Mikhail Bessmeltsev and Justin Solomon, Massachusetts Institute of Technology, 2018. 
Actual coding done by Mikhail Bessmeltsev, http://www-labs.iro.umontreal.ca/~bmpix/.

USAGE:
./polyvector(.exe) filename

OUTPUT:
Creates a filename.svg in the same folder. 

## Requirements

(Other versions might also work, but these are the ones I used)
- QT 5.7 (for the version with GUI, see line 2 of main.cpp)
- OpenCV 2.4
- Boost 1.65.1
- Eigen 3.3.1

## Citation

When using the code in your research work, please cite the following paper:

    @article{Bessm2019,
    author = {Bessmeltsev, Mikhail and Solomon, Justin},
    title = {Vectorization of Line Drawings via Polyvector Fields},
    journal = {ACM Trans. Graph.},
    issue_date = {February 2019},
    volume = {38},
    number = {1},
    month = jan,
    year = {2019},
    pages = {9:1--9:12},
    articleno = {9},
    numpages = {12},
    publisher = {ACM},
    address = {New York, NY, USA},
    keywords = {Vectorization, frame field, line drawing, polyvector field},
    } 


## Building

1. Create folder 'build' in the same folder as this README.  
2. In cmd (windows) or sh/bash (linux/mac):  
   > cd build  
   > cmake ..  
   > make  

(If you are on Windows, the last instruction will depend on the compiler you're using, e.g. 'nmake' for Visual Studio)

## Parameters

1. All the tunable parameters are in Params.h
2. In a couple of places in the code there are seemingly random constants like 10 - those are implementation-dependent constants and have nothing to do with algorithm parameters. Proper engineering would fix those, but before that happens - please do not change.

## Notes

1. This is not the optimized version we tested for performance. The optimized version is available upon request.
2. There are a few minor changes from the paper description:
	- Section 5.2 has been reimplemented since the submission for robustness. This might have introduced minor (~2 pixel)-differences to some of the results, mostly making them better.
	- Instead of outputting a raw non-smooth vectorization with many segments, which takes a while to output, we use Douglas-Peucker algorithm and then Laplacian smoothing on the result. Douglas-Peucker does not change anything significant in the result, except for the density of the control points. Laplacian smoothing was not used for the paper results, instead we used, as we noted, Adobe Illustrator's 'Simplify' feature. Those two steps were added to immediately output a decent .svg if Illustrator is not available.
3. If you want to build a command-line tool without GUI, comment out this line in the beginning of main.cpp:
     >  #define WITH_GUI 1
	
## Known Issues

1. Output might contain 'nan's. Your importer/viewer should ignore those points.