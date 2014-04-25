SEP
===

C library for Source Extraction and Photometry

*"... [it's] an SEP: Somebody Else's Problem."  
"Oh, good. I can relax then."*

[Source Extractor](http://www.astromatic.net/software/sextractor) is
great, but sometimes you just want to use some of the pieces from it
without running the entire executable. This library implements a few of
the algorithms used in SExtractor as stand-alone pieces. So far this
includes:

* background estimation
* source detection

Install
-------

The library has no dependencies outside the C standard library.

API
---

### Background estimation

```c
int makeback(float *im, float *mask, int w, int h,
             int bw, int bh, float maskthresh, int fbx, int fby,
             float fthresh, backmap **bkmap);
```
*create a representation of the image background and its variance*

Note that the returned pointer must be freed by calling `freeback()`.

* `im` : image data  
* `mask` : mask data (ignored if NULL)  
* `w`, `h` : width and height of image and variance arrays
  (width is the fast axis)  
* `bw`, `bh` : box size (width, height) used for estimating background
  [default in SExtractor: 64]  
* `maskthresh` : pixels with `mask > maskthresh` are ignored.  
* `fbx`, `fby` : Filter size in x and y directions used in median-filtering
  background meshes [SE default: 3]  
* `fthresh` : Threshold used in filtering [SE default: 0.0]
* `bkmap` : Returned background map.

Returns: `status` error code (0 if OK).

```c
void backline(backmap *bkmap, int y, float *line);
```

*Evaluate the background using bicubic spline interpolation at line y*

* `bkmap` : pointer to background map
* `y` : index of the line  
* `line` : array of size `bkmap->imnx` (image width)

```c
void backarray(backmap *bkmap, float *arr);
```

*Evaluate the background using bicubic spline interpolation for the entire
image*

* `bkmap` : pointer to background map   
* `arr` : array of size `bkmap->imnx * bkmap->imny` (original image pixels)

```c
void backrmsline(backmap *bkmap, int y, float *line);
```

*same as `backline()` but for background RMS*

```c
void backrmsarray(backmap *bkmap, float *arr);
```

*same as `backarray()` but for background RMS*

```c
float backpixlinear(backmap *bkmap, int x, int y);
```

*return background at position x,y using linear interpolation between
background map vertices*

```c
void freeback(backmap *bkmap);
```

*free memory associated with a background map*

### Source detection

```c
int extractobj(float *im, float *var, int w, int h,
               float thresh, int minarea,
               float *conv, int convw, int convh,
               int deblend_nthresh, double deblend_mincont,
               int clean_flag, double clean_param,
               int *nobj, sepobj **objects);
```

*Detect objects in an image*

* `im`: image array
* `var`: variance array (can be NULL), used for thresholding
* `w`, `h`: width and height of arrays (width is fast axis)
* `thresh`: detection threshold [SE default: 1.5*(background sigma)]
* `minarea`: Minimum number of pixels for detection [SE default: 5] 
* `conv`: convolution array (can be NULL) [SE default: {1 2 1 2 4 2 1 2 1}] 
* `convw`, `convh`: size of convolution array [SE default: 3, 3]
* `deblend_nthresh`: Number of thresholds in deblending [SE default: 32]
* `deblend_mincont`: Deblending parameter [SE default: 0.005]
* `clean_flag`: Perform cleaning? (1 = YES) [SE default: 1]
* `clean_param`: Cleaning parameter [SE default: 1.0]

If `var` is NULL, `thresh` is taken to be be an absolute threshold in ADU.
Otherwise, `thresh` is taken to be sigma and a variable threshold of
`thresh * sqrt(var)` is calculated for each pixel.

Return values:

* `nobj` : number of objects detected
* `objects` : array of `sepobj` structs of length `nobj`

An `sepobj` struct holds a collection of parameters characterizing the
size, shape, and brightness of the object. It also contains an array
of ints giving the pixels belonging to the object. Note that these arrays
must eventually be `free`d. 

Running Tests
-------------

In `test` directory, type `make` then run the executable `runtests`. 


Speed
-----

For a 2k x 4k image with ~2000 sources, `makeback` takes ~250ms and
`extractobjs` takes ~450ms, with SE default settings.
Tested on a 1.7 GHz Core i5 Ivybridge laptop.

License
-------

As Source Extractor's license is GPLv3 and this is a derived work,
the license for SEP is also GPLv3.

