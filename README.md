SEP
===

C library for Source Extraction and Photometry

*"... [it's] an SEP: Somebody Else's Problem."  
"Oh, good. I can relax then."*

[Source Extractor](http://www.astromatic.net/software/sextractor) is
great, but sometimes you just want to use some of the pieces from it
without running the entire executable. This library implements a few of
the algorithms used in SExtractor as stand-alone pieces.

Install
-------

The library has no external dependencies.

API
---

### Background estimation

```c
backmap *makeback(float *im, float *var, int w, int h,
                  int bw, int bh, float varthresh, int fbx, int fby,
                  float fthresh, *status);
```
*create a representation of the image background and its variance*

Note that the returned pointer must be freed by calling `freeback()`.

* `im` : image data  
* `var` : variance/mask data (ignored if NULL)  
* `w`, `h` : width and height of image and variance arrays
  (width is the fast axis)  
* `bw`, `bh` : box size (width, height) used for estimating background
  [default in SExtractor: 64]  
* `varthresh` : pixels with `var > varthresh` are ignored.  
* `fbx`, `fby` : Filter size in x and y directions used in median-filtering
  background meshes [SE default: 3]  
* `fthresh` : Threshold used in filtering [SE default: 0.0]
* `status` : Set to nonzero if error.

```c
void backline(backmap *bkmap, int y, float *line);
```

*Evaluate the background using bicubic spline interpolation at line y*

* `bkmap` : pointer to background map
* `y` : index of the line  
* `line` : array of size `bkmap->imnx` (image width)

```c
void backim(backmap *bkmap, float *arr);
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
void backrmsim(backmap *bkmap, float *arr);
```

*same as `backim()` but for background RMS*

```c
float backpixlinear(backmap *bkmap, int x, int y);
```

*return background at position x,y using linear interpolation between
background map vertices*

```c
void freeback(backmap *bkmap);
```

*free memory associated with a background spline*

### Object detection

```c
objliststruct *extract(PIXTYPE *im, PIXTYPE *var, int w, int h,
		       PIXTYPE dthresh, PIXTYPE athresh, PIXTYPE cdwthresh,
		       int threshabsolute, int minarea,
		       float *conv, int convw, int convh,
		       int deblend_nthresh, double deblend_mincont,
		       int clean_flag, double clean_param, int *status);
```

*Detect objects in an image*

* `im`: image array
* `var`: variance array (can be NULL)
* `w`, `h`: width and height of arrays (width is fast axis)
* `dthresh`: detection threshold [SE default: 1.5]
* `athresh`: analysis threshold [SE default: 1.5]
* `cdwthresh`: ?
* `threshabsolute`: Use abolute thresholding (0=relative) [SE default: 0] 
* `minarea`: Minimum number of pixels for detection [SE default: 5] 
* `conv`: convolution array (can be NULL) [SE default: {1 2 1 2 4 2 1 2 1}] 
* `convw`, `convh`: size of convolution array [SE default: 3, 3]
* `deblend_nthresh`: Number of thresholds in deblending [SE default: 32]
* `deblend_mincont`: Deblending parameter [SE default: 0.005]
* `clean_flag`: Perform cleaning? (1 = YES) [SE default: 1]
* `clean_param`: Cleaning parameter [SE default: 1.0]

Tests
-----



License
-------

As Source Extractor's license is GPLv3 and this is a derived work,
the license for SEP is also GPLv3.

