SEP
===

C library for Source Extraction and Photometry

*"... [it's] an SEP: Somebody Else's Problem."  
"Oh, good. I can relax then."*

[Source Extractor](http://www.astromatic.net/software/sextractor) is
great, but sometimes you just want to use some of the pieces from it
without running the entire executable. This library implements a few of
the algorithms used in SExtractor as stand-alone pieces.

API
---

### Background estimation

```c
backspline *makeback(float *im, float *weight, int w, int h,
                     int bw, int bh, float wthresh, int fbx, int fby,
                     float fthresh);
```
*create a spline representation of the image background and its variance*

`im` : image data  
`weight` : weight map (ignored if NULL)  
`w`, `h` : width and height of image and weight arrays
(width is the fast axis)  
`bw`, `bh` : box size (width, height) used for estimating background
(default in SExtractor is 64)  
`wthresh` : unintuitively, weights *above* wthresh are ignored **change this?**  
`fbx`, `fby` : Filter size in x and y directions used in median-filtering
background meshes (default in SExtractor is 3)  
`fthresh` : Threshold used in filtering (default in SExtractor is 0.0)

Note that the returned pointer must be freed using `freeback`.

```c
void backline(backspline *bkspl, int y, float *line);
```

*Evaluate the background using bicubic spline interpolation at line y*

`bkspl` : background spline  
`y` : index of the line  
`line` : array of size `bkspl->imnx` (image width)

```c
void backim(backspline *bkspl, float *arr);
```

*Evaluate the background using bicubic spline interpolation for the entire
image*

```c
void backrmsline(backspline *bkspl, int y, float *line);
```

*same as `backline()` but for background RMS*

```c
void backrmsim(backspline *bkspl, float *arr);
```

*same as `backim()` but for background RMS*

```c
float backpixlinear(backspline *bkspl, int x, int y);
```

*return background at position x,y using linear interpolation between
background map vertices*

```c
void freeback(backspline *bkspl);
```

*free memory associated with a background spline*

### Object detection

Nothing here yet.

### Photometry

Nothing here yet.

License
-------

As Source Extractor's license is GPLv3 and this is a derived work (though
heavily modified) the license for SEP is also GPLv3.

