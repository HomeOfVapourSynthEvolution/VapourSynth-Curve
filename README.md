Description
===========

Apply color adjustments using curves.

This filter is similar to Adobe Photoshop and GIMP curves tools. Each plane has its values defined by *N* key points tied from each other using a smooth curve. The x-axis represents the pixel values from the input frame, and the y-axis the new pixel values to be set for the output frame.

By default, a curve is defined by the two points *(0;0)* and *(1;1)*. This creates a straight line where each original pixel value is "adjusted" to its own value, which means no change to the image.

The filter allows you to redefine these two points and add some more. A new curve (using a natural cubic spline interpolation) will be defined to pass smoothly through all these new coordinates. The new defined points need to be strictly increasing over the x-axis, and their x and y values must be in the *[0;1]* interval. If the computed curves happened to go outside the vector spaces, the values will be clipped accordingly.


Usage
=====

    curve.Curve(clip clip[, int preset=0, float[] r=None, float[] g=None, float[] b=None, float[] master=None, string acv=None, int[] planes=[0, 1, 2]])

* clip: Clip to process. Any planar format with integer sample type of 8-16 bit depth is supported.

* preset: Selects one of the available color presets. This parameter can be used in addition to the `r`, `g`, `b` parameters; in this case, the later parameters takes priority on the preset values. Note that the values of preset 1, 2, 10 are defined in RGB and should not be applied to YUV clip. The other presets can be used on YUV clip, but only the first plane should be applied, except preset 8.
  * 0 = none
  * 1 = color negative
  * 2 = cross process
  * 3 = darker
  * 4 = increase contrast
  * 5 = lighter
  * 6 = linear contrast
  * 7 = medium contrast
  * 8 = negative
  * 9 = strong contrast
  * 10 = vintage

* r: Sets the key points for the red/y plane.

* g: Sets the key points for the green/u plane.

* b: Sets the key points for the blue/v plane.

* master: Sets the master key points. These points will define a second pass mapping. It is sometimes called a "luminance" or "value" mapping. It can be used with `r`, `g`, `b` since it acts like a post-processing LUT.

* acv: Specifies a Photoshop curves file (.acv) to import the settings from.

* planes: Sets which planes will be processed. Any unprocessed planes will be simply copied.


Examples
========

* Increase slightly the middle level of blue: ```curve.Curve(clip, b=[0,0, 0.5,0.58, 1,1])```

* Vintage effect: ```curve.Curve(clip, r=[0,0.11, 0.42,0.51, 1,0.95], g=[0,0, 0.5,0.48, 1,1], b=[0,0.22, 0.49,0.44, 1,0.8])```


Compilation
===========

```
meson build
ninja -C build
ninja -C build install
```
