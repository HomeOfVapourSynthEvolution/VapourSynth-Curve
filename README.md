Description
===========

Apply color adjustments using curves.

This filter is similar to Adobe Photoshop and GIMP curves tools. Each plane has its values defined by *N* key points tied from each other using a smooth curve. The x-axis represents the pixel values from the input frame, and the y-axis the new pixel values to be set for the output frame.

By default, a curve is defined by the two points *(0;0)* and *(1;1)*. This creates a straight line where each original pixel value is "adjusted" to its own value, which means no change to the image.

The filter allows you to redefine these two points and add some more. A new curve (using a natural cubic spline interpolation) will be defined to pass smoothly through all these new coordinates. The new defined points need to be strictly increasing over the x-axis, and their x and y values must be in the *[0;1]* interval. If the computed curves happened to go outside the vector spaces, the values will be clipped accordingly.


Usage
=====

    curve.Curve(clip clip[, int preset=0, string[] curve=['', '', ''], string master='', string acv='', int[] planes=[0, 1, 2]])

* clip: Clip to process. Any planar format with integer sample type of 8-16 bit depth is supported.

* preset: Selects one of the available color presets. Can be used in addition to the `curve` parameter. In this case, the later parameter takes priority on the preset values. Note that the values of preset 1, 2, 10 are defined in RGB and should not be applied to YUV clip. The other presets can be used on YUV clip, but only the first plane should be applied except preset 8.
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

* curve: Sets the key points for each plane. The syntax of key points list is `x0/y0 x1/y1 x2/y2 ...`. However, the delimiter character is arbitrary. `x0;y0 x1:y1 x2-y2 ...` works just fine.

* master: Sets the master key points. These points will define a second pass mapping. It is sometimes called a "luminance" or "value" mapping. It can be used with `curve` since it acts like a post-processing LUT.

* acv: Specifies a Photoshop curves file (.acv) to import the settings from.

* planes: Sets which planes will be processed. Any unprocessed planes will be simply copied.


Compilation
===========

```
meson build
ninja -C build
ninja -C build install
```
