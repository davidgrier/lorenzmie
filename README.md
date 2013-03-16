# lorenzmie

**IDL routines for analyzing in-line holographic microscopy images
using the Lorenz-Mie theory of light scattering.**

IDL is the Interactive Data Language, and is a product of
[Exelis Visual Information Solutions](http://www.exelisvis.com)

lorenzmie is licensed under the GPLv3.

## What it does

lorenzmie is intended for computing and analyzing holographic images
of micrometer-scale spheres.  It assumes that the illuminating laser beam
is collimated along the optical axis and that the light is monochromatic,
coherent, and uniformly polarized along the x (horizontal) axis.  A sphere
is specified by its three dimensional position, its radius and its 
complex refractive index.  Spheres are assumed to have isotropic 
and linear optical properties.  They may be radially stratified with layers
defined by their radii and refractive indexes.  The surrounding medium
is characterized by its complex refractive index, and is assumed to be
isotropic and homogeneous.

The code is organized from the most general components to the most
specialized:

* **sphericalfield**: Calculates the complex electric field at a displacement
(x,y,z) relative to the center of an illuminated object whose scattering
properties are defined by Lorenz-Mie scattering coefficients.

* **gpu_sphericalfield**: A hardware-accelerated implementation of
sphericalfield based on the [GPULib](http://www.txcorp.com/home/gpulib)
library of IDL bindings to
[CUDA](http://www.nvidia.com/object/cuda-parallel-computing-platform.html).

* **sphere_coefficients**: Calculates the Lorenz-Mie scattering
coefficients for a multilayered sphere of radius ap and refractive
index np immersed in a medium of refractive index nm at vacuum
wavelength lambda.  These coefficients then can be used in
sphericalfield to compute the field scattered by a stratified sphere.

* **spherefield**: Uses sphere_coefficients and either sphericalfield or
gpu_sphericalfield to compute the field scattered by a stratified sphere.

* **lmsphere/**: Computes normalized holograms of spheres based on the
field computed by spherefield.  Also fits experimentally measured
holograms of colloidal spheres to the predictions of Lorenz-Mie theory.

### References

1. S.-H. Lee, Y. Roichman, G.-R. Yi, S.-H. Kim, S.-M. Yang,
A. van Blaaderen, P. van Oostrum and D. G. Grier,
"Characterizing and tracking single colloidal particles with video
holographic microscopy," 
_Optics Express_ **15**, 18275-18282 (2007).

2. F. C. Cheong, B. Sun, R. Dreyfus, J. Amato-Grill, K. Xiao,
L. Dixon and D. G. Grier, "Flow visualization and flow cytometry with
holographic video microscopy," _Optics Express_ **17**
13071-13079 (2009).

3. F. C. Cheong, B. J. Krishnatreya and D. G. Grier,
"Strategies for three-dimensional particle tracking with
holographic video microscopy,"
_Optics Express_ **18**, 13563-13573 (2010).

4. H. Moyses, B. J. Krishnatreya and D. G. Grier,
_Optics Express_ **21** 5968-5973 (2013).

5. C. F. Bohren and D. R. Huffman, Absorption and Scattering of Light
by Small Particles (New York, Wiley 1983).

6. W. Yang, "Improved recurstive algorithm for light scattering
by a multilayered sphere," _Applied Optics_ **42**, 1710--1720 (2003).

7. O. Pena and U. Pal, "Scattering of electromagnetic radiation
by a multilayered sphere," _Computer Physics Communications_
**180**, 2348-2354 (2009).
