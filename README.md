# lorenzmie

**IDL routines for analyzing in-line holographic microscopy images
using the Lorenz-Mie theory of light scattering.**

IDL is the Interactive Data Language, and is a product of
[Exelis Visual Information Solutions](www.exelisvis.com)

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
sphericalfield based on the [GPULib](www.txcorp.com/home/gpulib)
library of IDL bindings to
[CUDA](www.nvidia.com/object/cuda-parallel-computing-platform.html).

* **sphere_coefficients**: Calculates the Lorenz-Mie scattering
coefficients for a multilayered sphere of radius ap and refractive
index np immersed in a medium of refractive index nm at vacuum
wavelength lambda.  These coefficients then can be used in
sphericalfield to compute the field scattered by a stratified sphere.

* **spherefield**: Uses sphere_coefficients and either sphericalfield or
gpu_sphericalfield to compute the field scattered by a stratified sphere.

