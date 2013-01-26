;+
; NAME:
;    lmsphere
;
; PURPOSE:
;    Computes Lorenz-Mie microscopy image of a sphere
;    immersed in a transparent medium.
;
; CATEGORY:
;    Holographic microscopy
;
; CALLING SEQUENCE:
;    holo = lmsphere(rp, ap, np, nm, lambda, mpp, dim)
;
; INPUTS:
;    rp  : [x,y,z] 3 dimensional position of sphere relative to center
;          of image.
;    ap  : radius of sphere [micrometers]
;    np  : (complex) refractive index of sphere
;    nm  : (complex) refractive index of medium
;    lambda: vacuum wavelength of light [micrometers]
;    mpp: micrometers per pixel
;    dim : [nx,ny] dimensions of image [pixels]
;
;    NOTE: The imaginary parts of the complex refractive indexes
;    should be positive for absorbing materials.  This follows the
;    convention used in SPHERE_COEFFICIENTS.
;
; OUTPUTS:
;    holo: [nx,ny] real-valued digital holographic image of sphere.
;
; KEYWORDS:
;    alpha: fraction of incident light scattered by particle.
;           Default: 1.
;
;    delta: wavefront distortion [wavelengths]
;           Default: 0.
;
;    precision: relative precision with which fields are calculated.
;
; KEYWORD FLAGS:
;    gpu: Use GPU accelerated calculation.  This only works on
;         systems with GPUlib installed.  Requires NVIDIA graphics
;         card with CUDA installed.
;
; PROCEDURE:
;    Calls SPHEREFIELD to compute the field.
;
; REFERENCE:
;    S. Lee, Y. Roichman, G. Yi, S. Kim, S. Yang, A. van Blaaderen,
;    P. van Oostrum and D. G. Grier,
;    Chararacterizing and tracking single colloidal particles with 
;    video holographic microscopy,
;    Optics Express 15, 18275-18282 (2007)
;
; EXAMPLE:
;    Display a DHM image of a 1.5 micrometer diameter polystyrene
;    sphere (np = 1.5) in water (nm = 1.33).
;
;    IDL> tvscl, lmsphere([0,0,200], 0.75, 1.5, 1.33, 0.532, 0.135, [201,201])
;
; NOTES:
;    Allow for multiple spheres
;
; MODIFICATION HISTORY:
; Written by David G. Grier, New York University, 3/2007
; 05/25/2007 DGG: Added ALPHA keyword.
; 02/08/2008 DGG: Adopted updated SPHEREFIELD syntax:
;    separate x, y, and z coordinates.
; 10/09/2008 DGG: Added LAMBDA and PRECISION keywords
; 10/14/2008 DGG: Added GPU keyword.
; 10/16/2008 DGG: Added LUT keyword.
; 03/15/2010 DGG: Documentation cleanups.
; 11/03/2010 DGG: ALPHA defaults to 1.  Clean up image calculation.
;    Documentation fixes and formatting.
; 06/23/2012 DGG flipped sign of xp and yp for consistency with
;    fitspheredhm and dggdhmspheredhm.
; 10/15/2012 DGG renamed to lmsphere.  Removed LUT option.
;    Reorganized parameters.  Added DELTA.
; 01/26/2013 DGG Documentation updates.
;
; Copyright (c) 2007-2013 David G. Grier
;-

function lmsphere, rp, ap, np, nm, lambda, mpp, dim, $
                   alpha = alpha, $
                   delta = delta, $
                   precision = precision, $
                   gpu = gpu

COMPILE_OPT IDL2

umsg = 'USAGE: lmsphere, rp, ap, np, nm, lambda, mpp, dim'
; rp: [xp,yp,zp] position of sphere relative to center of image [pixels]
;          NOTE: We take z to denote the sphere's
;                height above the imaging plane.  In general,
;                therefore, z should be input as a positive number.
; ap: radius of sphere [micrometers]
; np: (complex) index of particle
; nm: (complex) index of medium
; lambda: vacuum wavelength of light [micrometers]
; dim: [nx,ny] dimensions of image to create [pixels]

if n_params() ne 7 then begin
   message, umsg, /inf
   return, -1
endif

if ~isa(rp, /number, /array) then begin
   message, umsg, /inf
   message, 'RP should be a 3-element array of (x,y,z) coordinates', /inf
   return, -1
endif

if ~isa(ap, /number, /scalar) then begin
   message, umsg, /inf
   message, 'AP should be the radius of the sphere in micrometers', /inf
   return, -1
endif

if ~isa(np, /number, /scalar) then begin
   message, umsg, /inf
   message, 'NP should be the (complex) refractive index of the sphere', /inf
   return, -1
endif

if ~isa(nm, /number, /scalar) then begin
   message, umsg, /inf
   message, 'NM should be the (complex) refractive index of the medium', /inf
   return, -1
endif

if ~isa(lambda, /number, /scalar) then begin
   message, umsg, /inf
   message, 'LAMBDA should be the vacuum wavelength of the light in micrometers', /inf
   return, -1
endif

if ~isa(mpp, /number, /scalar) then begin
   message, umsg, /inf
   message, 'MPP should be the magnification in micrometers per pixel', /inf
   return, -1
endif

if ~isa(dim, /number, /array) then begin
   message, umsg, /inf
   message, 'DIM should be the dimensions of the desired image', /inf
   return, -1
endif

nx = float(dim[0])
ny = float(dim[1])
npts = nx * ny
x = findgen(nx) - (nx - 1.)/2. - float(rp[0])
y = findgen(1, ny) - (ny - 1.)/2. - float(rp[1])
x = rebin(x, nx, ny)
y = rebin(y, nx, ny)

zp = float(rp[2])

field = spherefield(x, y, zp, ap, np, nm, lambda, mpp, /cartesian, $
                    precision = precision, $
                    k = k, $
                    gpu = gpu)

if ~isa(alpha, /number, /scalar) then $
   alpha = 1.d

if ~isa(delta, /number, /scalar) then $
   delta = 0.d

fac = alpha * exp(dcomplex(0., -k*(zp + delta)))

field *= fac
field[0,*] += 1.d
a = total(real_part(field * conj(field)), 1)

return, reform(a, nx, ny)
end
