;+
; NAME:
;    spheredhm
;
; PURPOSE:
;    Computes holographic microscopy image of a sphere
;    immersed in a transparent medium.
;
; CATEGORY:
;    Holographic microscopy
;
; CALLING SEQUENCE:
;    holo = spheredhm(rp, ap, np, nm, dim)
;
; INPUTS:
;    rp  : [x,y,z] 3 dimensional position of sphere relative to center
;          of image.
;    ap  : radius of sphere [micrometers]
;    np  : (complex) refractive index of sphere
;    nm  : (complex) refractive index of medium
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
;    lambda: vacuum wavelength of light [micrometers]
;
;    mpp: micrometers per pixel
;
;    precision: relative precision with which fields are calculated.
;
; KEYWORD FLAGS:
;    gpu: Use GPU accelerated calculation.  This only works on
;         systems with GPUlib installed.  Requires NVIDIA graphics
;         card with CUDA installed.
;
;    lut: interpolate two-dimensional result from one-dimensional 
;         look-up table, with _substantial_ speed benefits, at the
;         cost of some precision.
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

; EXAMPLE:
;    Display a DHM image of a 1.5 micrometer diameter polystyrene
;    sphere (np = 1.5) in water (nm = 1.33).
;
;    IDL> tvscl, spheredhm([0,0,200],0.75,1.5,1.33,[201,201])
;
; MODIFICATION HISTORY:
;  Written by David G. Grier, New York University, 3/2007
;  05/25/2007 DGG: Added ALPHA keyword.
;  02/08/2008 DGG: Adopted updated SPHEREFIELD syntax:
;             separate x, y, and z coordinates.
;  10/09/2008 DGG: Added LAMBDA and PRECISION keywords
;  10/14/2008 DGG: Added GPU keyword.
;  10/16/2008 DGG: Added LUT keyword.
;  03/15/2010 DGG: Documentation cleanups.
;  11/03/2010 DGG: ALPHA defaults to 1.  Clean up image calculation.
;    Documentation fixes and formatting.
;  06/23/2012 DGG flipped sign of xp and yp for consistency with
;    fitspheredhm and dggdhmspheredhm.
;
; Copyright (c) 2007-2012 David G. Grier
;-

function spheredhm, rp, ap, np, nm, dim, $
                    alpha = alpha, $
                    mpp = mpp, $
                    lambda = lambda, $
                    precision = precision, $
                    gpu = gpu, $
                    lut = lut

; rp: [xp,yp,zp] position of sphere relative to center of image [pixels]
;          NOTE: We take z to denote the sphere's
;                height above the imaging plane.  In general,
;                therefore, z should be input as a positive number.
; ap: radius of sphere [microns]
; np: (complex) index of particle
; nm: (complex) index of medium
; dim: [nx,ny] dimensions of image to create [pixels]

nx = float(dim[0])
ny = float(dim[1])
npts = nx * ny
x = findgen(nx) - (nx - 1.)/2. - float(rp[0])
y = findgen(1, ny) - (ny - 1.)/2. - float(rp[1])
x = rebin(x, nx, ny)
y = rebin(y, nx, ny)

zp = float(rp[2])

if keyword_set(lut) then begin
    rho = sqrt(x^2 + y^2)
    x = findgen(1, fix(max(rho))+1)
    y = 0. * x
endif

field = spherefield(x, y, zp, ap, np, nm = nm, /cartesian, $
                    mpp = mpp, lambda = lambda, $
                    precision = precision, $
                    k = k, $
                    gpu = gpu)

fac = exp(dcomplex(0., -k*zp))
if isa(alpha, /number, /scalar) then $
   fac *= alpha
field *= fac
field[0,*] += 1.d
dhm = total(real_part(field * conj(field)), 1)
; scattered intensity
;dhm = total(real_part(field * conj(field)), 1)
; interference between scattered wave and incident plane wave
;dhm += 2.d * real_part(field[0,*] * exp(dcomplex(0,-k*zp)))
; intensity of plane wave
;dhm += 1.d

if keyword_set(lut) then $
  dhm = interpolate(dhm, rho, cubic=-0.5)

return, reform(dhm, nx, ny)
end
