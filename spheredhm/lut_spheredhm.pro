;+
; NAME:
;      spheredhm
;
; PURPOSE:
;      Computes holographic microscopy image of a dielectric sphere
;      immersed in a dielectric medium.
;
; CATEGORY:
;      Holographic microscopy
;
; CALLING SEQUENCE:
;      holo = spheredhm(r,a,np,nm,dim)
;
; INPUTS:
;      r : [x,y,z] 3 dimensional position of sphere relative to center
;          of image.
;      a : radius of sphere [micrometers]
;      np : refractive index of sphere
;      nm : refractive index of medium
;      dim : [nx,ny] dimensions of image [pixels]
;
; OUTPUTS:
;      holo: nx X ny real-valued digital holographic image of sphere.
;
; KEYWORDS:
;      alpha: fraction of incident light scattered by particle.
;
;      lambda: vacuum wavelength of light [micrometers]
;
;      mpp: micrometers per pixel
;
;      precision: relative precision with which fields are calculated.
;
;      lut: interpolate two-dimensional result from one-dimensional 
;           look-up table, with _substantial_ speed benefits, at the
;           cost of some precision.
;           Value of lut sets sub-pixel resolution for interpolation.
;           Example: LUT = 5 gives 1/5 pixel resolution and greatly
;           reduces interpolation errors.
;
; KEYWORD FLAGS:
;      gpu: Use GPU accelerated calculation.  This only works on
;           systems with GPUlib installed.  Requires NVIDIA graphics
;           card with CUDA installed.
;
; PROCEDURE:
;      Calls SPHEREFIELD to compute the field.
;
;      Described in S.H.Lee, et al., Optics Express (2007).
;
; EXAMPLE:
;      Display a DHM image of a 1.5 micrometer diameter polystyrene
;      sphere (np = 1.5) in water (nm = 1.33).
;
;      IDL> tvscl, spheredhm([0,0,200],0.75,1.5,1.33,[201,201])
;
; MODIFICATION HISTORY:
;      Written by David G. Grier, New York University, 3/2007
;      DGG, 5/25/2007: Added ALPHA keyword.
;      DGG, 2/8/2008: Adopted updated SPHEREFIELD syntax:
;           separate x, y, and z coordinates.
;      DGG, 10/09/2008: Added LAMBDA and PRECISION keywords
;      DGG, 10/14/2008: Added GPU keyword.
;      DGG, 10/16/2008: Added LUT keyword.
;-

function spheredhm, rsphere, a, np, nm, dim, alpha = alpha, $
                    mpp=mpp, lambda=lambda, precision=precision, $
                    gpu=gpu, lut=lut

; rsphere: [x,y,z] position of sphere relative to center of image [pixels]
;          NOTE: We take z to denote the sphere's
;                height above the imaging plane.  In general,
;                therefore, z should be input as a positive number.
; a: radius of sphere [microns]
; np: (complex) index of particle
; nm: (complex) index of medium
; dim: [nx,ny] dimensions of image to create [pixels]

nx = float(dim[0])
ny = float(dim[1])
npts = nx * ny
x = findgen(nx) - nx/2. - float(rsphere[0])
x = x # replicate(1.,ny)
y = findgen(ny) - float(ny)/2. - float(rsphere[1])
y = y ## replicate(1., nx)

zp = float(rsphere[2])

if keyword_set(lut) then begin
    zoomfactor = float(lut)
    rho = sqrt(x^2 + y^2)
    x = findgen(1,zoomfactor * fix(max(rho))+1) / zoomfactor
    y = 0. * x
endif

field = spherefield(x, y, zp, a, np, nm=nm, $
                  mpp=mpp, lambda=lambda, precision=precision, $
                  k=k, /cartesian, gpu=gpu)

; interference between light scattered by particle
; and a plane wave propagating along z
dhm = field[0,*] * exp(dcomplex(0,-k*zp))

if n_elements(alpha) eq 1 then begin
    fsq = total(field*conj(field), 1)
    dhm = 1.d + 2.d*alpha * temporary(dhm) + alpha^2 * fsq
endif

dhm = real_part(dhm)

if keyword_set(lut) then $
  dhm = interpolate(dhm, zoomfactor * rho, cubic=-0.5)

return, reform(dhm,nx,ny)
end
