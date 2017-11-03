;+
; NAME:
;    lmspheres
;
; PURPOSE:
;    Computes Lorenz-Mie microscopy image of one or more spheres
;    immersed in a transparent medium.
;
; CATEGORY:
;    Holographic microscopy
;
; CALLING SEQUENCE:
;    holo = lmspheres(sample, instrument)
;
; INPUTS:
;    sample : dictionary of parameters describing a particle:
;      rp  : [x,y,z] 3 dimensional position of sphere relative to
;            lower-left corner of image.
;      ap  : radius of sphere [micrometers]
;      np  : (complex) refractive index of sphere
;      nm  : (complex) refractive index of medium
;    sample also can be an array of dictionaries
;
;    instrument : dictionary of parameters describing the instrument:
;      lambda: vacuum wavelength of light [micrometers]
;      mpp: micrometers per pixel
;      dim : [nx,ny] dimensions of image [pixels]
;
;    NOTE: The imaginary parts of the complex refractive indexes
;    should be positive for absorbing materials.  This follows the
;    convention used in SPHERE_COEFFICIENTS.
;
; OUTPUTS:
;    holo: [nx,ny] real-valued digital holographic image of spheres.
;
; KEYWORDS:
;    alpha: fraction of incident light scattered by particle.
;           Default: 1.
;
;    delta: wavefront distortion [wavelengths]
;           Default: 0.
;
;    field: On output, field in imaging plane
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
;    IDL> sample = dictionary()
;    IDL> sample.rp = [320,256,100]
;    IDL> sample.ap = 1.
;    IDL> sample.np = 1.5
;    IDL> sample.nm = 1.3
;    IDL> instrument = dictionary()
;    IDL> instrument.lambda = 0.532
;    IDL> instrument.mpp = 0.135
;    IDL> instrument.dim = [201,201]
;    IDL> holo1 = lmspheres(sample, instrument)
;
;    IDL> p = sample[*] ; make a copy
;    IDL> p.rp = [100,100,150]
;    IDL> sample2 = [sample, p]
;    IDL> holo2 = lmspheres(sample2, instrument)
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
;    Major change: Coordinates computed relative to lower-left corner,
;    rather than relative to center.  Simplifies interpretation of
;    coordinates and eliminates off-by-one and off-by-half errors.
; 03/06/2013 DGG Implementation of stratified spheres.  Added FIELD
;    keyword.
; 03/22/2013 DGG rebin(/sample) is more efficient.
; 10/03/2014 DGG simplify field calculation
; 11/02/2017 DGG introduce sphere parameters as a hash or an array of
;     hashes.
;
; Copyright (c) 2007-2017 David G. Grier
;-

function lmspheres, sample_, instrument, $
                    alpha = alpha, $
                    delta = delta, $
                    field = field, $
                    precision = precision, $
                    gpu = gpu

  COMPILE_OPT IDL2

  umsg = 'USAGE: lmspheres(sample, instrument)'

  if n_params() ne 2 then begin
     message, umsg, /inf
     return, -1
  endif

  sample = isa(sample_, 'dictionary') ? [sample_] : $
           isa(sample_, /array) ? sample_ : []
  if ~isa(sample[0], 'dictionary') then begin
     message, umsg, /inf
     message, 'SAMPLE should be a dictionary of parameters', /inf
     return, -1
  endif

  if ~isa(instrument, 'dictionary') then begin
     message, umsg, /inf
     message, 'INSTRUMENT should be a dictionary of parameters', /inf
     return, -1
  endif
  
  nx = float(instrument.dim[0])
  ny = float(instrument.dim[1])
  npts = nx * ny
  x = findgen(nx)
  y = findgen(1, ny)

  field = dcomplexarr(3, nx*ny)

  foreach p, sample do begin
     xp = rebin(x - p.rp[0], nx, ny, /sample)
     yp = rebin(y - p.rp[1], nx, ny, /sample)
     zp = float(p.rp[2])
     this = spherefield(xp, yp, zp, p.ap, p.np, p.nm, $
                        instrument.lambda, instrument.mpp, $
                        /cartesian, precision = precision, $
                        k = k)
     this *= exp(dcomplex(0., -k*zp)) ; scattered field
     field += this
  endforeach

  field[0,*] += 1.d                            ; incident field
  a = total(real_part(field * conj(field)), 1) ; intensity

  return, reform(a, nx, ny)
end
