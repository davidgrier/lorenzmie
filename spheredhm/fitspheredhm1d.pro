;+
; NAME:
;      fitspheredhm1d
;
; PURPOSE:
;      Measure the radius, refractive index, and axial
;      position of a colloidal sphere immersed in a dielectric 
;      medium by fitting its azimuthally averaged
;      digital holographic microscopy (DHM)
;      image to Mie scattering theory.
;
; CATEGORY:
;      Holographic microscopy
;
; CALLING SEQUENCE:
;      params = fitspheredhm1d(a,p)
;
; INPUTS:
;      a : two-dimensional real-valued DHM image of sphere.
;
;      p : initial guess for fitting parameters.
;      p[0] : zp : z-coordinate [pixels]
;      p[1] : ap : sphere's radius [micrometers]
;      p[2] : np : sphere's refractive index
;      p[3] : amplitude : arbitrary.  1 is a reasonable value
;      p[4] : nm : refractive index of medium
;
; KEYWORD PARAMETERS:
;      lambda: vacuum wavelength of illumination [micrometers].
;           Default: 0.632816
;      mpp: Length-scale calibration factor [micrometers/pixel].
;           Default: 0.135
;      precision: Accuracy with which scattering coefficients are
;           summed.  PRECISION=0.001 noticably speeds up 
;           CPU-based fits at the expense of some precision
;           in the fitting parameters.  
;           Default: PRECISION=0, retain all terms.
;
;      aplimits: [minap, maxap] limits on ap [micrometers]
;           Default: [0.05, 10.]
;
;      nplimits: [minnp, maxnp] limits on np
;           Default: [1.01 nm, 3.]
;
; KEYWORD FLAGS:
;      fixap: If set, do not allow ap to vary.
;      fixnp: If set, do not allow np to vary.
;      fixnm: If set, do not allow nm to vary.
;      fixzp: If set, do not allow zp to vary.
;      fixalpha: If set, do not allow alpha to vary.
;
;      gpu: If set, use GPU acceleration to calculate fields on
;           systems with GPUlib installed.  Requires NVIDIA graphics
;           accelerator with CUDA support.
;
;      quiet: If set, do not show results of intermediate calculations.
;
; OUTPUTS:
;      params: Least-squares fits for the values estimated in P.
;              params[0,*]: Fit values.
;              params[1,*]: Error estimates.
;              NOTE: errors are set to 0 for parameters held constant
;              with the FIX keywords.
;
; RESTRICTIONS:
;      Becomes slower and more sensitive to accuracy of initial
;      guesses as spheres become larger.
;
;      Refractive indices are restricted to be real-valued.
;
; PROCEDURE:
;      Uses MPFIT by Craig Marquardt to minimize the difference
;      between the measured DHM image and the image computed by
;      SPHEREDHM.
;
; REFERENCE:
;      S. Lee, Y. Roichman, G. Yi, S. Kim, S. Yang, A. van Blaaderen,
;      P. van Oostrum and D. G. Grier,
;      Chararacterizing and tracking
;      single colloidal particles with video holographic microscopy,
;      Optics Express 15, 18275-18282 (2007)
;
; MODIFICATION HISTORY:
; Written by David G. Grier, New York University, 01/16/2009.
; Based heavily on FITSPHEREDHM.
; 02/14/2009 DGG Added APLIMITS and NPLIMITS keywords.
; 07/15/2012 DGG Added COMPILE_OPT, code cleanups
;
; Copyright (c) 2009-2012 David G. Grier
;-

function spheredhm1d_f, rho, p, $
                        lambda = lambda, $
                        mpp = mpp, $
                        precision = precision, $
                        gpu = gpu

COMPILE_OPT IDL2, HIDDEN

; p[0] : zp
; p[1] : ap
; p[2] : np
; p[3] : amplitude
; p[4] : nm

zp = p[0]
ap = p[1]
np = p[2]
alpha = p[3]
nm = p[4]

field = spherefield(rho, fltarr(n_elements(rho)), zp, $
                    ap, np, $
                    k=k, nm=nm, lambda=lambda, mpp=mpp, /cartesian, $
                    precision=precision, gpu=gpu)

; interference between light scattered by the particle
; and a plane wave polarized along x and propagating along z
;dhm = 1.d + 2.d * alpha * field[0,*] * exp(dcomplex(0,-k*zp))
;dhm += alpha^2 * total(field * conj(field), 1)
;dhm = real_part(dhm)

field *= alpha * exp(dcomplex(0, -k*zp)) ; amplitude and phase factors
field[0, *] += 1.                        ; \hat{x}
dhm = total(real_part(field*conj(field)), 1)

return, dhm
end

function fitspheredhm1d, a, $               ; image
                         p0, $              ; starting estimates for parameters
                         lambda = lambda, $ ; wavelength of light [micrometers]
                         mpp = mpp, $       ; micrometers per pixel
                         aplimits = aplimits, $
                         nplimits = nplimits, $
                         fixnp = fixnp, $       ; fix particle refractive index
                         fixnm = fixnm, $       ; fix medium refractive index
                         fixap = fixap, $       ; fix particle radius
                         fixzp = fixzp, $       ; fix particle axial position
                         fixalpha = fixalpha, $ ; fix illumination
                         precision = precision, $ ; precision of coefficients
                         gpu = gpu, $             ; use gpu acceleration
                         quiet = quiet            ; don't print diagnostics

COMPILE_OPT IDL2

npts = n_elements(a)

err = replicate(5., npts)       ; FIXME this works but could be made rigorous

rho = dindgen(npts)             ; coordinates of pixels

if ~isa(lambda, /number, /scalar) then $
   lambda = 0.632816d           ; HeNe wavelength (micrometers)

if ~isa(mpp, /number, /scalar) then $
   mpp = 0.101                  ; Zeiss rig

if ~isa(precision, /number, scalar) then $
   precision = 0.               ; Keep all scattering coefficients

gpu = keyword_set(gpu)

nparams = n_elements(p0)
parinfo = replicate({limited:[0,0], limits:[0.,0.], fixed:0}, nparams)
; zp: No limit
; ap: Radius of particle.  Must be positive
parinfo[1].limited[0] = 1
parinfo[1].limits[0] = 0.05     ; smallest particle
parinfo[1].limited[1] = 1
parinfo[1].limits[1] = 10.      ; largest particle
if n_elements(aplimits) eq 2 then $
   parinfo[1].limits = aplimits
; np: Refractive index of particle
parinfo[2].limited[0] = 1
parinfo[2].limits[0] = 1.01*p0[4] ; FIXME what about low-index particles?
parinfo[2].limited[1] = 1
parinfo[2].limits[1] = 3.0      ; a bit more than titania
if n_elements(nplimits) eq 2 then $
   parinfo[2].limits = nplimits
; alpha: Illumination at particle
parinfo[3].limited[0] = 1
parinfo[3].limits[0] = 0.       ; cannot be negative
; nm: Refractive index of medium: No restrictions
; Flags to prevent adjusting parameter values
if keyword_set(fixzp)    then parinfo[0].fixed=1
if keyword_set(fixap)    then parinfo[1].fixed=1
if keyword_set(fixnp)    then parinfo[2].fixed=1
if keyword_set(fixalpha) then parinfo[3].fixed=1
if keyword_set(fixnm)    then parinfo[4].fixed=1

; parameters passed to the fitting function
argv = {lambda:lambda, mpp:mpp, precision:precision, gpu:gpu}

; errors from fit
perror = fltarr(nparams)

; perform fit
p = mpfitfun("spheredhm1d_f", rho, a, err, p0, functargs=argv, $
               parinfo = parinfo, /fastnorm, $
               perror=perror, bestnorm=bestnorm, dof=dof, quiet=quiet)


; failure?
if n_elements(p) eq 1 then return, -1

; success
; rescale fit uncertainties into error estimates
dp = perror * sqrt(bestnorm/dof)

return, [transpose(p),transpose(dp)]
end
