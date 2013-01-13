;+
; NAME:
;    fitlmsphere1d
;
; PURPOSE:
;    Measure the radius, refractive index, and axial
;    position of a colloidal sphere immersed in a dielectric 
;    medium by fitting its azimuthally averaged
;    holographic video microscopy (HVM)
;    image to Lorenz-Mie scattering theory.
;
; CATEGORY:
;    Holographic microscopy
;
; CALLING SEQUENCE:
;    params = fitlmsphere1d(a, p, lambda, mpp)
;
; INPUTS:
;    a : two-dimensional real-valued DHM image of sphere.
;
;    p : initial guess for fitting parameters.
;      p[0] : zp : z-coordinate [pixels]
;      p[1] : ap : sphere's radius [micrometers]
;      p[2] : np : sphere's refractive index
;      p[3] : kp : sphere's extinction coefficient
;      p[4] : nm : medium's refractive index of medium
;      p[5] : km : medium's extinction coefficient
;      p[6] : alpha : relative amplitude of illumination at particle
;      p[7] : delta : wavefront distortion of illumination
;
;    lambda: vacuum wavelength of illumination [micrometers].
;    mpp: Length-scale calibration factor [micrometers/pixel].
;
; KEYWORD PARAMETERS:
;    precision: Accuracy with which scattering coefficients are
;      summed.  PRECISION=0.001 noticably speeds up 
;      CPU-based fits at the expense of some precision
;      in the fitting parameters.  
;      Default: PRECISION=0, retain all terms.
;
;    aplimits: [minap, maxap] limits on ap [micrometers]
;      Default: [0.05, 10.]
;
;    nplimits: [minnp, maxnp] limits on np
;      Default: [1.01 nm, 3.]
;
; KEYWORD OUTPUTS:
;    fit: best fit to data
;
; KEYWORD FLAGS:
;    fixap: If set, do not allow ap to vary.
;    fixnp: If set, do not allow np to vary.
;    fixkp: If set, do not allow kp to vary.  Default: 1
;    fixnm: If set, do not allow nm to vary.  Default: 1
;    fixkm: If set, do not allow km to vary.  Default: 1
;    fixzp: If set, do not allow zp to vary.
;    fixalpha: If set, do not allow alpha to vary.
;    fixdelta: If set, do not allow delta to vary.
;
;    gpu: If set, use GPU acceleration to calculate fields on
;      systems with GPUlib installed.  Requires NVIDIA graphics
;      accelerator with CUDA support.
;
;    quiet: If set, do not show results of intermediate calculations.
;
; OUTPUTS:
;    params: Least-squares fits for the values estimated in P.
;      params[0,*]: Fit values.
;      params[1,*]: Error estimates.
;      NOTE: errors are set to 0 for parameters held constant
;        with the FIX keywords.
;
; RESTRICTIONS:
;    Becomes slower and more sensitive to accuracy of initial
;    guesses as spheres become larger.
;
; PROCEDURE:
;    Uses MPFIT by Craig Marquardt to minimize the difference
;    between the measured HVM image and the image computed by
;    LMSPHERE.
;
; REFERENCE:
;    S. Lee, Y. Roichman, G. Yi, S. Kim, S. Yang, A. van Blaaderen,
;    P. van Oostrum and D. G. Grier,
;    Chararacterizing and tracking
;    single colloidal particles with video holographic microscopy,
;    Optics Express 15, 18275-18282 (2007)
;
; MODIFICATION HISTORY:
; Written by David G. Grier, New York University, 01/16/2009.
; Based heavily on FITSPHEREDHM.
; 02/14/2009 DGG Added APLIMITS and NPLIMITS keywords.
; 07/15/2012 DGG Added COMPILE_OPT, code cleanups
; 10/15/2012 DGG Major revision to account for changes to spheredhm.
;    Renamed to fitlmsphere1d.
; 11/10/2012 DGG Added FIT keyword to return fitted function.
;    Constrain upper limit of alpha.  Update formatting documentation.
;
; Copyright (c) 2009-2012 David G. Grier
;-

function lmsphere1d_f, rho, p, $
                       lambda = lambda, $
                       mpp = mpp, $
                       precision = precision, $
                       gpu = gpu

COMPILE_OPT IDL2, HIDDEN

; p[0] : zp
; p[1] : ap
; p[2] : np
; p[3] : kp
; p[4] : nm
; p[5] : km
; p[6] : alpha
; p[7] : delta

zp = p[0]
ap = p[1]
np = dcomplex(p[2], p[3])
nm = dcomplex(p[4], p[5])
alpha = p[6]
delta = p[7]

field = spherefield(rho, fltarr(n_elements(rho)), zp, $
                    ap, np, nm, lambda, mpp, $
                    k = k, /cartesian, $
                    precision = precision, gpu = gpu)

; interference between light scattered by the particle
; and a plane wave polarized along x and propagating along z

field *= alpha * exp(dcomplex(0, -k*(zp + delta))) ; amplitude and phase factors
field[0, *] += 1.                        ; \hat{x}
dhm = total(real_part(field*conj(field)), 1)

return, dhm
end

function fitlmsphere1d, a, $                     ; image
                        p0, $                    ; starting estimates for parameters
                        lambda, $                ; wavelength of light [micrometers]
                        mpp, $                   ; micrometers per pixel
                        aplimits = aplimits, $
                        nplimits = nplimits, $
                        fixap = fixap, $          ; fix particle radius
                        fixnp = fixnp, $          ; fix particle refractive index
                        fixkp = fixkp, $          ; fix particle extinction coefficient
                        fixnm = fixnm, $          ; fix medium refractive index
                        fixkm = fixkm, $          ; fix medium extinction coefficient
                        fixzp = fixzp, $          ; fix particle axial position
                        fixalpha = fixalpha, $    ; fix illumination
                        fixdelta = fixdelta, $    ; fix wavefront distortion
                        precision = precision, $  ; precision of coefficients
                        fit = fit, $              ; best fit to data
                        gpu = gpu, $              ; use gpu acceleration
                        quiet = quiet             ; don't print diagnostics
  
COMPILE_OPT IDL2

umsg = 'USAGE: p = fitlmsphere1d(a, p, lambda, mpp)'

if n_params() ne 4 then begin
   message, umsg, /inf
   return, -1
endif

if ~isa(a, /number, /array) then begin
   message, umsg, /inf
   message, 'A must be a numerical array', /inf
   return, -1
endif

if ~isa(p0, /number, /array) then begin
   message, umsg, /inf
   message, 'P must be a numerical array', /inf
   return, -1
endif else if n_elements(p0) ne 8 then begin
   message, umsg, /inf
   message, 'P must have 8 elements', /inf
   return, -1
endif

if ~isa(lambda, /number, /scalar) then begin
   message, umsg, /inf
   message, 'LAMBDA must be a numerical scalar', /inf
   return, -1
endif

if ~isa(mpp, /number, /scalar) then begin
   message, umsg, /inf
   message, 'MPP must be a numerical scalar', /inf
   return, -1
endif

if ~isa(precision, /number, /scalar) then $
   precision = 0.               ; Keep all scattering coefficients

npts = n_elements(a)
err = replicate(0.1, npts)      ; FIXME this works but could be made rigorous
rho = dindgen(npts)             ; coordinates of pixels

gpu = keyword_set(gpu)

;;; Constraints on fitting parameters
nparams = n_elements(p0)
parinfo = replicate({limited:[0,0], $
                     limits:[0.,0.], $
                     fixed:0}, nparams)

;; Particle position
; zp: No constraints

;; Particle properties
; ap: Radius of particle.  Must be positive
parinfo[1].limited = 1
parinfo[1].limits = [0.05d, 10.d]
if n_elements(aplimits) eq 2 then $
   parinfo[1].limits = aplimits
parinfo[1].fixed = keyword_set(fixap)

; np: Refractive index of particle
parinfo[2].limited = 1
parinfo[2].limits = [1.01*p0[4], 3.d] ; FIXME what about low-index particles?
if n_elements(nplimits) eq 2 then $
   parinfo[2].limits = nplimits
; kp: Extinction coefficient of particle
parinfo[3].limited = 1
parinfo[3].limits = [0.d, 10.d]
parinfo[3].fixed = (n_elements(fixkp) eq 1) ? keyword_set(fixkp) : 1 ; Fixed by default

;; Medium properties
; nm: Refractive index of medium
parinfo[4].limited = 1
parinfo[4].limits = [1.d, 3.d]
parinfo[4].fixed = (n_elements(fixnm) eq 1) ? keyword_set(fixnm) : 1 ; Fixed by default
; km: Extinction coefficient of medium
parinfo[5].limited = 1
parinfo[5].limits = [0.d, 10.d]
parinfo[5].fixed = (n_elements(fixkm) eq 1) ? keyword_set(fixkm) : 1 ; Fixed by default

;; Illumination properties
; alpha: Illumination at particle
parinfo[6].limited = 1
parinfo[6].limits = [0.d, 2.d]
parinfo[6].fixed = keyword_set(fixalpha)
; delta:
parinfo[7].fixed = keyword_set(fixdelta)

; parameters passed to the fitting function
argv = {lambda:lambda, mpp:mpp, precision:precision, gpu:gpu}

; errors from fit
perror = fltarr(nparams)

; perform fit
p = mpfitfun('lmsphere1d_f', rho, a, err, p0, functargs = argv, $
             parinfo = parinfo, /fastnorm, $
             perror = perror, bestnorm = bestnorm, dof = dof, $
             yfit = fit, $
             quiet = quiet)

; failure?
if n_elements(p) eq 1 then return, -1

; success
; rescale fit uncertainties into error estimates
dp = perror * sqrt(bestnorm/dof)

return, [transpose(p),transpose(dp)]
end
