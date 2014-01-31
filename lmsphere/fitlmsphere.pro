;+
; NAME:
;    fitlmsphere
;
; PURPOSE:
;    Measure the radius, refractive index, and three-dimensional
;    position of a colloidal sphere immersed in a dielectric 
;    medium by fitting its holographic video microscopy (HVM)
;    image to predictions of Lorenz-Mie scattering theory.
;
; CATEGORY:
;    Holographic microscopy
;
; CALLING SEQUENCE:
;    params = fitlmsphere(a, p, lambda, mpp)
;
; INPUTS:
;    a : two-dimensional real-valued DHM image of sphere.
;
;    p : initial guess for fitting parameters.
;      p[0] : xp : x-coordinate of sphere's center [pixels]
;      p[1] : yp : y-coordinate [pixels]
;      p[2] : zp : z-coordinate [pixels]
;      p[3] : ap : sphere's radius [meicrometers]
;      p[4] : np : sphere's refractive index
;      p[5] : kp : sphere's extinction coefficient: Default: fixed at 0.
;      p[6] : nm : medium's refractive index
;      p[7] : km : medium's extinction coefficient: Default: fixed at 0.
;      p[8] : alpha: relative amplitude of illumination at
;             sphere's position, typically around 1.
;      p[9] : delta: wavefront distortion of illumination at
;             sphere's position [pixels]
;
;    NOTE: The extinction coefficients are assumed to be non-negative
;    which is appropriate for absorbing materials in the convention
;    used by SPHERE_COEFFICIENTS.
;
;    lambda: vacuum wavelength of illumination [micrometers].
;    mpp: Length-scale calibration factor [micrometers/pixel].
;
; KEYWORD INPUTS:
;    errors: Estimate for error at each pixel.  Either an array with
;      the same dimensions as A, or a single digit for an overall
;      error estimate.
;      Default: 1.D
;
;    precision: Convergence tolerance of nonlinear least-squares fit.
;      Default: 5d-5.
;
;    aplimits: [minap, maxap] limits on ap [micrometers]
;      Default: [0.05, 10.]
;
;    nplimits: [minnp, maxnp] limits on np
;      Default: [1.01, 3.] relative to nm
;
;    resolution: Set resolution to which Lorenz-Mie coefficients
;      should be calculated.
;      Default: 0 -- retain full resolution.
;      Setting resolution = 1e-5 speeds the calculation substantially
;      without appreciably influencing results.
;
; KEYWORD FLAGS:
;    deinterlace: Only fit to odd (DEINTERLACE = 1) 
;      or even (DEINTERLACE = 2) scan lines.  This is useful for analyzing
;      holograms acquired with interlaced cameras.  Input hologram and
;      optional error estimates are assumed to be deinterlaced.  Setting this
;      keyword then ensures that the y coordinates are scaled and
;      offset appropriately.
;
;    fixap: If set, do not allow ap to vary.
;    fixnp: If set, do not allow np to vary.
;    fixkp: If set, do not allow kp to vary: Set to 1 by default
;    fixnm: If set, do not allow nm to vary: Set to 1 by default
;    fixkm: If set, do not allow km to vary: Set to 1 by default
;    fixzp: If set, do not allow zp to vary.
;    fixalpha: If set, do not allow alpha to vary.
;    fixdelta: If set, do not allow delta to vary.
;
;    gpu: If set, use GPU acceleration to calculate fields on
;         systems with GPULib installed.
;         Requires NVIDIA GPU with CUDA support.
;
;    quiet: If set, do not show results of intermediate calculations.
;
; KEYWORD OUTPUTS:
;    chisq: Chi-squared value of a (successful) fit.
;
;    yfit:  Fit obtained from best parameters
;
;    residuals: Residuals from the best fit.
;
; OUTPUTS:
;    params: Least-squares fits for the values estimated in P.
;      params[0,*]: Fit values.
;      params[1,*]: Error estimates.
;      NOTE: errors are set to 0 for parameters held constant
;            with the FIX* keyword flags.
;
; RESTRICTIONS:
;    Becomes slower and more sensitive to accuracy of initial
;    guesses as spheres become larger.
;
; PROCEDURE:
;    Uses MPFIT by Craig Marquardt (http://purl.com/net/mpfit/)
;    to minimize the difference between the measured DHM image and 
;    the image computed by SPHEREDHM.
;
; REFERENCES:
; 1. S. Lee, Y. Roichman, G. Yi, S. Kim, S. Yang, A. van Blaaderen,
;    P. van Oostrum and D. G. Grier,
;    Chararacterizing and tracking single colloidal particles with 
;    video holographic microscopy,
;    Optics Express 15, 18275-18282 (2007)
;
; 2. C. B. Markwardt,
;    Non-linear least squares fitting in IDL with MPFIT,
;    in Astronomical Data Analysis and Systems XVIII,
;    D. Bohlender, P. Dowler and D. Durand, eds.
;    (Astronomical Society of the Pacific, San Francisco, 2008).
;
; MODIFICATION HISTORY:
; Written by David G. Grier, New York University, 4/2007.
; 05/22/2007 DGG Added LAMBDA keyword.
; 05/26/2007 DGG Revised to use Bohren and Huffman version of
;   SPHEREFIELD.
; 06/10/2007 DGG Updated for more accurate BH code.
; 09/11/2007 DGG Made nm a fitting parameter and removed NM keyword.  
;   Replaced FIXINDEX keword with FIXNM.
;   Added FIXNP keyword.  
; 11/03/2007 DGG Changed FIXRADIUS to FIXAP.  Added FIXZP and FIXALPHA.
; 02/08/2008 DGG Treat coordinates as one-dimensional arrays internally
;   to eliminate repeated calls to REFORM.
;   Adopt updated syntax for SPHEREFIELD: separate x, y and z coordinates.
;   Y coordinates were incorrectly cast to float rather than double.
; 02/10/2008 DGG Added DEINTERLACE. Small documentation fixes.
; 04/16/2008 DGG Added MPP keyword.  Small documentation fixes.
; 10/13/2008 DGG Added PRECISION and GPU keywords to make use of new
;   capabilities in SPHEREFIELD.
; 10/17/2008 DGG Added LUT keyword to accelerate CPU-based fits.  
;   This required setting .STEP = 0.0001 pixel restrictions on
;   the x and y centroids in PARINFO.
; 01/15/2009 DGG Documentation clean-ups.
; 02/14/2009 DGG Added APLIMITS and NPLIMITS keywords.
; 03/17/2009 DGG Added support for complex refractive indexes by
;    accounting for the particle and medium extinction coefficients,
;    kp and km.
; 03/26/2009 Fook Chiong Cheong, NYU: np and nm should be cast as
;    dcomplex rather than complex when kp or km are non-zero.
; 06/18/2010 DGG Added COMPILE_OPT.
; 10/20/2010 DGG Cleaned up alpha code in spheredhm_f.
; 11/30/2010 DGG & FCC Vary kp and km independently.  Documentation
;    and formatting.
; 11/08/2011 DGG Removed LUT option: could not guarantee precision.
;    Added OBJECT keyword to compute fits with DGGdhmSphereDHM object
;    for improved efficiency.  Documentation upgrades.
; 11/09/2011 DGG PRECISION keyword now corresponds to FTOL in MPFIT.
; 04/17/2012 DGG Fixed deinterlace code for object-based fits for
;    centers aligned with grid but outside of field of view.
; 05/03/2012 DGG Updated parameter checking.
; 07/16/2012 DGG Used Paige Hasebe's faster approach to calculating
;    intensity in spheredhm_f. 
; 10/12/2012 DGG Major overhaul of parameter handling to incorporate
;    DELTA.  LAMBDA and MPP now are required inputs, rather than
;    optional keywords.  Renamed to fitlmsphere.
; 01/14/2013 DGG Added CHISQ keyword.
; 01/26/2013 DGG Compute coordinates relative to lower-left corner of
;    image rather than center.  Correctly ignore deinterlace = 0.
; 02/24/2013 DGG sample ERR when deinterlacing.
; 03/13/2013 DGG correct deinterlace code for object fitting.
; 03/22/2013 DGG rebin(/sample) is more efficient.
; 05/28/2013 DGG optionally return residuals from fit.
; 05/29/2013 DGG optionally return fit itself.  Honor FIXNP flag.
; 06/04/2013 DGG use ERRORS rather than WEIGHTS.  Document ERRORS
;    keyword.  Overhaul sampling code when called with deinterlace.
; 07/19/2013 DGG object now handles CPU code.  Use object throughout.
;    Remove OBJECT keyword.
; 07/24/2013 DGG Added RESOLUTION keyword.
; 08/11/2013 DGG Simplify test for failed fit.
;
; Copyright (c) 2007-2013, David G. Grier, Fook Chiong Cheong and
;    Paige Hasebe.
;-
function lmsphere_f, obj, p

COMPILE_OPT IDL2, HIDDEN

; p[0] : xp         x position of sphere center
; p[1] : yp         y position of sphere center
; p[2] : zp         z position of sphere center
; p[3] : ap         radius of sphere
; p[4] : np         real part of sphere's refractive index
; p[5] : kp         imaginary part of sphere's refractive index
; p[6] : nm         real part of medium's refractive index
; p[7] : km         imaginary part of medium's refractive index
; p[8] : alpha      amplitude of illumination at particle's position
; p[9] : delta      wavefront distortion at particle's position

obj.setproperty, rp = p[0:2], ap = p[3], np = dcomplex(p[4], p[5]), $
                 nm = dcomplex(p[6], p[7]), $
                 alpha = p[8], delta = p[9]  

return, obj.hologram
end

function fitlmsphere, a, $                     ; image
                      p0, $                    ; starting estimates for parameters
                      lambda, $                ; wavelength of light [micrometers]
                      mpp, $                   ; micrometers per pixel
                      errors = errors, $       ; error estimate at each pixel
                      precision = precision, $ ; precision of convergence
                      chisq = chisq, $         ; chi-squared value of fit
                      yfit = yfit, $           ; best fit to data
                      residuals = residuals, $ ; residuals from fit
                      aplimits = aplimits, $   ; limits on ap [micrometers]
                      nplimits = nplimits, $   ; limits on np
                      fixnp = fixnp, $         ; fix particle refractive index
                      fixkp = fixkp, $         ; fix particle extinction coefficient
                      fixnm = fixnm, $         ; fix medium refractive index
                      fixkm = fixkm, $         ; fix medium extinction coefficient
                      fixap = fixap, $         ; fix particle radius
                      fixzp = fixzp, $         ; fix particle axial position
                      fixalpha = fixalpha, $   ; fix illumination amplitude
                      fixdelta = fixdelta, $   ; fix wavefront distortion
                      deinterlace = deinterlace, $
                      resolution = resolution, $
                      gpu = gpu, $             ; use GPU acceleration
                      quiet = quiet            ; don't print diagnostics

COMPILE_OPT IDL2

;;; Command-line arguments
umsg = 'USAGE: p = fitslmphere(a, p0, lambda, mpp)'

if n_params() ne 4 then begin
   message, umsg, /inf
   return, -1
endif

if ~isa(a, /number, /array) then begin
   message, umsg, /inf
   message, 'A must be a two-dimensional numerical array', /inf
   return, -1
endif

if ~isa(lambda, /scalar, /number) then begin
   message, umsg, /inf
   message, 'LAMBDA must be a real-valued number', /inf
   return, -1
endif

if ~isa(mpp, /scalar, /number) then begin
   message, umsg, /inf
   message, 'MPP must be a real-valued number', /inf
   return, -1
endif

if ~isa(precision, /scalar, /number) then $
   precision = 5d-5             ; Keep all scattering coefficients

if n_elements(gpu) ne 1 then $
   gpu = 0.

sz = size(a, /dimensions)
nx = sz[0]
ny = sz[1]
npts = nx*ny

aa = double(a)
err = isa(errors, /number) ? double(errors) : 1.D

;;; Constraints on fitting parameters
nparams = n_elements(p0)
parinfo = replicate({limited: [0, 0], $
                     limits : [0.d, 0.d], $
                     fixed  : 0, $
                     step   : 0.d}, nparams)

;; Particle position
; xp, yp: 
parinfo[0:1].step = 1d-4 ; overly small steps sometimes prevent convergence
if keyword_set(deinterlace) then $
   parinfo[1] *= 2.
; zp: 
parinfo[2].step = 1d-4
parinfo[2].fixed = keyword_set(fixzp)

;; Particle properties
; ap: particle radius
parinfo[3].limited = 1
parinfo[3].limits = [0.05d, 10.d]
if n_elements(aplimits) eq 2 then $
   parinfo[3].limits = aplimits
parinfo[3].fixed = keyword_set(fixap)

; np: Refractive index of particle
parinfo[4].limited = 1
parinfo[4].limits = [1.001d*p0[6], 3.d] ; FIXME what about low-index particles?
if n_elements(nplimits) eq 2 then $
   parinfo[4].limits = nplimits
parinfo[4].fixed = keyword_set(fixnp)
; kp: Extinction coefficient of particle
parinfo[5].limited = 1
parinfo[5].limits = [0.d, 10.d]
parinfo[5].fixed = (n_elements(fixkp) eq 1) ? keyword_set(fixkp) : 1 ; Fixed by default

;; Medium properties
; nm: Refractive index of medium
parinfo[6].limited = 1
parinfo[6].limits = [1.d, 3.d]
parinfo[6].fixed = (n_elements(fixnm) eq 1) ? keyword_set(fixnm) : 1 ; Fixed by default
; km: Refractive index of medium
parinfo[7].limited = 1
parinfo[7].limits = [0.d, 10.d]
parinfo[7].fixed = (n_elements(fixkm) eq 1) ? keyword_set(fixkm) : 1 ; Fixed by default

;; Illumination properties
; alpha: Illumination at particle
parinfo[8].limited = 1
parinfo[8].limits = [0.d, 2.d]
parinfo[8].fixed = keyword_set(fixalpha)
; delta: 
parinfo[9].fixed = keyword_set(fixdelta)

; errors from fit
perror = fltarr(nparams)

obj = DGGdhmLMSphere(dim = [nx, ny], $
                     lambda = lambda, $
                     mpp = mpp, $
                     rp = p0[0:2], $
                     ap = p0[3], $
                     nm = dcomplex(p0[4], p0[5]), $
                     np = dcomplex(p0[6], p0[7]), $
                     alpha = p0[8], $
                     delta = p0[9], $
                     deinterlace = deinterlace, $
                     resolution = resolution, $
                     gpu = gpu $
                    )

if ~isa(obj, 'DGGdhmLMSphere') then begin
   message, 'could not create a DGGdhmLMSphere object', /inf
   return, -1
endif

p = mpfitfun('lmsphere_f', obj, aa, err, p0, $
             ftol = precision, $
             parinfo = parinfo, /fastnorm, $
             perror = perror, bestnorm = chisq, dof = dof, $
             status = status, errmsg = errmsg, quiet = quiet, $
             yfit = yfit, best_resid = residuals)
residuals = reform(residuals, nx, n_elements(residuals)/nx)

; failure
if status le 0 then $
   return, -1

; success
; rescale fit uncertainties into error estimates
dp = perror * sqrt(chisq/dof)

return, [transpose(p), transpose(dp)]
end
