;+
; NAME:
;    fitspheredhm
;
; PURPOSE:
;    Measure the radius, refractive index, and three-dimensional
;    position of a colloidal sphere immersed in a dielectric 
;    medium by fitting its digital holographic microscopy (DHM)
;    image to Mie scattering theory.
;
; CATEGORY:
;    Holographic microscopy
;
; CALLING SEQUENCE:
;    params = fitspheredhm(a, p)
;
; INPUTS:
;    a : two-dimensional real-valued DHM image of sphere.
;
;    p : initial guess for fitting parameters.
;      p[0] : xp : x-coordinate of sphere's center [pixels]
;      p[1] : yp : y-coordinate [pixels]
;      p[2] : zp : z-coordinate [pixels]
;      p[3] : ap : sphere's radius [micrometers]
;      p[4] : np : sphere's refractive index
;      p[5] : amplitude : arbitrary.  1 is a reasonable value
;      p[6] : nm : refractive index of medium
;
;      Optional:
;      p[7] : kp : sphere's extinction coefficient
;      p[8] : km : medium's extinction coefficient
;    NOTE: The extinction coefficients are assumed to be positive
;    which is appropriate for absorbing materials in the convention
;    used by SPHERE_COEFFICIENTS.
;
; KEYWORD PARAMETERS:
;    lambda: vacuum wavelength of illumination [micrometers].
;      Default: 0.632816
;
;    mpp: Length-scale calibration factor [micrometers/pixel].
;      Default: 0.135
;
;    precision: Convergence tolerance of nonlinear least-squares fit.
;      Default: 5d-5.
;
;    aplimits: [minap, maxap] limits on ap [micrometers]
;      Default: [0.05, 10.]
;
;    nplimits: [minnp, maxnp] limits on np
;      Default: [1.01 nm, 3.]
;
; KEYWORD FLAGS:
;    deinterlace: Only fit to odd (DEINTERLACE = 1) 
;      or even (DEINTERLACE = 2) scan lines.  This is useful for analyzing
;      holograms acquired with interlaced cameras.
;
;    fixap: If set, do not allow ap to vary.
;    fixnp: If set, do not allow np or kp to vary.
;    fixnm: If set, do not allow nm or km to vary.
;    fixzp: If set, do not allow zp to vary.
;    fixalpha: If set, do not allow alpha to vary.
;
;    gpu: If set, use GPU acceleration to calculate fields on
;         systems with GPULib installed.
;         Requires NVIDIA GPU with CUDA support.
;
;    object: If set, use a DGGdhmSphereDHM object to compute holograms.
;         This requires GPULib.
;
;    quiet: If set, do not show results of intermediate calculations.
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
;
; Copyright (c) 2007-2012, David G. Grier, Fook Chiong Cheong and
;    Paige Hasebe.
;-
function spheredhm_objf, obj, p

COMPILE_OPT IDL2, HIDDEN

; p[0] : xp         x position of sphere center
; p[1] : yp         y position of sphere center
; p[2] : zp         z position of sphere center
; p[3] : ap         radius of sphere
; p[4] : np         real part of sphere's refractive index
; p[5] : alpha  
; p[6] : nm         real part of medium's refractive index
;
; Optional:
; p[7] : kp         imaginary part of sphere's refractive index
; p[8] : km         imaginary part of medium's refractive index

nparams = n_elements(p)
obj.setproperty, rp = p[0:2], ap = p[3], alpha = p[5], $
                      np = dcomplex(p[4], (nparams ge 8) ? p[7] : 0), $
                      nm = dcomplex(p[6], (nparams ge 9) ? p[8] : 0)

return, obj.hologram
end

function spheredhm_f, x, y, p, $
                      lambda = lambda, $
                      mpp = mpp, $
                      gpu = gpu                      

COMPILE_OPT IDL2, HIDDEN

; p[0] : xp         x position of sphere center
; p[1] : yp         y position of sphere center
; p[2] : zp         z position of sphere center
; p[3] : ap         radius of sphere
; p[4] : np         real part of sphere's refractive index
; p[5] : alpha  
; p[6] : nm         real part of medium's refractive index
;
; Optional:
; p[7] : kp         imaginary part of sphere's refractive index
; p[8] : km         imaginary part of medium's refractive index

xx = x - p[0]
yy = y - p[1]
zp = p[2]
ap = p[3]
np = p[4]
alpha = p[5]
nm = p[6]
if n_elements(p) ge 8 then $
   np = dcomplex(np, p[7])      ; sphere's complex refractive index
if n_elements(p) ge 9 then $ 
   nm = dcomplex(nm, p[8])      ; medium's complex refractive index

field = spherefield(xx, yy, zp, ap, np, nm, lambda, mpp, $
                    /cartesian, $
                    gpu = gpu, $
                    k = k)

; interference between light scattered by the particle
; and a plane wave polarized along x and propagating along z
field *= alpha * exp(dcomplex(0, -k*zp)) ; amplitude and phase factors
field[0, *] += 1.                        ; \hat{x}

return, total(real_part(field * conj(field)), 1)
end

function fitspheredhm, a, $               ; image
                       p0, $              ; starting estimates for parameters
                       aplimits = aplimits, $ ; limits on ap [micrometers]
                       nplimits = nplimits, $ ; limits on np
                       lambda = lambda, $ ; wavelength of light [micrometers]
                       mpp = mpp, $       ; micrometers per pixel
                       fixnp = fixnp, $   ; fix particle refractive index
                       fixnm = fixnm, $   ; fix medium refractive index
                       fixap = fixap, $   ; fix particle radius
                       fixzp = fixzp, $   ; fix particle axial position
                       fixalpha = fixalpha, $ ; fix illumination
                       deinterlace = deinterlace, $
                       precision = precision, $ ; precision of convergence
                       gpu = gpu, $             ; use GPU acceleration
                       object = object, $       ; use DGGdhmSphereDHM object
                       quiet = quiet            ; don't print diagnostics

COMPILE_OPT IDL2

sz = size(a, /dimensions)
nx = sz[0]
ny = sz[1]
npts = nx*ny

if ~isa(lambda, /scalar, /number) then $
   lambda = 0.632816d           ; HeNe wavelength (micrometers)

if ~isa(mpp, /scalar, /number) then $
   mpp = 0.135d

if ~isa(precision, /scalar, /number) then $
   precision = 5d-5             ; Keep all scattering coefficients

if n_elements(gpu) ne 1 then $
   gpu = 0.

nparams = n_elements(p0)
parinfo = replicate({limited: [0,0], $
                     limits : [0.d, 0.d], $
                     fixed  : 0, $
                     step   : 0.}, nparams)
;; Restrictions on fitting parameters
; xp and yp: overly small steps sometimes prevent convergence
parinfo[0:2].step = 0.0001
; zp: No restrictions
; 
; ap: Radius must be positive
parinfo[3].limited[0] = 1
parinfo[3].limits[0] = 0.05
parinfo[3].limited[1] = 1
parinfo[3].limits[1] = 10.
if n_elements(aplimits) eq 2 then $
   parinfo[3].limits = aplimits
; np: Refractive index of particle
parinfo[4].limited[0] = 1
parinfo[4].limits[0] = 1.01*p0[6] ; FIXME what about low-index particles?
parinfo[4].limited[1] = 1
parinfo[4].limits[1] = 3.0      ; a bit more than titania
if n_elements(nplimits) eq 2 then $
   parinfo[4].limits = nplimits
; alpha: Illumination at particle
parinfo[5].limited[0] = 1
parinfo[5].limits[0] = 0.       ; cannot be negative
; nm: Refractive index of medium: No restrictions
; Flags to prevent parameters from being adjusted
parinfo[2].fixed = keyword_set(fixzp)
parinfo[3].fixed = keyword_set(fixap)
parinfo[4].fixed = keyword_set(fixnp)
parinfo[5].fixed = keyword_set(fixalpha)
parinfo[6].fixed = keyword_set(fixnm)
if nparams ge 8 then begin
   parinfo[7].fixed = keyword_set(fixnp)
   parinfo[7].limited[*] = 1
   parinfo[7].limits = [0.d, 5.d]
endif
if nparams ge 9 then begin
   parinfo[8].fixed = keyword_set(fixnm)
   parinfo[8].limited[*] = 1
   parinfo[8].limits = [0.d, 5.d]
endif

; errors from fit
perror = fltarr(nparams)

if keyword_set(object) then begin
   np = dcomplex(p0[4], (nparams ge 8) ? p0[7] : 0)
   nm = dcomplex(p0[6], (nparams ge 8) ? p0[8] : 0)
   obj = DGGdhmSphereDHM(dim = [nx, ny], $
                         lambda = lambda, $
                         mpp = mpp, $
                         rp = p0[0:2], $
                         ap = p0[3], $
                         nm = nm, $
                         np = np, $
                         alpha = p0[5], $
                         deinterlace = deinterlace $
                         )

   if ~isa(obj, 'DGGdhmSphereDHM') then begin
	message, 'could not create a DGGdhmSphereDHM object', /inf
	return, -1
   endif

   aa = double(a)
   if n_elements(deinterlace) gt 0 then begin
      w = where((lindgen(ny) mod 2) eq (deinterlace mod 2), ny)
      aa = aa[*, w]
   endif
   err = replicate(5., nx * ny) ; FIXME this works but could be made rigorous
      
   p = mpfitfun("spheredhm_objf", obj, aa, err, p0, $
                parinfo = parinfo, /fastnorm, $
                perror = perror, bestnorm = bestnorm, dof = dof, $
                status = status, errmsg = errmsg, quiet = quiet, $
                ftol = precision)
endif else begin
   x = dindgen(npts) mod nx     ; coordinates of pixels
   y = double(floor(dindgen(npts) / nx))

   aa = double(reform(a, npts))

   if n_elements(deinterlace) gt 0 then begin
      w = where((y mod 2) eq (deinterlace mod 2), npts)
      x = x[w]
      y = y[w]
      aa = aa[w]
   endif
   err = replicate(5., npts)    ; FIXME this works but could be made rigorous

   x -= double(nx - 1.d) / 2.d
   y -= double(ny - 1.d) / 2.d

; parameters passed to the fitting function
   argv = {lambda:lambda, mpp:mpp, precision:precision, gpu:gpu}

; perform fit
   p = mpfit2dfun("spheredhm_f", x, y, aa, err, p0, functargs = argv, $
                  parinfo = parinfo, /fastnorm, $
                  perror = perror, bestnorm = bestnorm, dof = dof, $
                  status = status, errmsg = errmsg, quiet = quiet, $
                  ftol = precision)

endelse

if status le 0 then begin 
   message, errmsg, /inf
   return, -1
endif

; failure?
if n_elements(p) eq 1 then begin
   message, "MPFIT2DFUN did not return a result",  /inf
   return, -1
endif

; success
; rescale fit uncertainties into error estimates
dp = perror * sqrt(bestnorm/dof)

return, [transpose(p),transpose(dp)]
end
