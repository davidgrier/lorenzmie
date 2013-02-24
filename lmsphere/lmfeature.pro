;+
; NAME:
;    lmfeature
;
; PURPOSE:
;    Identify, locate and characterize spheres in in-line
;    holographic video microscopy images using
;    Lorenz-Mie microscopy.
;
; CATEGORY:
;    Image analysis, holographic video microscopy, feature detection
;
; CALLING SEQUENCE:
;    features = lmfeature(hologram, lambda, mpp)
;
; INPUTS:
;    hologram: two-dimensional normalized holographic microscopy image.
;    lambda: vacuum wavelength of light [micrometers]
;    mpp: length calibration [micrometers per pixel]
;
; KEYWORD PARAMETERS:
;    Parameters for 2D feature location
;    noise: Estimate for the RMS additive noise at each pixel.
;        Default: estimated from image
;
;    pickn: Maximum number of features to consider.
;        Default: Process all features
;
;    Parameters for 3D feature location and characterization:
;    ap: ballpark radius of sphere [micrometers]
;        Default: estimated from image
;
;    np: ballpark complex refractive index of particle.
;        Default: 1.5
;
;    nm: complex refractive index of medium.
;        Default: water at room temperature for the given wavelength.
;
; KEYWORD FLAGS:
;    deinterlace: If set to an even number, consider only even
;        lines in feature extraction.  Similarly, if set to an
;        odd number.
;        NOTE: DEINTERLACE is not currently implemented.
;
;    gpu: If set, use hardware acceleration through GPULIB.
;        Requires appropriate hardware and installation of GPULIB.
;
;    graphics: If set, creates a graphical window showing progress of
;        analysis.  If set to a named variable, the same window will
;        be used in subsequent calls.
;
;    quiet: If set, minimize output
;
; OUTPUTS:
;    features: list() of parameters describing each of the features
;       in the hologram.  For the n-th feature:
;    p = features[n]:
;    p[0:1, 0]: [xp, dxp] : x location and error [pixels]
;    p[0:1, 1]: [yp, dyp] : y location and error [pixels]
;    p[0:1, 2]: [zp, dzp] : z location and error [pixels]
;    p[0:1, 3]: [ap, dap] : radius and error [micrometers]
;    p[0:1, 4]: [np, dnp] : refractive index and error of sphere
;    p[0:1, 5]: [kp, dnp] : extinction coefficient and error of sphere
;    p[0:1, 6]: [nm, dnm] : refractive index of medium
;    p[0:1, 7]: [km, dkm] : extinction coefficient of medium
;    p[0:1, 8]: [alpha, dalpha] : relative illumination amplitude
;    p[0:1, 9]: [delta, ddelta] : wavefront distortion
;    p[0,  10]: range of region of interest
;    p[1,  10]: chi-squared value of fit
;      
; PROCEDURE:
;   CT_FEATURE identifies candidate features using CIRCLETRANSFORM.
;   CT_RANGE estimates range around each feature.
;   RS1D estimates zp.
;   Poisson-spot model estimates ap and improves zp.
;   FITSPHERELM1D refines estimate for zp, ap and np.
;   FITSPHERELM yields final estimates for parameters.
;
; NOTES:
; 01/23/2013 DGG Tried scanning initial estimates for zp and ap in
;   fitlmsphere1d.  Fits consistently locked on to same solution over
;   entire range.  Z-scan not necessary, therefore.
;
; MODIFICATION HISTORY:
; 10/07/08: Written by FC Cheong and David G. Grier, New York University.
; 11/07/08: Modified by FC Cheong to include keyword FIT
; 12/07/08: Modified by FC Cheong to include Rayleigh Sommerfeld
;    approximation
; 01/27/2009 DGG Major overhaul.  Use fitspheredhm1d rather than
;    fresneldhm to estimate zp, ap and np.  Implemented DEINTERLACE
;    First draft of algorithm to filter against "bad" features.
; 02/08/2009 DGG Avoid fitting to badly estiamted features.
;    This prevents lock-ups in dhmtool.  Major documentation overhaul.
; 02/14/2009 DGG Added LIMITAP and LIMITNP keywords.
; 07/13/2012 DGG Major overhaul:
;    Added COMPILE_OPT.  Removed SMOOTH_FACTOR, limts variables and LUT.
;    Added PICKN, GRAPHICS.  Complete rewrite of variable estimation code.
;    GPU now uses object code for fitting.
; 07/18/2012 DGG Improved algorithm for estimating RAD using
;   CT_HITORMISS.
; 10/12/2012 DGG Updated for revisions to fitspheredhm.  LAMBDA and
;   MPP are now required inputs.
; 11/11/2012 DGG Treat nm and np as complex values.  Eliminate km and kp.
; 11/30/2012 DGG Updated for revisions to CTFEATURE and related
;   routines.  Implemented FIXALPHA.
; 01/16/2013 DGG Estimate rad using CT_RANGE.  Remove RAD keyword.
; 02/04/2013 DGG Added SMOOTH keyword to improve performance with
;   noisy images.  Will retire this once CT routines have better 
;   noise performance.
; 02/09/2013 DGG Fixed noise performance; retired SMOOTH.  Revised
;   coordinate code to locate features relative to lower-left corner
;   rather than relative to center of cropped image.  Removed
;   THRESHOLD keyword.  Features returned as list.  Return empty list
;   on failure.
; 02/14/2013 DGG Fixed deinterlace for cropped images.
; 02/15/2013 DGG Crop to odd dimensions.
; 02/16/2013 DGG Introduced lmf_report to clean up reporting.  Use
;   information about radius to estimate axial position.
; 02/24/2013 DGG Use LMF_RANGE to compute range.  Return RAD and CHISQ
;   with parameters for each feature.
;
; Copyright (c) 2008-2013 David G. Grier and Fook Chiong Cheong
;-
function lmfeature, a, lambda, mpp, $
                    noise = noise, $
                    pickn = pickn, $
                    ap = ap0, $
                    np = np0, $
                    nm = nm0, $
                    fixalpha = fixalpha, $
                    fixdelta = fixdelta, $
                    gpu = gpu, $
                    deinterlace = deinterlace, $
                    graphics = graphics, $
                    quiet = quiet, $
                    debug = debug

COMPILE_OPT IDL2

umsg = 'USAGE: p = lmfeature(a, lambda, mpp)'

features = list()

;;; Process command line parameters
if n_params() ne 3 then begin
   message, umsg, /inf
   return, features
endif

sz = size(a)
if ~isa(a, /number, /array) or sz[0] ne 2 then begin
   message, umsg, /inf
   message, 'A must be a two-dimensional normalized hologram', /inf
   return, features
endif
width = sz[1]
height = sz[2]

if ~isa(lambda, /number, /scalar) then begin
   message, umsg, /inf
   message, 'LAMBDA must be a number', /inf
   return, features
endif

if ~isa(mpp, /number, /scalar) then begin
   message, umsg, /inf
   message, 'MPP must be a number', /inf
   return, features
endif

if ~isa(np0, /number, /scalar) then $
   np0 = dcomplex(1.5d, 0.d)

if ~isa(nm0, /number, /scalar) then $
   nm0 = refractiveindex(lambda, 24.)

k = 2.d * !dpi * real_part(nm0) / lambda ; wavenumber [radians/um]

doap = ~isa(ap0, /number, /scalar)  ; estimate radius of each feature

fixalpha = keyword_set(fixalpha)
fixdelta = keyword_set(fixdelta)

gpu = keyword_set(gpu)
doreport = ~keyword_set(quiet)
dographics = arg_present(graphics) || keyword_set(graphics)
debug = keyword_set(debug)

;;; Find candidate features
rp = ctfeature(a, noise = noise, deinterlace = deinterlace, $
               pickn = pickn, count = count)

if doreport then begin
   print
   message, 'noise: ' + strtrim(noise, 2), /inf
   if keyword_set(deinterlace) then $
      message, 'analyzing ' + ((deinterlace mod 2) ? 'even' : 'odd') + ' field', /inf
   message, 'features found: ' + strtrim(count, 2), /inf
endif

if dographics then begin
   if isa(graphics, 'graphicswin') then $
      graphics.erase $
   else $
      graphics = window()
   im = image(a, axis_style = 2, /current)
   pf = plot(rp[0,*], rp[1,*], /over, $
             linestyle = '', symbol = 'o', color = 'red')
endif

;;; Loop over features to process each feature
for ndx = 0L, count - 1 do begin
   if doreport then message, 'feature ' + strtrim(ndx, 2), /inf

   ;; Initial estimate for xp and yp
   rp0 = rp[0:1, ndx]           ; particle position in original image

   ;; Crop image to region of interest
   thisrad = lmf_range(a, rp0, deinterlace = deinterlace)
   x0 = round(rp0[0] - thisrad) > 0L
   x1 = round(rp0[0] + thisrad + 1L) < width - 1L
   y0 = round(rp0[1] - thisrad) > 0L
   y1 = round(rp0[1] + thisrad + 1L) < height - 1L
   aa = a[x0:x1, y0:y1]         ; cropped image
   r0 = double([x0, y0])        ; lower-left corner
   rc = rp0 - r0                ; particle position in cropped image

   if dographics then begin
      poly = [[x0, y0], [x1, y0], [x1, y1], [x0, y1], [x0, y0]]
      roi = plot(poly, over = graphics, linestyle = '--', color = 'light green')
   endif

  ;; Use radial profile to estimate parameters
   b = aziavg(aa, center = rc, rad = thisrad, $
              deinterlace = deinterlace) ; radial profile around center

   ;; Model observed interference pattern as Poisson's spot to
   ;; relate axial position, zp, and sphere radius, ap
   c = fix(b ge 1)
   z0 = float(where(abs(c - c[1:*]) gt 0)) + 1. ; FIXME: zero crossings
   j0n = [2.4048, 5.5201, 8.6537]               ; zeros of J0(x)
   fac = 2.*k*mean(z0/j0n) - 1.
   if doap then begin
      ;; Estimate axial position as position of peak brightness in the
      ;; Rayleigh-Sommerfeld reconstruction: zp
      z = dindgen(50) * 4.d + 10.d ; FIXME: set range more intelligently
      res = rs1d(aa, z, rc, lambda = lambda/nm0, mpp = mpp)
      m = max(abs(res), loc)
      zp = z[loc]                       ; axial position [pixels]
      ;; Estimate radius using model
      ap = zp /fac ; estimated radius [um]
   endif else begin
      ;; Use provided radius
      ap = ap0
      ;; Estimate axial position using model
      zp = ap * fac
   endelse

   ;; Estimate np: FIXME: currently uses input value: np0

   ;; Estimate alpha and delta
   alpha = 1.d
   delta = 0.d

   p0 = [zp, ap, real_part(np0), imaginary(np0), $
         real_part(nm0), imaginary(nm0), alpha, delta] ; initial estimates

   if debug then lmf_report, 'starting estimates:', [rp0, p0]

   ;; Improve estimates by fitting to radial profile
   p1 = fitlmsphere1d(b, p0, lambda, mpp, fixdelta = fixdelta, chisq = thischisq, /quiet)
   if ~finite(thischisq) then begin
      message, 'parameter estimate failed -- trying again with fixed delta', /inf
      p1 = fitlmsphere1d(b, p0, lambda, mpp, /fixdelta, chisq = thischisq, /quiet)
   endif else begin
      ;; Some fits fail with alpha = 2.0; have to fixdelta.
      peggedalpha = (p1[0, 6] ge 1.9) && ~fixdelta
      if peggedalpha then begin
         message, 'estimate for alpha exceeds bounds -- trying again with fixed delta', /inf
         p1 = fitlmsphere1d(b, p0, lambda, mpp, /fixdelta, chisq = thischisq, /quiet)
      endif
   endelse
   if ~finite(thischisq) then begin
      message, 'parameter estimate failed -- skipping this feature', /inf
      continue
   endif
   
   if debug then lmf_report, 'improved estimates: ' + string(thischisq), $
                                [rp0, reform(p1[0, *])]

   if dographics then begin
      rp[0:1,ndx] = rp0
      pf.putdata, rp[0,*], rp[1,*]
   endif

   ;; 2D fit to refine estimates
   p2 = [rc, reform(p1[0, *])] ; initial parameters from 1D fit
   thisfeature = fitlmsphere(aa, p2, lambda, mpp, $
                             chisq = thischisq, $
                             fixalpha = fixalpha, $
                             fixdelta = fixdelta || peggedalpha, $
                             deinterlace = keyword_set(deinterlace) ? deinterlace + y0 : 0, $
                             object = gpu, $
                             quiet = ~debug)
   if n_elements(thisfeature) eq 1 then begin ; fit failed
      message, 'fit failed -- skipping this feature', /inf
      continue
   endif

   thisfeature[0,0:1] += r0     ; center in original image
   
   if doreport then lmf_report, 'chisq: ' + string(thischisq), thisfeature[0, *]

   if dographics then begin
      roi.setproperty, color = 'green'
      rp[0:1,ndx] = thisfeature[0,0:1]
      pf.putdata, rp[0,*], rp[1,*]
   endif

   features.add, [[thisfeature], [thisrad, thisrad]]
endfor

return, features
end
