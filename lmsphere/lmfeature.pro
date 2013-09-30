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
;    nfringes: Number of fringes to use to set range.
;        Default: 20
;
;    maxrefits: Maxmimum number of refits due to large displacements.
;        Default: 1
;
;    resolution: Resolution to which Lorenz-Mie coefficients are to be
;        computed.  Should be a value between 0 and 1.
;        Default: 0 -- retain full resolution.
;        Setting resolution = 1e-5 speeds up computations enormously,
;        without appreciably influencing results.
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
; KEYWORD OUTPUTS:
;    residuals: Residuals of best fit.
;
;    count: number of features actually returned
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
; 03/12/2013 DGG FIXME: CTFEATURE returns too many small features,
;   necessitating use of PICKN to achieve stable results.
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
; 03/12/2013 DGG Replace LMF_RANGE with LMF_CROP, which also
;   calculates statistical weights. Incorporated LMF_REPORT, LMF_CROP.
;   Refit if initial (x,y) estimate is too far off the mark.
; 03/22/2013 DGG rebin(/sample) is more efficient.
; 05/06/2013 DGG and David Ruffner: Faster and more accurate estimate
;   of zp using spherical-wave model rather than back propagation.
; 05/13/2013 DGG Fixed 'odd'/'even' reporting for deinterlaced images.
;   Fixed cropping bug.
; 06/01/2013 DGG Identify fringes as extrema in the azimuthal average.
;   Use fringes to set range for fits.
; 06/02/2013 DGG Use deviates from azimuthal average to compute
;   weightings for fits.  Added NFRINGES and MAXREFITS keywords.
; 06/04/2013 DGG Added COUNT keyword.  Estiamte per-pixel errors for
;   two-dimensional fit.
; 07/19/2013 DGG Update for CPU object implementation.
; 07/25/2013 DGG Added RESOLUTION keyword.
; 08/09/2013 DGG Return array rather than LIST() for compatibility
;   with IDL_IDLBridge.
; 09/24/2013 DGG Use provided ap to estimate zp.
; 09/30/2013 DGG and Bhaskar Jyoti Krishnatreya: J0(x) estimate for ap.
;
; Copyright (c) 2008-2013 David G. Grier, Bhaskar Jyoti Krishnatreya,
;    David Ruffner and Fook Chiong Cheong
;-
;;;;;
;
; LMF_REPORT
;
; Print information about fit parameters
;
pro lmf_report, msg, p

COMPILE_OPT IDL2, HIDDEN

message, msg, /inf
message, string(p[0:2], $
                format = '("  rp = (",F0.2,", ",F0.2,", ",F0.2,")")'), /inf
message, string(p[3:4], $
                format = '("  ap = ",F0.3," um, np = ",F0.3)'), /inf
message, string(p[8:9], $
                format = '("  alpha = ",F0.3,", delta = ",F0.3)'), /inf
message, /inf

return
end

;;;;;
;
; LMFEATURE
;
; Main routine
;
function lmfeature, a, lambda, mpp, $
                    noise = noise, $
                    maxrefits = maxrefits, $
                    nfringes = nfringes, $
                    pickn = pickn, $
                    count = count, $
                    ap = ap, $
                    fixap = fixap, $
                    np = np0, $
                    fixnp = fixnp, $
                    nm = nm0, $
                    fixalpha = fixalpha, $
                    fixdelta = fixdelta, $
                    roi = roi, $
                    residuals = residuals, $
                    resolution = resolution, $
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
maxx = sz[1]-1
maxy = sz[2]-1

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

dozp = 1                                          ; estimate axial position of each feature
doap = ~isa(ap, /number, /scalar)                 ; estimate radius of each feature
j0n = [2.4048, 5.5201, 8.6537, 11.7915, 14.9309] ; zeros of J0(x) used for estimating ap
;j1n = [3.8317, 7.0156, 10.1735, 13.3237, 16.4706] ; zeros of J1(x) used for estimating ap

fixalpha = keyword_set(fixalpha)
fixdelta = keyword_set(fixdelta)

gpu = keyword_set(gpu)
doreport = ~keyword_set(quiet)
dographics = arg_present(graphics) || keyword_set(graphics)
if dographics then begin
   if isa(graphics, 'lmfgraphics') then $
      graphics.im -> putdata, deinterlace(a, deinterlace) $
   else begin
      im = image(deinterlace(a, deinterlace), axis_style = 2)
      pf = plot([0], [0], over = im, $
                linestyle = '', symbol = 'o', color = 'red')
      graphics = {lmfgraphics, im:im, pf:pf}
   endelse
endif

quiet = ~keyword_set(debug)
   
;;; Find candidate features
rp = ctfeature(a, noise = noise, deinterlace = deinterlace, $
               pickn = pickn, count = nfeatures)

if nfeatures le 0 then $
   return, features

if dographics then $
   graphics.pf -> putdata, rp[0, *], rp[1, *]

if doreport then begin
   print
   message, 'noise: ' + strtrim(noise, 2), /inf
   if keyword_set(deinterlace) then $
      message, 'analyzing ' + ((deinterlace mod 2) ? 'odd' : 'even') + ' field', /inf
   message, 'features found: ' + strtrim(nfeatures, 2), /inf
endif

if ~isa(nfringes, /scalar, /number) then nfringes = 20
if ~isa(maxrefits, /scalar, /number) then maxrefits = 1

;;; Loop over features
count = 0L
for ndx = 0L, nfeatures - 1 do begin
   if doreport then message, 'feature ' + strtrim(ndx, 2), /inf

   ;; Use starting estimate for particle position
   rc = rp[0:1, ndx]

   dorefit = 0
refit:

   aa = aziavg(a, center = rc, deinterlace = deinterlace, deviates = dev)

   ;;; analyze radial profile
   rn = extrema(aa, ismin = ismin, count = nfound) ; coordinates of maxima and minima
   if nfound le 1 then continue                    ; flatline: nothing here
   w = where(rn gt 2, nfound)                      ; ignore noise near center
   if nfound le 1 then continue
   rn = rn[w]
   ismin = ismin[w]

   ;;; region of interest
   n = nfringes < (nfound - 1)     ; range set by fringe number.
   range = rn[n]                
   xc = rc[0]
   yc = rc[1]
   x0 = round(xc - range) > 0 < maxx
   x1 = round(xc + range) > 0 < maxx
   y0 = round(yc - range) > 0 < maxy
   y1 = round(yc + range) > 0 < maxy
   r0 = [x0, y0]                ; origin of ROI
   roi = [x0, x1, y0, y1]
   if dographics then begin
      poly = [[x0, y0], [x1, y0], [x1, y1], [x0, y1], [x0, y0]]
      proi = plot(poly, over = graphics.im, linestyle = '--', color = 'light green')
   endif

   ;;; cropped and deinterlaced iamge
   aa = aa[0:range]
   if keyword_set(deinterlace) then begin
      ac = a[x0:x1, y0:y1:2]
      n0 = ((y0 + deinterlace) mod 2) ? ceil(y0/2) : floor(y0/2)
      n1 = n0 + n_elements(ac[0, *]) - 1
      dev = dev[x0:x1, n0:n1]
   endif else begin
      ac = a[x0:x1, y0:y1]
      dev = dev[x0:x1, y0:y1]
   endelse

   ;;; significance estimate within ROI
   err = abs(dev)/noise > 1.

   ;; Use radial profile to estimate axial position, zp,
   ;; David Ruffner's method based on spherical wave model.
   if dozp then begin
      if doap then begin
         rho = findgen(range) + 0.5
         lap = deriv(aa)        ; Laplacian of azimuthal profile
         lap = deriv(lap) + lap/rho
         w = where(abs(aa-1.) gt 1e-2)
         qsq = -lap[w]/(aa[w] - 1.)
         zsq = rho[w]^2 * ((k*mpp)^2/qsq - 1.)
         w = where(zsq gt 0., ngood)
         if ngood lt 2 then begin
            message, 'could not estimate zp -- skipping this feature', /inf
            continue
         endif
         sigmatrim, zsq[w], mzsq
         zp = sqrt(mzsq)        ; estimated axial position [pixel]
      endif else $
         zp = ap*(4.*k)/median(j0n/rn[where(ismin)])
   endif

   ;; Model observed interference pattern as Poisson's spot to
   ;; obtain radius, ap, from axial position, zp
   if doap then $
      ap = zp/(4.*k)*median(j0n/rn[where(ismin)])

   ;; Estimate np: FIXME: currently uses input value: np0

   ;; Estimate alpha and delta
   alpha = 1.d
   delta = 0.d

   ;; Starting Estimates
   p0 = [zp, ap, real_part(np0), imaginary(np0), $
         real_part(nm0), imaginary(nm0), alpha, delta] ; initial estimates

   if ~quiet then lmf_report, 'starting estimates:', [rc, p0]

   ;; Improve estimates by fitting to radial profile
   p1 = fitlmsphere1d(aa, p0, lambda, mpp, $
                      fixap = fixap, fixnp = fixnp, fixdelta = fixdelta, $
                      chisq = thischisq, /quiet)
   if ~finite(thischisq) then begin
      message, 'parameter estimate failed -- trying again with fixed delta', /inf, $
               noprint = quiet
      p1 = fitlmsphere1d(aa, p0, lambda, mpp, $
                         fixap = fixap, fixnp = fixnp, /fixdelta, $
                         chisq = thischisq, /quiet)
   endif else begin
      ;; Some fits fail with alpha = 2.0; have to fixdelta.
      peggedalpha = (p1[0, 6] ge 1.9) && ~fixdelta
      if peggedalpha then begin
         message, 'alpha exceeds bounds -- trying again with fixed delta', /inf, $
                  noprint = quiet
         p1 = fitlmsphere1d(aa, p0, lambda, mpp, $
                            fixap = fixap, fixnp = fixnp, /fixdelta, $
                            chisq = thischisq, /quiet)
      endif
   endelse
   if ~finite(thischisq) then begin
      message, 'parameter estimate failed -- skipping this feature', /inf
      if dographics then proi.delete
      continue
   endif
   if max(abs(p1[0, 0:2] - p0[0:2])/p0[0:2]) gt 0.5 then begin
      message, '1D estimate performed badly -- reverting to initial estimates', /inf, $
               noprint = quiet
      p1[0, *] = p0
   endif

   if ~quiet then lmf_report, 'improved estimates: ' + string(thischisq), $
                              [rc, reform(p1[0, *])]

   ;; 2D fit to refine estimates
   p2 = [rc-r0, reform(p1[0, *])]  ; initial parameters from 1D fit

   thisfixdelta = fixdelta or peggedalpha
   thisdeinterlace = keyword_set(deinterlace) ? deinterlace + r0[1] : 0
   thisfeature = fitlmsphere(ac, p2, lambda, mpp, $
                             errors = err, $
                             chisq = thischisq, $
                             residuals = residuals, $
                             fixap = fixap, $
                             fixnp = fixnp, $
                             fixalpha = fixalpha, $
                             fixdelta = thisfixdelta, $
                             deinterlace = thisdeinterlace, $
                             resolution = resolution, $
                             gpu = gpu, $
                             quiet = quiet)

   ;;; Handle failure
   ;; 1. Fit did not converge
   if n_elements(thisfeature) eq 1 then begin ; fit failed
      message, 'fit failed -- skipping this feature', /inf
      if dographics then proi.delete
      continue
   endif
   ;; 2. Center shifted by more than half a pixel in plane
   ;;    or large discrepancy in estimated size.
   ;;    Accept new center and refit -- only do this once
   thisfeature[0, 0:1] += r0    ; center in original image
   if (abs(ap - thisfeature[0, 3])/ap gt 0.5) then begin
      if dorefit++ lt maxrefits then begin
         rc = reform(thisfeature[0, 0:1])
         message, 'large displacement -- refitting', /inf
         if dographics then begin
            rp[0:1, ndx] = rc
            graphics.pf -> putdata, rp[0, *], rp[1, *]
            proi.delete
         endif
         goto, refit
      endif
   endif

   if doreport then lmf_report, 'chisq: ' + string(thischisq), thisfeature[0, *]

   if dographics then begin
      proi.delete
      rp[0:1, ndx] = thisfeature[0,0:1]
      graphics.pf -> putdata, rp[0,*], rp[1,*]
   endif

   thisfeature = [[thisfeature], [thischisq, range]]
   features.add, thisfeature
   count++
endfor

return, features.toarray()
end
