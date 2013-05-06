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
;
; Copyright (c) 2008-2013 David G. Grier, David Ruffner and Fook Chiong Cheong
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
; LMF_CROP
;
; Use the azimuthal average and the azimuthal standard deviation to 
; estimate the signal to noise ratio as a function of distance from a
; specified center.  Use this to establish a region of interest to
; crop from a hologram and also the statistical weighting to apply to
; each pixel in that region.
;
pro lmf_crop, a, rc, aa, roi, ac, weight, deinterlace = deinterlace

COMPILE_OPT IDL2, HIDDEN

sz = size(a, /dimensions)

as = azistd(a, aa, center = rc, deinterlace = deinterlace) ; azimuthal standard deviation
aa -= 1.                                                   ; azimuthal average

;;; region of interest
range = max(where(abs(aa) ge as)) > 30 ; estimate for range of useful signal
x0 = round(rc[0] - range) > 0 < (sz[0]-1)
x1 = (x0 + 2*range + 1) > 0 < (sz[0]-1)
y0 = round(rc[1] - range) > 0 < (sz[1]-1)
y1 = (y0 + 2*range + 1) > 0 < (sz[1]-1)
roi = [[x0, y0], [x1, y1]]      ; corners of ROI

;;; cropped image
ac = a[x0:x1, y0:y1]

;;; significance estimate within ROI
; dr = savgol(5, 5, 1, 3) ; Savitzky-Golay derivative filter
dr = [0.0582750, -0.0571095, -0.103341, -0.0977078, -0.0574980, 0.00000, $
      0.0574980,  0.0977078,  0.103341,  0.0571095,  -0.0582750]
sgn = fix(convol(aa, dr, /edge_truncate) ge 0) ; sign of derivative: 1 if non-negative
r = where(abs(sgn-sgn[1:*]) gt 0.5, npeaks)    ; zero crossing when sign changes
snr = abs(aa[r])/as[r]                         ; signal-to-noise ratio
;ratio at extrema
nx = x1 - x0 + 1
ny = y1 - y0 + 1
xsq = rebin((findgen(nx)    + x0 - rc[0])^2, nx, ny, /sample)
ysq = rebin((findgen(1, ny) + y0 - rc[1])^2, nx, ny, /sample)
err = reform(interpol(snr, r, sqrt(xsq + ysq), /spline), nx, ny) > min(snr)
end

;;;;;
;
; LMFEATURE
;
; Main routine
;
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
quiet = ~keyword_set(debug)

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

;;; Loop over features
for ndx = 0L, count - 1 do begin
   if doreport then message, 'feature ' + strtrim(ndx, 2), /inf

   ;; Initial estimate for xp and yp
   rp0 = rp[0:1, ndx]           ; particle position in original image

   dorefit = 0
refit:
   ;; Crop image to region of interest
   lmf_crop, a, rp0, aa, roi, ac, weight, deinterlace = deinterlace
   rc = rp0 - roi[0:1]          ; particle position in cropped image

   if dographics then begin
      poly = roi[[[0,1],[2,1],[2,3],[0,3],[0,1]]]
      proi = plot(poly, over = graphics, linestyle = '--', color = 'light green')
   endif

   ;; Use radial profile to estimate axial position, zp,
   ;; David Ruffner's method based on spherical wave model.
   rho = findgen(n_elements(aa)) + 0.5
   lap = deriv(aa)
   lap = deriv(lap) + lap/rho
   w = where(abs(aa) gt 1e-2)
   qsq = -lap[w]/aa[w]
   zsq = rho[w]^2 * ((k*mpp)^2/qsq - 1.)
   sigmatrim, zsq, mzsq
   zp = sqrt(mzsq)              ; estimated axial position [pixel]

   ;; Model observed interference pattern as Poisson's spot to
   ;; obtain radius, ap, from axial position, zp
   if doap then begin
      c = fix(aa ge 0.)
      z0 = float(where(abs(c - c[1:*]) gt 0)) + 1. ; FIXME: zero crossings
      j0n = [2.4048, 5.5201, 8.6537]               ; zeros of J0(x)
      fac = 2.*k*mean(z0/j0n) - 1.
      ap = zp /fac              ; estimated radius [um]
   endif else $
      ap = ap0

   ;; Estimate np: FIXME: currently uses input value: np0

   ;; Estimate alpha and delta
   alpha = 1.d
   delta = 0.d

   p0 = [zp, ap, real_part(np0), imaginary(np0), $
         real_part(nm0), imaginary(nm0), alpha, delta] ; initial estimates

   if ~quiet then lmf_report, 'starting estimates:', [rp0, p0]

   ;; Improve estimates by fitting to radial profile
   p1 = fitlmsphere1d(aa+1., p0, lambda, mpp, fixdelta = fixdelta, chisq = thischisq, /quiet)
   if ~finite(thischisq) then begin
      message, 'parameter estimate failed -- trying again with fixed delta', /inf, noprint = quiet
      p1 = fitlmsphere1d(aa+1, p0, lambda, mpp, /fixdelta, chisq = thischisq, /quiet)
   endif else begin
      ;; Some fits fail with alpha = 2.0; have to fixdelta.
      peggedalpha = (p1[0, 6] ge 1.9) && ~fixdelta
      if peggedalpha then begin
         message, 'alpha exceeds bounds -- trying again with fixed delta', /inf, noprint = quiet
         p1 = fitlmsphere1d(aa+1, p0, lambda, mpp, /fixdelta, chisq = thischisq, /quiet)
      endif
   endelse
   if ~finite(thischisq) then begin
      message, 'parameter estimate failed -- skipping this feature', /inf
      continue
   endif
   
   if ~quiet then lmf_report, 'improved estimates: ' + string(thischisq), $
                              [rp0, reform(p1[0, *])]

   if dographics then begin
      rp[0:1,ndx] = rp0
      pf.putdata, rp[0,*], rp[1,*]
   endif

   ;; 2D fit to refine estimates
   p2 = [rc, reform(p1[0, *])]  ; initial parameters from 1D fit
   thisfeature = fitlmsphere(ac, p2, lambda, mpp, $
                             weight = weight, $
                             chisq = thischisq, $
                             fixalpha = fixalpha, $
                             fixdelta = fixdelta || peggedalpha, $
                             deinterlace = keyword_set(deinterlace) ? deinterlace + roi[1] : 0, $
                             object = gpu, $
                             quiet = quiet)
   if n_elements(thisfeature) eq 1 then begin ; fit failed
      message, 'fit failed -- skipping this feature', /inf
      continue
   endif
   if max(abs(rc - thisfeature[0, 0:1])) gt 0.6 then begin
      if dorefit eq 0 then begin
         dorefit = 1
         rp0 = reform(thisfeature[0, 0:1]) + roi[0:1]
         message, 'large displacement -- refitting', /inf
         goto, refit
      endif
   endif

   thisfeature[0,0:1] += roi[0:1] ; center in original image
   
   if doreport then lmf_report, 'chisq: ' + string(thischisq), thisfeature[0, *]

   if dographics then begin
      proi.setproperty, color = 'green'
      rp[0:1,ndx] = thisfeature[0,0:1]
      pf.putdata, rp[0,*], rp[1,*]
   endif

   thisfeature = [[thisfeature], [thischisq, n_elements(aa)]]
   features.add, thisfeature
endfor

return, features
end
