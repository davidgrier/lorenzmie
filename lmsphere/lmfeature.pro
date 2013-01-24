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
;    p = lmfeature(image, lambda, mpp)
;
; INPUTS:
;    image: two-dimensional normalized holographic microscopy image.
;    lambda: vacuum wavelength of light [micrometers]
;    mpp: length calibration [micrometers per pixel]
;
; KEYWORD PARAMETERS:
;    Parameters for 2D feature location
;    noise: Estimate for the RMS additive noise at each pixel.
;        Default: estimated from image
;
;    threshold: Minimum weighting for a pixel to be considered part
;        of a feature.
;        Default: estimated from image
;
;    pickn: Maximum number of features to consider.
;        Default: Process all features above threshold
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
;    p: [10,2,npts] array of data for the features found in the image.
;       npts: Number of features found
;       p[0,*,n]: [xp, dxp] : x location and error [pixels]
;       p[1,*,n]: [yp, dyp] : y location and error [pixels]
;       p[2,*,n]: [zp, dzp] : z location and error [pixels]
;       p[3,*,n]: [ap, dap] : radius and error [micrometers]
;       p[4,*,n]: [np, dnp] : refractive index and error of sphere
;       p[5,*,n]: [kp, dnp] : extinction coefficient and error of sphere
;       p[6,*,n]: [nm, dnm] : refractive index of medium
;       p[7,*,n]: [km, dkm] : extinction coefficient of medium
;       p[8,*,n]: [alpha, dalpha] : relative illumination amplitude
;       p[9,*,n]: [delta, ddelta] : wavefront distortion
;      
; PROCEDURE:
;   CT_FEATURE identifies candidate features using CIRCLETRANSFORM.
;   CT_RANGE estimates range around each feature.
;   RS1D estimates zp.
;   Poisson-spot model estimates ap and improves zp.
;   FITSPHERELM1D refines estimate for zp, ap and np.
;   FITSPHERELM yields final estimates for parameters
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
;
; Copyright (c) 2008-2012 David G. Grier and Fook Chiong Cheong
;-

function lmfeature, a, lambda, mpp, $
                    noise = noise, $
                    threshold = threshold, $
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

;;; Process command line parameters
if n_params() ne 3 then begin
   message, umsg, /inf
   return, -1
endif

sz = size(a)
if ~isa(a, /number, /array) or sz[0] ne 2 then begin
   message, umsg, /inf
   message, 'A must be a two-dimensional normalized hologram', /inf
   return, -1
endif
width = sz[1]
height = sz[2]

if ~isa(lambda, /number, /scalar) then begin
   message, umsg, /inf
   message, 'LAMBDA must be a number', /inf
   return, -1
endif

if ~isa(mpp, /number, /scalar) then begin
   message, umsg, /inf
   message, 'MPP must be a number', /inf
   return, -1
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
rp = ctfeature(a, noise = noise, threshold = threshold, $
               pickn = pickn, count = count, deinterlace = deinterlace)

if doreport then $
   message, string(count, (count ne 1) ? 's' : '', $
                   format = '(I0," feature",A," found")'), /inf

if dographics then begin
   if isa(graphics, 'graphicswin') then $
      graphics.erase $
   else $
      graphics = window()
   im = image(a, axis_style = 2, /current)
   pf = plot(rp[0,*], rp[1,*], /over, $
             linestyle = '', symbol = 'o', color = 'red')
endif

if count ge 1 then $
   rad = ct_range(a, rp, noise = noise, deinterlace = deinterlace)

;;; Loop over features to process each feature
p = []
for ndx = 0L, count - 1 do begin
   ;; Initial estimate for xp and yp
   rp0 = rp[0:1, ndx]

   ;; Crop image to region of interest
   thisrad = rad[ndx]
   x0 = round(rp0[0] - thisrad) > 0L
   x1 = round(rp0[0] + thisrad) < width - 1L
   y0 = round(rp0[1] - thisrad) > 0L
   y1 = round(rp0[1] + thisrad) < height - 1L
   if dographics then begin
      poly = [[x0, y0], [x1, y0], [x1, y1], [x0, y1], [x0, y0]]
      roi = plot(poly, /over, linestyle = '--', color = 'light green')
   endif
   aa = a[x0:x1, y0:y1]         ; cropped image
   r0 = double([x0, y0])
   r1 = double([x1, y1])
   rc = rp0 - r0                ; particle position in cropped image
   origin = (r0 + r1)/2.d       ; center of cropped image
   drc = rp0 - origin           ; particle displacement from center of cropped image

   ;; Estimate axial position as position of peak brightness in the
   ;; Rayleigh-Sommerfeld reconstruction: zp
   z = dindgen(50) * 4.d + 10.d ; NOTE: set range more intelligently
   res = rs1d(aa, z, rc, lambda = lambda/nm0, mpp = mpp)
   m = max(abs(res), loc)
   zmax = z[loc]                ; axial position of intensity maximum [pixels]

   b = aziavg(aa, center = rc, rad = thisrad, $
              deinterlace = deinterlace) ; radial profile around center
   
   ;; Model observed interference pattern as Poisson's spot
   ;; to estimate radius of particle: ap
   if doap then begin
      ;; Find radii of "zero" crossings
      c = fix(b ge 1)
      z0 = float(where(abs(c - c[1:*]) gt 0)) + 1. ; zero crossings [pixels]
      ;; Compare radii to those of Bessel function
      x0 = [2.4048, 5.5201, 8.6537]                ; zeros of J0(x)
      ap = zmax /(2.*k*mean(z0/x0) - 1.)           ; estimated radius [um]
   endif else $
      ap = ap0

   zp = zmax + ap/mpp           ; improved axial position estimate [pixels]

   ;; Estimate np: NOTE: currently uses input value: np0

   ;; Estimate alpha and delta
   alpha = 1.d
   delta = 0.d

   ;; Starting estimates
   if doreport then begin
      message, 'feature: '+strtrim(ndx,2), /inf
      message, 'starting estimates:', /inf
      message, string(rp0, zp, thisrad, $
                      format = '("  rp = (",I0,", ",I0,", ",I0,") +/- ",I0)'), /inf
      message, string(ap, real_part(np0), $
                      format = '("  ap = ",F0.3," um, np = ",F0.3)'), /inf
   endif      

   ;; Improve estimates by fitting to azimuthally averaged image
   p0 = [zp, ap, real_part(np0), imaginary(np0), $
         real_part(nm0), imaginary(nm0), alpha, delta] ; initial estimates
   p1 = fitlmsphere1d(b, p0, lambda, mpp, fixdelta = fixdelta, chisq = chisq, /quiet)

   ;; Some fits fail with alpha = 2.0; have to fixdelta.
   peggedalpha = (p1[0, 6] ge 1.9) && ~fixdelta
   if peggedalpha then begin
      message, 'alpha pegged, fixing delta', /inf
      p1 = fitlmsphere1d(b, p0, lambda, mpp, /fixdelta, /quiet)
   endif

   ;; Improved estimates
   if doreport then begin
      message, 'improved estimates:'+string(chisq), /inf
      message, string(rp0, p1[0, 0], thisrad, $
                      format = '("  rp = (",I0,", ",I0,", ",I0,") +/- ",I0)'), /inf
      message, string(p1[0, 1:2], $
                      format = '("  ap = ",F0.3," um, np = ",F0.3)'), /inf
      message, string(p1[0, 6:7], $
                      format = '("  alpha = ",F0.3,", delta = ",F0.3)'), /inf
   endif

   if dographics then begin
      rp[0:1,ndx] = rp0
      pf.putdata, rp[0,*], rp[1,*]
   endif

   ;; 2D fit to refine estimates
   p2 = [drc, reform(p1[0, *])] ; initial parameters from 1D fit
   thisp = fitlmsphere(aa, p2, lambda, mpp, $
                       chisq = chisq, $
                       fixalpha = fixalpha, $
                       fixdelta = fixdelta || peggedalpha, $
                       deinterlace = deinterlace, $
                       object = gpu, $
                       quiet = ~debug)
   if n_elements(thisp) eq 1 then begin ; fit failed
      message, 'fit failed -- continuing', /inf
      continue
   endif

   thisp[0,0:1] += origin       ; center in original image
   
   if doreport then begin
      message, 'refined estimates:'+string(chisq), /inf
      message, string(thisp[0, 0:2], $
                      format = '("  rp = (",F0.2,", ",F0.2,", ",F0.2,")")'), /inf
      message, string(thisp[0, 3:4], $
                      format = '("  ap = ",F0.3," um, np = ",F0.3)'), /inf
      message, string(thisp[0, 8:9], $
                      format = '("  alpha = ",F0.3,", delta = ",F0.3)'), /inf
      message, /inf
   endif

   if dographics then begin
      roi.setproperty, color = 'green'
      rp[0:1,ndx] = thisp[0,0:1]
      pf.putdata, rp[0,*], rp[1,*]
   endif

   p = [[[p]], [[thisp]]]
endfor

return, p
end
