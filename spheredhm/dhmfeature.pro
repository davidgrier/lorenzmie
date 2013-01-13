;+
; NAME:
;    dhmfeature
;
; PURPOSE:
;    Identify, locate and characterize spheres in
;    normalized holographic video microscopy images.
;
; CATEGORY:
;    Image analysis, holographic video microscopy, feature detection
;
; CALLING SEQUENCE:
;    p = dhmfeature(image)
;
; INPUTS:
;    image: two-dimensional normalized holographic microscopy image.
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
;    rad: Radius over which to fit a given feature [pixels]
;        Default: estimated from image
;
;    rp: zp, or (xp, yp), or (xp, yp, zp)
;        Initial estimates for the centers of the particle [pixels]
;        Default: estimated from image
;
;    ap: ballpark radius of sphere [micrometers]
;        Default: estimated from image
;
;    np: ballpark refractive index of particle.
;        Default: 1.5
;
;    nm: refractive index of medium.
;        Default: water at room temperature for the given wavelength.
;
;    lambda: vacuum wavelength of light [micrometers]
;        Default: 0.532
;
;    mpp: length calibration [micrometers per pixel]
;        Default: 0.135
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
;    quiet: If set, report only errors
;
; OUTPUTS:
;    p: [7,2,npts] array of data for the features found in the image.
;       npts: Number of features found
;       p[0,*,n]: [xp, dxp] : x location and error [pixels]
;       p[1,*,n]: [yp, dyp] : y location and error [pixels]
;       p[2,*,n]: [zp, dzp] : z location and error [pixels]
;       p[3,*,n]: [ap, dap] : radius and error [micrometers]
;       p[4,*,n]: [np, dnp] : refractive index and error
;       p[5,*,n]: [alpha, dalpha] : relative illumination at rp, and error
;       p[6,0,n]: nm: refractive index of medium
;      
; PROCEDURE:
; CIRCLETRANSFORM highlights positions of sphere in the plane.
; FASTFEATURE identifies initial particle positions in the plane.
; For each feature:
;   CT_HITORMISS estimates range around each feature.
;   RS1D estimates zp.
;   Poisson-spot model estimates ap and improves zp.
;   FITSPHEREDHM1D refines estimate for zp, ap and np.
;   FITSPHEREDHM yields final estimates for parameters
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
; 07/18/2012 DGG Improved algorithm for estimating RAD using CT_HITORMISS.
;
; Copyright (c) 2008-2012 David G. Grier and Fook Chiong Cheong
;-

function dhmfeature, a, $
                     noise = noise, $
                     threshold = threshold, $
                     pickn = pickn, $
                     rad = rad, $
                     rp = rp, $
                     ap = ap, $
                     np = np, $
                     nm = nm, $
                     lambda = lambda, $
                     mpp = mpp, $
                     gpu = gpu, $
                     deinterlace = deinterlace, $
                     graphics = graphics, $
                     quiet = quiet, $
                     debug = debug

COMPILE_OPT IDL2

;;; Process command line parameters

umsg = 'USAGE: p = dhmfeature(a)'

if n_params() ne 1 then begin
   message, umsg, /inf
   return, -1
endif

sz = size(a)
if ~isa(a, /number, /array) or sz[0] ne 2 then begin
   message, umsg, /inf
   message, 'A must be a two-dimensional normalized hologram', /inf
   return, -1
endif

if ~isa(noise, /number, /scalar) then $
   noise = 0.05d*mean(a)

dorad = ~isa(rad, /number, /scalar) ; estimate radius for each feature
doap = ~isa(ap, /number, /scalar)
if ~isa(np, /number, /scalar) then $
   np = 1.5d

if ~isa(lambda, /number, /scalar) then $
   lambda = 0.532d

if ~isa(mpp, /number, /scalar) then $
   mpp = 0.135d

if ~isa(nm, /number, /scalar) then $
   nm = refractiveindex(lambda, 24.)

k = 2.d*!dpi*nm/lambda            ; wavenumber [radians/um]

debug = keyword_set(debug)
gpu = keyword_set(gpu)
doreport = ~keyword_set(quiet)
dographics = arg_present(graphics)

;;; Find candidate features
;; Transform ring-like patterns into spots
b = circletransform(a, noise = noise, deinterlace = deinterlace)

;; Centers of spots are estimates for particle centers: (xp, yp)
if ~isa(threshold, /number, /scalar) then begin
   threshold = 0.5d*max(b)
   print, 'setting threshold', threshold
endif

feature = fastfeature(b, threshold, pickn = pickn) ; find peaks

if n_elements(feature) lt 3 then begin
   message, string(threshold, $
                   format = '("no features found above threshold: ",F0)'), /inf
   return, -1
endif

npts = n_elements(feature[0,*])

if isa(pickn,/scalar,/number) && npts lt pickn then begin
   message, 'PICKN: not enough features found', /inf
endif

if doreport then $
   message, string(npts, format = '(I0," features found")'), /inf

if dographics then begin
   if isa(graphics, 'graphicswin') then $
      graphics.erase $
   else $
      graphics = window()
   im = image(a, axis_style = 2, /current)
   pf = plot(feature[0,*], feature[1,*], /over, $
             linestyle = '', symbol = 'o', color = 'red')
endif

;;; Loop over features to process each feature
p = []
precision = 2                   ; NOTE: make into a parameter?
for ndx = 0L, npts - 1 do begin
   ;; Initial estimate for xp and yp
   rp = reform(feature[0:1, ndx])

   ;; Crop image to isolate feature
   if dorad then begin          ; estimate radius
      hm = ct_hitormiss(a, rp, noise = noise, precision = precision, $
                        deinterlace = deinterlace, /coordinates)
      rad = round(100.*sqrt(total(hm[2,*])/n_elements(hm[2,*])))
   endif
   x0 = round(rp[0] - rad) > 0L
   x1 = round(rp[0] + rad) < sz[1] - 1L
   y0 = round(rp[1] - rad) > 0L
   y1 = round(rp[1] + rad) < sz[2] - 1L
   if dographics then begin
      poly = [[x0, y0], [x1, y0], [x1, y1], [x0, y1], [x0, y0]]
      roi = plot(poly, /over, linestyle = '--', color = 'light green')
   endif
   aa = a[x0:x1, y0:y1]         ; cropped image
   r0 = double([x0, y0])
   r1 = double([x1, y1])
   rc = rp - r0                 ; particle position in cropped image
   origin = (r0 + r1)/2.d       ; center of cropped image
   drc = rp - origin            ; deviation from middle of cropped image

   ;; Estimate axial position as position of peak brightness in the
   ;; Rayleigh-Sommerfeld reconstruction: zp
   z = dindgen(50) * 4.d + 10.d ; NOTE: set range more intelligently
   res = rs1d(aa, z, rc, lambda = lambda/nm, mpp = mpp)
   m = max(abs(res), loc)
   zmax = z[loc]                ; axial position of intensity maximum [pixels]

   b = aziavg(aa, center = rc, rad = rad, $
              deinterlace = deinterlace) ; radial profile around center
   
   ;; Model observed interference pattern as Poisson's spot
   ;; to estimate radius of particle: ap
   if doap then begin
      ;; Find radii of "zero" crossings
      c = fix(b ge 1)
      z0 = float(where(abs(c - c[1:*]) gt 0)) + 1. ; zero crossings [pixels]
      ;; Compare radii to those of Bessel function
      x0 = [2.4048, 5.5201]                        ; first zeros of J0(x)
      ap = zmax /(2.*k*mean(z0/x0) - 1.)           ; estimated radius [um]
   endif

   zp = zmax + ap/mpp           ; improved axial position estimate [pixels]

   ;; Estimate np NOTE: FIXME!

   ;; Estimate alpha
   alpha = 1.d

   ;; Improve estimates by fitting to azimuthally averaged image
   p0 = [zp, ap, np, alpha, nm]
   p1d = fitspheredhm1d(b, p0, lambda = lambda, mpp = mpp, /fixnm, quiet = ~debug)
   zp = p1d[0,0]
   ap = p1d[0,1]
   np = p1d[0,2]
   alpha = p1d[0,3]

   if doreport then begin
      message, 'feature: '+strtrim(ndx,2), /inf
      message, string(rp, zp, rad, $
                      format = '("  rp = (",I0,",",I0,",",I0,") +/- ",I0)'), /inf
      message, string(ap, np, $
                      format = '("  ap = ",F0," um, np = ",F0)'), /inf
   endif

   if dographics then begin
      feature[0:1,ndx] = rp
      pf.putdata, feature[0,*], feature[1,*]
   endif

   ;; Initial parameters
   p0 = [drc, zp, ap, np, alpha, nm]

   ;; Fit to refine estimates
   thisp = fitspheredhm(aa, p0, lambda = lambda, mpp = mpp, /fixnm, $
                        deinterlace = deinterlace, $
                        object = gpu, quiet = ~debug)
   thisp[0,0:1] += origin       ; center in original image
   
   if dographics then begin
      roi.setproperty, color = 'green'
      feature[0:1,ndx] = thisp[0,0:1]
      pf.putdata, feature[0,*], feature[1,*]
   endif

   p = [[[p]], [[thisp]]]
endfor

return, p
end
