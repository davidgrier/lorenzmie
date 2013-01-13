;+
; NAME:
;    dhmtool
;
; PURPOSE:
;    Given a holographic video microscope snapshot containing
;    the image of at least one sphere, perform a Lorenz-Mie fit 
;    the region around the sphere.  Return the refined fits for 
;    the sphere's radius, refractive index and position relative 
;    to the image's origin.
;
; CATEGORY:
;    Holographic video microscopy
;
; CALLING SEQUENCE:
;    p = dhmtool(a, zp, ap, np, nm, alpha, delta, lambda, mpp)
;
; INPUTS:
;    a: holographic snapshot containing images of spheres
;    ap: sphere's radius [micrometers]
;    np: sphere's refractive index.
;    nm: refractive index of medium.
;    alpha: estimate for relative illumination amplitude at sphere.
;        1 is a good starting point.
;    delta: estimate for wavefront distortion at sphere.
;        0 is a good starting point.
;    lambda: vacuum wavelength of light [micrometers]
;    mpp: micrometers per pixel.
;
; KEYWORD PARAMETERS:
;    rp: [xp, yp, zp] sphere's centroid [pixels].  
;        Default: Found algorithmically
;
;    rad: range of interest [pixels].
;        Default: Determined algorithmically
;
;    deinterlace: Fit only even field (DEINTERLACE = 0) or 
;        odd field (DEINTERLACE = 1) of (interlaced) hologram.
;        Default: No deinterlacing
;
;    aplimits: [minap, maxap] limits on ap [micrometers]
;           Default: [0.05, 10.]
;
;    nplimits: [minnp, maxnp] limits on np
;           Default: [1.01 nm, 3.]
;
; KEYWORD FLAGS:
;    fixnp: do not adjust np; p[1,4] = 0
;    fixalpha: do not adjust alpha; p[1,5] = 0
;    fixnm: do not adjust nm; p[1,6] = 0
;
;    gpu: use GPU acceleration, if available.
;
;    quiet: do not provide diagnostic output.
;
; OUTPUTS:
;    p: [2,7] array of fitting parameters and error estimates
;    p[0,0]: xp [pixels]
;    p[0,1]: yp [pixels]
;    p[0,2]: zp [pixels]
;    p[0,3]: ap [micrometers]
;    p[0,4]: np : refractive index of sphere
;    p[0,5]: kp : extinction coefficient of sphere
;    p[0,6]: nm : refractive index of medium
;    p[0,7]: km : extinction coefficient of medium
;    p[0,8]: alpha
;    p[0,9]: delta
;
;    p[1,*] contains associated error estimates.
;
; SIDE EFFECTS:
;    If RC is not provided, then DHMTOOL presents an image
;    and requires the user to click on the sphere's center.
;    Unless QUIET is set, DHMTOOL provides diagnostic output.
;
; RESTRICTIONS:
;    Requires a normalized hologram.
;
; PROCEDURE:
;    Crop region around the estimated centroid.  Perform fit
;    on cropped image using FITSPHEREDHM.  Adjust fit parameters 
;    to account for cropping.
;
; MODIFICATION HISTORY:
; 01/15/2009 Written by David G. Grier, New York University.
;   This formalizes a routine that has been in use since 2007.
; 02/14/2009 DGG. Added APLIMITS and NPLIMITS keywords.
;   Documentation fixes.
; 05/04/2011 Ke Xiao Corrected a bug in deinterlace code and
;   in estimating particles' centers.
; 09/
; Copyright (c) 2007-2011 David G. Grier and Ke Xiao
;-

function dhmtool, a, $               ; hologram
                  zp, $              ; z offset of particle [pixel]
                  ap, $              ; radius of particle [micrometers]
                  np, $              ; refractive index of particle
                  alpha, $           ; fraction of scattered light
                  nm, $              ; refractive index of medium
                  rc = rc, $         ; centroid [pixels]
                  rad = rad, $       ; range of interest [pixels]
                  lambda = lambda, $ ; wavelength of light [micrometers]
                  mpp = mpp, $       ; micrometers per pixel
                  aplimits = aplimits, $
                  nplimits = nplimits, $
                  fixnp = fixnp, $
                  fixnm = fixnm, $
                  fixalpha = fixalpha, $
                  deinterlace = deinter, $
                  lut = lut, gpu = gpu, $ ; acceleration
                  quiet = quiet           ; minimize output

if n_elements(rad) lt 1 then rad = 50

; estimate center of scattering pattern, if not provided
if n_elements(rc) ne 2 then begin
    plotimage, bytscl(a), /iso
    print, "Click on center of sphere"
    cursor, x, y, /data
    xc = x
    yc = y
    rc = [xc, yc]
endif else begin
    xc = rc[0]
    yc = rc[1]
endelse

chatty = ~keyword_set(quiet)

if n_elements(deinter) ge 1 then $
   a = deinterlace(a, odd = deinter)

; clip hologram to region around center

sz = size(a, /dimensions)
x0 = round(xc - rad) > 0
x1 = round(xc + rad) < sz[0]-1
y0 = round(yc - rad) > 0
y1 = round(yc + rad) < sz[1]-1

ac = a[x0:x1, y0:y1]            ; cropped image

if chatty then $
   plotimage, bytscl(ac), /iso

; dimensions of cropped image
nx = x1 - x0 + 1
ny = y1 - y0 + 1

; particle position relative to center of cropped image
xp = xc - x0 - nx/2.
yp = yc - y0 - ny/2.

; show initial estimates
if chatty then begin $
   tvscl, ac
   b = spheredhm([xp, yp, zp], ap, np, nm, [nx, ny], $
                 alpha = alpha, lambda = lambda, mpp = mpp)
   tvscl, b, nx+1, 0
endif
; perform fit
pin = [xp, yp, zp, ap, np, alpha, nm]

if n_elements(deinter) ge 1 then $
   p = fitspheredhm(ac, pin, $
                    lambda = lambda, mpp = mpp, $
                    aplimits = aplimits, nplimits = nplimits, $
                    fixnm = fixnm, fixnp = fixnp, $
                    fixap = fixap, fixalpha = fixalpha, $
                    lut = lut, gpu = gpu, $
                    deinterlace =  deinter, quiet = quiet) $
else $
   p = fitspheredhm(ac, pin, $
                    lambda = lambda, mpp = mpp, $
                    aplimits = aplimits, nplimits = nplimits, $
                    fixnm = fixnm, fixnp = fixnp, $
                    fixap = fixap, fixalpha = fixalpha, $
                    lut = lut, gpu = gpu, quiet = quiet)

if n_elements(p) eq 1 then return, -1

if chatty then begin
   b = spheredhm(p[0, 0:2], p[0, 3], p[0, 4], nm, [nx, ny], $
                   alpha = p[0, 5], lambda = lambda, mpp = mpp)
   tvscl, b, nx+1, 0
   print, "Height: "+strtrim(p[0, 2], 2)+" +/- "+strtrim(p[1, 2], 2)+" [pixels]"
   print, "Radius: "+strtrim(p[0, 3], 2)+" +/- "+strtrim(p[1, 3], 2)+" [micron]"
   print, "Index : "+strtrim(p[0, 4], 2)+" +/- "+strtrim(p[1, 4], 2)
endif

p[0,0:1] += [x0+nx/2., y0+ny/2.]

return, p
end
