;+
; NAME:
;    feature2dhm
;
; PURPOSE:
;    Calculate a digital holographic microscopy image of one or more 
;    spheres given parameters obtained from Lorenz-Mie fitting.
;
; CATEGORY:
;    Holographic video microscopy
;
; CALLING SEQUENCE:
;    hologram = feature2dhm(features, dim)
;
; INPUTS:
;    features: Array of Lorenz-Mie fitting parameters in the format
;        returned by DHMFEATURE
;
;    dim: [w,h] Dimensions of hologram to be calculated [pixel]
;
; KEYWORD PARAMETERS:
;    lambda: Wavelength of light [micrometers]
;
;    mpp: Micrometers per pixel
;
; KEYWORD FLAGS:
;    gpu: If set, use GPU to calculate hologram
;
; OUTPUTS:
;    hologram: Floating point array with dimensions given by dim.
;
; MODIFICATION HISTORY:
; 07/17/2012 Written by David G. Grier, New York University
;
; Copyright (c) 2012 David G. Grier
;-

function feature2dhm, features, dim, lambda = lambda, mpp = mpp, gpu = gpu

npts = n_elements(features[0,0,*])

if keyword_set(gpu) then $
   dhm = dggdhmspheredhm(dim=dim,lambda=lambda,mpp=mpp)

a = fltarr(dim)

for n = 0, npts-1 do begin
   r0 = (dim - 1.)/2.
   rp = reform(features[0,0:2,n])
   rp[0:1] -= r0
   ap = features[0,3,n]
   np = features[0,4,n]
   alpha = features[0,5,n]
   nm = features[0,6,n]
   if keyword_set(gpu) then begin
      dhm.setproperty, rp = rp, ap = ap, np = np, alpha = alpha, nm = nm
      a += dhm.hologram
   endif else $
      a += spheredhm(rp, ap, np, nm, dim, alpha = alpha, $
                     lambda = lambda, mpp = mpp)
endfor
return, a - npts + 1.
end
