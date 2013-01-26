;+
; NAME:
;      lmsphereprofile
;
; PURPOSE:
;      Uses Lorenz-Mie theory to calculate the radial profile 
;      of the in-line hologram of a sphere, as obtained with
;      holographic video microscopy.
;
; CATEGORY:
;      Digital holography, microscopy
;
; CALLING SEQUENCE:
;      p = lmsphereprofile(rho,zp,ap,np,nm,alpha,delta,lambda,mpp)
;
; INPUTS:
;      rho: [npts] array of radii [pixels]
;      zp: Axial position of the sphere relative to the microscope's
;          focal plane [pixels]
;      ap: Particle radius [micrometers]
;      np: Complex refractive index of particle
;      nm: Complex refractive index of medium
;      alpha: Relative amplitude of illumination
;      delta: Wavefront distortion [pixels]
;      lambda: vacuum wavelength of light [micrometers]
;      mpp: length-scale calibration [micrometers/pixel]
;
; OUTPUTS:
;      p: computed hologram profile.
;
; PROCEDURE:
;      Calls SPHEREFIELD to compute Lorenz-Mie scattering pattern
;      for sphere.
;
; REFERENCE:
;      S. Lee, Y. Roichman, G. Yi, S. Kim, S. Yang, A. van Blaaderen,
;      P. van Oostrum and D. G. Grier,
;      Chararacterizing and tracking
;      single colloidal particles with video holographic microscopy,
;      Optics Express 15, 18275-18282 (2007)
;
; MODIFICATION HISTORY:
; 11/04/2007 Written by David G. Grier, New York University.
;           Based on earlier version called HSA.
; 02/08/2008 DGG revised to work with updated SPHEREFIELD syntax.
; 04/15/2008 DGG Added MPP keyword.  Associated documentation fixes.
; 10/12/2012 DGG Revised interference code.  Added support for DELTA.
;    Renamed to lmsphereprofile.
; 01/26/2013 DGG Correct documentation
;
; Copyright (c) 2007-2013 David G. Grier
;-

function lmsphereprofile, rho, zp, ap, np, nm, alpha, delta, lambda, mpp

COMPILE_OPT IDL2

npts = n_elements(rho)

x = reform(rho, 1, npts)
y = dblarr(1, npts)

dhm = spherefield(x, y, zp, ap, np, nm, lambda, mpp, k = k, /cartesian)
dhm *= alpha * exp(dcomplex(0,-k*(zp + delta)))
dhm[0, *] += 1.

return, real_part(total(dhm * conj(dhm), 1))
end
