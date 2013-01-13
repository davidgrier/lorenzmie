;+
; NAME:
;    spheredhmprofile
;
; PURPOSE:
;    Calculates the radial profile of the in-line hologram of a
;    sphere, as obtained with digital holographic microscopy.
;
; CATEGORY:
;    Digital holography, microscopy
;
; CALLING SEQUENCE:
;    p = spheredhmprofile(rho,zp,ap,np,nm,alpha)
;
; INPUTS:
;    rho: [npts] array of radii [pixels]
;    zp: Axial position of the sphere relative to the microscope's
;      focal plane [pixels]
;    ap: Particle radius [micrometers]
;    np: Complex refractive index of particle
;    nm: Complex refractive index of medium
;    alpha: Relative brightness of illumination (typically around 1).
;    delta: Wavefront distortion [micrometers]
;    lambda: vacuum wavelength of light [micrometers]
;    mpp: length-scale calibration [micrometers/pixel]
;
; OUTPUTS:
;    p: computed hologram profile.
;
; PROCEDURE:
;    Calls SPHEREFIELD to compute Lorenz-Mie scattering pattern
;      for sphere.
;
; REFERENCE:
;    S. Lee, Y. Roichman, G. Yi, S. Kim, S. Yang, A. van Blaaderen,
;    P. van Oostrum and D. G. Grier,
;    Chararacterizing and tracking
;    single colloidal particles with video holographic microscopy,
;    Optics Express 15, 18275-18282 (2007)
;
; MODIFICATION HISTORY:
; 11/04/2007 Written by David G. Grier, New York University.
;   Based on earlier version called HSA.
; 02/08/2008 DGG revised to work with updated SPHEREFIELD syntax.
; 04/15/2008 DGG Added MPP keyword.  Associated documentation fixes.
; 11/10/2012 DGG Updated command syntax. Updated documentation.  Added
;   COMPILE_OPT.
;
; Copyright (c) 2007-2012 David G. Grier
;-

function spheredhmprofile, rho, zp, ap, np, nm, alpha, delta, lambda, mpp

COMPILE_OPT IDL2

npts = n_elements(rho)

x = reform(rho, 1, npts)
y = dblarr(1, npts)

field = spherefield(x, y, zp, ap, np, nm, lambda, mpp, $
                    /cartesian, k = k)
; interference between light scattered by the particle
; and a plane wave polarized along x and propagating along z
field *=  alpha * exp(dcomplex(0, -k*(zp + delta))) ; amplitude and phase factors
field[0, *] += 1.                                   ; \hat{x}

return, total(real_part(field * conj(field)), 1)
end
