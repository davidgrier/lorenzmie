;+
; NAME:
;    polarizability
;
; PURPOSE:
;    Calculates the dipole polarizabilities
;    of a sphere illuminated by a monochromatic plane wave.
;
; CATEGORY:
;    Light scattering
;
; CALLING SEQUENCE:
;    alpha = polarizability(ap, np, nm, lambda)
;
; INPUTS:
;    ap: radius of sphere [micrometers]
;    np: (complex) refractive index of sphere
;    nm: (complex) refractive index of medium
;    lambda: wavelength of light [micrometers]
;
; OUTPUTS:
;    alpha: complex electric and magnetic dipole polarizabilities.
;        alpha[0]: electric polarizability
;        alpha[1]: magnetic polarizability
;
; REFERENCES:
; 1. C. F. Bohren and D. R. Huffman,
;    Absorption and Scattering of Light by Small Particles,
;    (New York, Wiley, 1983).
;
; 2. M. I. Mishchenko, L. D. Travis and A. A. Lacis,
;    Scattering, Absorption and Emission of Light by Small Particles,
;    (Cambridge University Press, Cambridge, 2001).
;
; 3. J. Chen, J. Ng, Z. Lin and C. T. Chan,
;    Optical pulling force,
;    Nature Photonics 5, 531-534 (2011).
;
; MODIFICATION HISTORY:
; 08/27/2015 Written by David G. Grier, New York University
;
; Copyright (c) 2015 David G. Grier
;-
function polarizability, ap, np, nm, lambda

  COMPILE_OPT IDL2

  ab = sphere_coefficients(ap, np, nm, lambda)

  k = 2.*!pi*nm/lambda
  alphae = !const.i * 6. * !const.pi * !const.eps0 * nm^2 * ab[0, 0] / k^3
  alpham = !const.i * 6. * !const.pi * ab[1, 0] / (!const.mu0 * k^3)

  return, [alphae, alpham]
end

