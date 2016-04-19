;+
; NAME:
;    sphere_coefficients
;
; PURPOSE:
;    Calculates the Mie scattering coefficients
;    for a multilayered sphere illuminated
;    by a coherent plane wave linearly polarized in the x direction.
;
; CATEGORY:
;    Holography, light scattering, microscopy
;
; CALLING SEQUENCE:
;    ab = sphere_coefficients(ap, np, nm, lambda)
;
; INPUTS:
;    ap: [nlayers] radii of layered sphere [micrometers]
;        NOTE: ap and np are reordered automatically so that
;        ap is in ascending order.
;
;    np: [nlayers] (complex) refractive indexes of sphere's layers
;
;    nm: (complex) refractive index of medium
;
;    lambda: wavelength of light [micrometers]
;
; OUTPUTS:
;    ab: [2,nc] complex a and b scattering coefficients.
;
; KEYWORDS:
;    resolution: minimum magnitude of Lorenz-Mie coefficients to retain.
;         Default: See references.
;
; REFERENCES:
; 1. Adapted from Chapter 8 in
;    C. F. Bohren and D. R. Huffman, 
;    Absorption and Scattering of Light by Small Particles,
;    (New York, Wiley, 1983).
;
; 2. W. Yang, 
;    Improved recursive algorithm for 
;    light scattering by a multilayered sphere,
;    Applied Optics 42, 1710--1720 (2003).
;
; 3. O. Pena, U. Pal,
;    Scattering of electromagnetic radiation by a multilayered sphere,
;    Computer Physics Communications 180, 2348--2354 (2009).
;    NB: Equation numbering follows this reference.
;
; 4. W. J. Wiscombe,
;    Improved Mie scattering algorithms,
;    Applied Optics 19, 1505-1509 (1980).
;
; EXAMPLE: 
;  ab = sphere_coefficients([1.0, 1.5], [1.45, 1.56], 1.3326, 0.6328)
;   
; MODIFICATION HISTORY:
; Replaces sphere_coefficients.pro, which calculated scattering
;    coefficients for a homogeneous sphere.  This earlier version
;    was written by Bo Sun and David G. Grier, and was based on
;    algorithms by W. J. Wiscombe (1980) and W. J. Lentz (1976).
;
; 10/31/2010 Written by F. C. Cheong, New York University
; 11/02/2010 David G. Grier (NYU) Formatting.  Explicit casts 
;    to double precision.  Use complex functions.
; 04/10/2011 DGG Cleaned up Nstop code.  Use array math rather than
;    iteration where convenient.  Eliminate an and bn subroutines.
;    Simplify function names.  Use IDL 8 array notation.  Fixed
;    bug in bn coefficients for multilayer spheres.
; 05/31/2011 DGG Fixed corner condition in Nstop code.
; 09/04/2011 DGG Added RESOLUTION keyword
; 01/26/2013 DGG Calculate x with real part of nm
; 07/22/2013 DGG RESOLUTION retains up to last coefficient
;    with sufficiently large magnitude.
;
; Copyright (c) 2010-2016 F. C. Cheong and D. G. Grier
; 

;;;;;
;
; Nstop
;
; Number of terms to keep in the partial wave expansion
;
function Nstop, x, m

  COMPILE_OPT IDL2, HIDDEN

  ;; Wiscombe (1980)
  xl = x[-1]
  if xl lt 8. then $
     ns = floor(xl + 4. * xl^(1./3.) + 1.) $
  else if xl lt 4200. then $
     ns = floor(xl + 4.05 * xl^(1./3.) + 2.) $
  else if xl gt 4199. then $
     ns = floor(xl + 4. * xl^(1./3.) + 2.)

  ;; Yang (2003) Eq. (30)
  return, double(floor(max([ns, abs(x*m), abs(shift(x, -1)*m)])) + 15)
end

;;;;;
;
; sphere_coefficients
;
function sphere_coefficients, ap, np, nm, lambda, $
                              resolution = resolution

  COMPILE_OPT IDL2

  nlayers = n_elements(ap)

  if n_elements(np) ne nlayers then $
     message, "ap and np must have the same number of elements"

  ;; arrange shells in size order
  if nlayers gt 1 then begin
     order = sort(ap)
     ap = ap[order]
     np = np[order]
  endif

  x = 2.d * !dpi * real_part(nm) * ap / lambda ; size parameter
  m = dcomplex(np/nm)                          ; relative refractive index

  nmax = Nstop(x, m)            ; number of terms in partial-wave expansion

  ci = dcomplex(0, 1)

  ;; storage for results
  ab = dcomplexarr(2, nmax+1, /nozero)

  D1     = dcomplexarr(nmax+2)
  D1_a   = dcomplexarr(nlayers, nmax+2)
  D1_am1 = dcomplexarr(nlayers, nmax+2)

  D3     = dcomplexarr(nmax+1)
  D3_a   = dcomplexarr(nlayers, nmax+1)
  D3_am1 = dcomplexarr(nlayers, nmax+1)
  
  Psi         = dcomplexarr(nmax+1) 
  Zeta        = dcomplexarr(nmax+1) 
  PsiZeta     = dcomplexarr(nmax+1) 
  PsiZeta_a   = dcomplexarr(nlayers, nmax+1) 
  PsiZeta_am1 = dcomplexarr(nlayers, nmax+1) 

  Q  = dcomplexarr(nlayers, nmax+1) 
  Ha = dcomplexarr(nlayers, nmax+1) 
  Hb = dcomplexarr(nlayers, nmax+1) 

  ;;; Iterate outward from the sphere's core
  z1 = x[0] * m[0]

  ; D1_a[0, nmax + 1] = dcomplex(0) ; Eq. (16a)
  for n = nmax + 1, 1, -1 do begin ; downward recurrence
     dn = double(n)
     D1_a[0, n-1] = dn/z1 - 1.d/(D1_a[0, n] + dn/z1) ; Eq. (16b)
  endfor

  PsiZeta_a[0, 0] = 0.5d * (1.d - exp(2.d * ci * z1)) ; Eq. (18a)
  D3_a[0, 0] = ci                                     ; Eq. (18b)
  for n = 1, nmax do begin                            ; upward recurrence
     dn = double(n)
     PsiZeta_a[0, n] = PsiZeta_a[0, n-1] * $ ; Eq. (18c)
                       (dn/z1 - D1_a[0, n-1]) * (dn/z1 - D3_a[0, n-1])
     D3_a[0, n] = D1_a[0, n] + ci/PsiZeta_a[0, n] ; Eq. (18d)
  endfor 

  Ha[0, *] = D1_a[0, 0:-2]      ; Eq. (7a)
  Hb[0, *] = D1_a[0, 0:-2]      ; Eq. (8a)

  ;; Iterate from layer 2 to outermost layer L
  for ii = 1, nlayers - 1 do begin 
     z1 = x[ii] * m[ii]
     z2 = x[ii-1] * m[ii]                                
     for n = nmax + 1, 1, -1 do begin ; Downward recurrence for D1
        dn = double(n)
        D1_a[ii, n-1]   = dn/z1 - 1.d/(D1_a[ii, n]   + dn/z1) ; Eq. (16b)
        D1_am1[ii, n-1] = dn/z2 - 1.d/(D1_am1[ii, n] + dn/z2)
     endfor 
                                
     PsiZeta_a[ii, 0]   = 0.5d * (1.d - exp(2.d * ci * z1)) ; Eq. (18a)
     PsiZeta_am1[ii, 0] = 0.5d * (1.d - exp(2.d * ci * z2))
     D3_a[ii, 0]   = ci         ; Eq. (18b)     
     D3_am1[ii, 0] = ci           
     for n = 1, nmax do begin   ; Upward recurrence for PsiZeta and D3
        dn = double(n)
        PsiZeta_a[ii, n]   = PsiZeta_a[ii, n-1] * $
                             (dn/z1 - D1_a[ii, n-1]) * $
                             (dn/z1 - D3_a[ii, n-1]) ; Eq. (18c)
        PsiZeta_am1[ii, n] = PsiZeta_am1[ii, n-1] * $
                             (dn/z2 - D1_am1[ii, n-1]) * $
                             (dn/z2 - D3_am1[ii, n-1])
        D3_a[ii, n]   = D1_a[ii, n]   + ci/PsiZeta_a[ii, n]
        D3_am1[ii, n] = D1_am1[ii, n] + ci/PsiZeta_am1[ii, n]
     endfor 

     Q[ii, 0] = (exp(-2.d * ci * z2) - 1.d) / $
                (exp(-2.d * ci * z1) - 1.d) ; Eq. (19a)
     for n = 1, nmax do begin               ; Upward recurrence for Q
        dn = double(n)
        Q[ii, n] = Q[ii, n-1] * (x[ii-1]/x[ii])^2 * $
                   (z2 * D1_am1[ii, n] + dn)/(z1 * D1_a[ii, n] + dn) * $
                   (dn - z2 * D3_am1[ii, n])/(dn - z1 * D3_a[ii, n]) ; Eq. (19b)
     endfor 

     for n = 1, nmax do begin   ; Upward recurrence for Ha and Hb
        G1 = m[ii] * Ha[ii-1, n] - m[ii-1] * D1_am1[ii, n] ; Eq. (12)
        G2 = m[ii] * Ha[ii-1, n] - m[ii-1] * D3_am1[ii, n] ; Eq. (13)
        Ha[ii, n] = (G2 * D1_a[ii, n] - G1 * Q[ii, n] * D3_a[ii, n]) / $
                    (G2 - G1 * Q[ii, n]) ; Eq. (7b)

        G1 = m[ii-1] * Hb[ii-1, n] - m[ii] * D1_am1[ii, n] ; Eq. (14)
        G2 = m[ii-1] * Hb[ii-1, n] - m[ii] * D3_am1[ii, n] ; Eq. (15)
        Hb[ii, n] = (G2 * D1_a[ii, n] - G1 * Q1[ii, n] * D3_a[ii, n]) / $
                    (G2 - G1 * Q[ii, n]) ; Eq. (8b)
     endfor
  endfor                        ;ii (layers)

  xL = dcomplex(x[-1])
  for n = nmax, 1, -1 do begin  ; Downward recurrence for D1
     dn = double(n)
     D1[n-1] = dn/xL - (1.d/(D1[n] + dn/xL)) ; Eq. (16b)
  endfor

  Psi[0]     = sin(xL)                           ; Eq. (20a)
  Zeta[0]    = -ci * exp(ci * xL)                ; Eq. (21a)
  PsiZeta[0] = 0.5d * (1.d - exp(2.d * ci * xL)) ; Eq. (18a)
  D3[0] = ci                                     ; Eq. (18b)
  for n = 1, nmax do begin      ; Upward recurrence for Psi, Zeta, PsiZeta and D3
     dn = double(n)
     Psi[n]  = Psi[n-1]  * (dn/xL - D1[n-1])                           ; Eq. (20b)
     Zeta[n] = Zeta[n-1] * (dn/xL - D3[n-1])                           ; Eq. (21b)
     PsiZeta[n] = PsiZeta[n-1] * (dn/xL - D1[n-1]) * (dn/xL - D3[n-1]) ; Eq. (18c)
     D3[n] = D1[n] + ci/PsiZeta[n]                                     ; Eq. (18d)
  endfor

  ;; Scattering coefficients
  n = dindgen(nmax + 1)
  ab[0, *] = ((Ha[-1, *]/m[-1] + n/x[-1]) * Psi  - shift(Psi,  1)) / $
             ((Ha[-1, *]/m[-1] + n/x[-1]) * Zeta - shift(Zeta, 1)) ; Eq. (5)
  ab[1, *] = ((Hb[-1, *]*m[-1] + n/x[-1]) * Psi  - shift(Psi,  1)) / $
             ((Hb[-1, *]*m[-1] + n/x[-1]) * Zeta - shift(Zeta, 1)) ; Eq. (6)
  ab[*, 0] = dcomplex(0)

  if keyword_set(resolution) then begin
     w = where(total(abs(ab), 1) gt resolution, ngood)
     if ngood ge 1 then $
        ab = ab[*,0:w[-1]]
  endif

  return, ab
end
