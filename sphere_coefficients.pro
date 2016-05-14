; docformat = 'rst'

;+
; Calculate the Lorenz-Mie scattering coefficients
; for a multilayered sphere illuminated
; by a coherent plane wave linearly polarized in the x direction.
;
; :Examples:
;    IDL> ab = sphere_coefficients(ap, np, nm, lambda)
;
;    IDL> ab = sphere_coefficients([1.0, 1.5], [1.45, 1.56], 1.3326, 0.6328)
;
; :Params:
;    ap : in, required, type=`double|doublearr(nlayers)`
;        Radii of sphere's layers [um]
;        NOTE: ap and np are reordered automatically so that
;        ap is in ascending order.
;
;    np : in, required, type=`dcomplex|dcomplexarr(nlayers)`
;        Refractive indexes of sphere's layers
;
;    nm : in, required, type=dcomplex
;        Refractive index of medium
;
;    lambda : in, required, type=double
;        Wavelength of light [micrometers]
;
; :Returns:
;    ab: [2,nc] complex a and b scattering coefficients.
;
; :Keywords:
;    resolution : in, optional, type=double, default=`See references`
;        Minimum magnitude of Lorenz-Mie coefficients to retain.
;
; :References:
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
; 5. A. A. R. Neves and D. Pisignano,
;    Effect of finite terms on the truncation error of Mie series,
;    Optics Letters 37, 2481-2420 (2012).
;
; :History:
; Replaces sphere_coefficients.pro, which calculated scattering
;    coefficients for a homogeneous sphere.  This earlier version
;    was written by Bo Sun and David G. Grier, and was based on
;    algorithms by W. J. Wiscombe (1980) and W. J. Lentz (1976).
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
; 04/19/2016 DGG Substantial overhaul to streamline implementation.
; 04/21/2016 DGG Corrections to b coefficients.  Commit Neves stopping
;    criterion.
;
; :Author:
;    Fook C. Cheong and David G. Grier
;
; :Copyright:
;    Copyright (c) 2010-2016 F. C. Cheong and D. G. Grier
;- 

;+
; Number of terms to keep in the partial wave expansion
; according to Neves and Pisignano (2012), Eq. (9).
;
; :Hidden:
;-
function neves_pisignano, x, epsilon

  COMPILE_OPT IDL2, HIDDEN

  return, floor(x + 0.44*(x * (alog(epsilon))^2)^(1./3) - 4.1)
end

;+
; Number of terms to keep in the partial wave expansion
; according to Wiscombe (1980) and Yang (2003).
;
; :Hidden:
;-
function wiscombe_yang, x, m

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
  return, floor(max([ns, abs(x*m), abs(shift(x, -1)*m)])) + 15
end

;+
; sphere_coefficients
;-
function sphere_coefficients, ap, np, nm, lambda, $
                              resolution = resolution_, $
                              wiscombe = wiscombe_

  COMPILE_OPT IDL2

  nlayers = n_elements(ap)

  if n_elements(np) ne nlayers then $
     message, "ap and np must have the same number of elements"

  resolution = isa(resolution_, /number, /scalar) ? float(resolution_) : 0.
  wiscombe = keyword_set(wiscombe_) || (resolution le 0)

  ;; arrange shells in size order
  if nlayers gt 1 then begin
     order = sort(ap)
     ap = ap[order]
     np = np[order]
  endif

  x = 2.d * !dpi * real_part(nm) * ap / lambda ; size parameter
  m = dcomplex(np/nm)                          ; relative refractive index

  ;; number of terms in partial-wave expansion
  nmax = (wiscombe) ? wiscombe_yang(x, m) : neves_pisignano(x, resolution)

  ci = dcomplex(0, 1)

  ;; storage for results
  ab = dcomplexarr(2, nmax+1, /nozero)

  D1_z1 = dcomplexarr(nmax+2, /nozero)
  D1_z2 = dcomplexarr(nmax+2, /nozero)
  D3_z1 = dcomplexarr(nmax+1, /nozero)
  D3_z2 = dcomplexarr(nmax+1, /nozero)
  Psi   = dcomplexarr(nmax+1, /nozero) 
  Zeta  = dcomplexarr(nmax+1, /nozero) 
  Q = dcomplexarr(nmax+1, /nozero)

  ;; one-time initializations
  D1_z1[-1] = dcomplex(0)       ; Eq. (16a)
  D1_z2[-1] = dcomplex(0)
  D3_z1[0] = ci                 ; Eq. (18b)     
  D3_z2[0] = ci
  
  ;;; Iterate outward from the sphere's core
  ;; Initialize Ha and Hb in core
  z1 = x[0] * m[0]
  for n = nmax+1, 1, -1 do begin ; Downward recurrence for D1
     dn = double(n)
     D1_z1[n-1] = dn/z1 - 1.d/(D1_z1[n] + dn/z1) ; Eq. (16b)
  endfor
  Ha = D1_z1[0:-2]              ; Eq. (7a)
  Hb = D1_z1[0:-2]              ; Eq. (8a)

  ;; Iterate outward from layer 2 to outermost layer L
  for ii = 1, nlayers - 1 do begin 
     z1 = x[ii] * m[ii]
     z2 = x[ii-1] * m[ii]
     
     for n = nmax+1, 1, -1 do begin ; Downward recurrence for D1
        dn = double(n)
        D1_z1[n-1] = dn/z1 - 1.d/(D1_z1[n] + dn/z1) ; Eq. (16b)
        D1_z2[n-1] = dn/z2 - 1.d/(D1_z2[n] + dn/z2)
     endfor 
                                
     PsiZeta_z1 = 0.5d * (1.d - exp(2.d * ci * z1)) ; Eq. (18a)
     PsiZeta_z2 = 0.5d * (1.d - exp(2.d * ci * z2))
     Q[0] = (exp(-2.d * ci * z2) - 1.d) / $
            (exp(-2.d * ci * z1) - 1.d) ; Eq. (19a)
     for n = 1, nmax do begin   ; Upward recurrence for PsiZeta, D3, Q
        dn = double(n)
        PsiZeta_z1 *= (dn/z1 - D1_z1[n-1]) * (dn/z1 - D3_z1[n-1]) ; Eq. (18c)
        PsiZeta_z2 *= (dn/z2 - D1_z2[n-1]) * (dn/z2 - D3_z2[n-1])
        D3_z1[n] = D1_z1[n] + ci/PsiZeta_z1 ; Eq. (18d)
        D3_z2[n] = D1_z2[n] + ci/PsiZeta_z2
        Q[n] = Q[n-1] * (x[ii-1]/x[ii])^2 * $
               (z2 * D1_z2[n] + dn)/(z1 * D1_z1[n] + dn) * $
               (dn - z2 * D3_z2[n])/(dn - z1 * D3_z1[n]) ; Eq. (19b)
     endfor

     ;; Update Ha and Hb
     G1 = m[ii] * Ha - m[ii-1] * D1_z2                  ; Eq. (12)
     G2 = m[ii] * Ha - m[ii-1] * D3_z2                  ; Eq. (13)
     Ha = (G2 * D1_z1 - Q * G1 * D3_z1) / (G2 - G1 * Q) ; Eq. (7b)
     
     G1 = m[ii-1] * Hb - m[ii] * D1_z2                  ; Eq. (14)
     G2 = m[ii-1] * Hb - m[ii] * D3_z2                  ; Eq. (15)
     Hb = (G2 * D1_z1 - Q * G1 * D3_z1) / (G2 - G1 * Q) ; Eq. (8b)
  endfor

  ;; Iterate into medium
  z1 = dcomplex(x[-1])
  for n = nmax+1, 1, -1 do begin ; Downward recurrence for D1
     dn = double(n)
     D1_z1[n-1] = dn/z1 - (1.d/(D1_z1[n] + dn/z1)) ; Eq. (16b)
  endfor

  Psi[0]  = sin(z1)                           ; Eq. (20a)
  Zeta[0] = -ci * exp(ci * z1)                ; Eq. (21a)
  PsiZeta = 0.5d * (1.d - exp(2.d * ci * z1)) ; Eq. (18a)
  for n = 1, nmax do begin ; Upward recurrence for Psi, Zeta and D3
     dn = double(n)
     Psi[n]  = Psi[n-1]  * (dn/z1 - D1_z1[n-1])             ; Eq. (20b)
     Zeta[n] = Zeta[n-1] * (dn/z1 - D3_z1[n-1])             ; Eq. (21b)
     PsiZeta *= (dn/z1 - D1_z1[n-1]) * (dn/z1 - D3_z1[n-1]) ; Eq. (18c)
     D3_z1[n] = D1_z1[n] + ci/PsiZeta                       ; Eq. (18d)
  endfor

  ;; Scattering coefficients
  n = dindgen(nmax + 1)
  fac = Ha/m[-1] + n/x[-1]
  ab[0, *] = (fac * Psi  - shift(Psi,  1)) / $
             (fac * Zeta - shift(Zeta, 1)) ; Eq. (5)
  
  fac = m[-1]*Hb + n/x[-1]
  ab[1, *] = (fac * Psi  - shift(Psi,  1)) / $
             (fac * Zeta - shift(Zeta, 1)) ; Eq. (6)
  ab[*, 0] = dcomplex(0)

  if (wiscombe) && (resolution gt 0.) then begin
     w = where(total(abs(ab), 1) gt resolution, ngood)
     if ngood ge 1 then $
        ab = ab[*,0:w[-1]]
  endif

  return, ab
end
