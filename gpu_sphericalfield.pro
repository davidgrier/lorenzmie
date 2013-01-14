;+
; NAME:
;    gpu_sphericalfield
;
; PURPOSE:
;    Calculates the electric field in a light scattering
;    pattern defined by a set of Lorenz-Mie scattering
;    coefficients.  Uses gpulib for hardware acceleration.
;
; CATEGORY:
;    Light scattering, holographic microscopy
;
; CALLING SEQUENCE:
;    field = gpu_sphericalfield(x, y, z, ab, lambda, k = k)
;
; INPUTS:
;    x: [npts] array of pixel coordinates [pixels]
;    y: [npts] array of pixel coordinates [pixels]
;    z: If field is required in a single plane, then
;       z is the plane's distance from the sphere's center
;       [pixels].
;       Otherwise, z is an [npts] array of coordinates.
;
;    ab: [2,nc] array of a and b scattering coefficients, where
;        nc is the number of terms required for convergence.
;
;    lambda: wavelength of light in medium [pixels]
;
; KEYWORD FLAGS:
;    cartesian: If set, return field components in Cartesian
;               coordinates.  Default: Spherical polar coordinates
;
; OUTPUT KEYWORDS:
;    k: scaled wavenumber of light in medium [inverse pixels]
;
; OUTPUTS:
;    field: [3,npts] complex values of field at the positions r.
;           [0,*]: r component
;           [1,*]: theta component
;           [2,*]: phi component
;
;           If CARTESIAN is set:
;           [0,*]: x component (incident polarization)
;           [1,*]: y component (transverse component)
;           [2,*]: z component (axial component, relative to
;           incident beam).
;
; REFERENCES:
; 1. Adapted from Chapter 4 in
;    C. F. Bohren and D. R. Huffman, 
;    Absorption and Scattering of Light by Small Particles,
;    (New York, Wiley, 1983).
;
; 2. W. J. Wiscombe, "Improved Mie scattering algorithms,"
;    Appl. Opt. 19, 1505-1509 (1980).
;
; 3. W. J. Lentz, "Generating Bessel function in Mie scattering
;    calculations using continued fractions," Appl. Opt. 15,
;    668-671 (1976).
;
; 4. S. H. Lee, Y. Roichman, G. R. Yi, S. H. Kim, S. M. Yang,
;    A. van Blaaderen, P. van Oostrum and D. G. Grier,
;    "Characterizing and tracking single colloidal particles with
;    video holographic microscopy," Opt. Express 15, 18275-18282
;    (2007).
;
; 5. GPULIB: http://www.txcorp.com/products/GPULib/
;
; 6. F. C. Cheong, B. Sun, R. Dreyfus, J. Amato-Grill, K. Xiao,
;    L. Dixon and D. G. Grier,
;    "Flow visualization and flow cytometry with holographic video
;    microscopy," Opt. Express 17, 13071-13079 (2009).
;
; 7. P. Messmer, P. J. Mullowney and B. E. Granger,
;    "GPULib: GPU computing in high-level languages,"
;    Computer Science and Engineering 10, 70-73 (2008).
;
; MODIFICATION HISTORY:
; Written by David G. Grier, New York University, 5/2007
; 06/09/2007 DGG finally read Section 4.8 in Bohren and Huffman about
;    numerical stability of the recursions used to compute the scattering
;    coefficients.  Feh.  Result is a total rewrite.
; 06/20/2007 DGG Calculate \tau_n(\cos\theta) and \pi_n(\cos\theta)
;    according to recurrence relations in 
;    W. J. Wiscombe, Appl. Opt. 19, 1505-1509 (1980).
;    This is supposed to improve numerical accuracy.
; 02/08/2008 DGG. Replaced single [3,npts] array of input coordinates
;    with two [npts] arrays for x and y, and a separate input for z.
;    Eliminated double() call for coordinates.  Z may have 1 element or
;    npts elements. Small documentation fixes.
; 04/03/2008 Bo Sun (Sephiroth), NYU Calculate Lorenz-Mie a and b
;    coefficients using continued fractions rather than recursion.
;    Osman Akcakir from Arryx pointed out that the results are
;    more accurate in extreme cases.  Method described in
;    William J. Lentz, "Generating Bessel functions in Mie scattering
;    calculations using continued fractions," Appl. Opt. 15, 668-671
;    (1976).
; 04/04/2008 DGG small code clean-ups and documentation.  Added
;    RECURSIVE keyword for backward compatibility in computing a and b
;    coefficients.
; 04/11/2008 Sephiroth Corrected small error in jump code for
;    repeated fractions in Mie coefficients.
; 06/25/2008 DGG Don't clobber x input coordinates.
; 09/26/2008 Rewritten for use with CUDA via GPULIB by DGG and
;     Jesse Amato-Grill.
; 10/09/2008 Bo Sun (Sephiroth) Rewritten to reduce all complex
;     operation into real ones
; 10/11/2008 Bo Sun (Sephiroth) Debug and fixed VRAM leaking.
; 10/13/2008 DGG adapted from GPU_SPHEREFIELD to conform to updated
;     DHM code layout.  Eliminating temporary GPU variables increases
;     speed by nearly factor of 2 :).
;     Eliminated RECURSIVE keyword because scattering coefficients are
;     now calculated elsewhere.
; 02/05/2009 DGG convert to low-level GPUlib API calls for improved speed.
;     Explicitly cast constants to single-precision.
; 02/07/2009 DGG eliminate extraneous calls to
;     cudaThreadSynchronize(), which appear to be very costly.
;     Report CUDA errors, which are reported by cudaThreadSynchronize().
; 08/23/2009 DGG updated references.
; 06/18/2010 DGG Added COMPILE_OPT
; 10/18/2010 DGG Updated for GPUlib 1.4.0 with IDL 8.0 syntax.
;     Initial use of complex GPU routines, which greatly reduces code size
;     and increases speed.  Reverted to high-level GPUlib calls for readability,
;     version compatibility, and to simplify support for double precision.
; 10/19/2010 DGG Experimental support for double precision.  Fixed
;     projection bug in Es1 that affected only GPU code.  GPU and CPU
;     versions now agree to better than 3e-5 in single precision.
;     Added reference to GPUlib publication.
; 06/28/2011 DGG identified bug: returns NAN when x = y = 0.  Problem
;     line identified by FIXME below.  Current workaround: make sure
;     not to call with x = y = 0.  Even a small offset will work.
; 01/14/2013 DGG documentation fixes.
;
; Copyright (c) 2007-2013 Bo Sun, Jesse Amato-Grill and David G. Grier
;-

function gpu_sphericalfield, x_, y_, z_, ab_, lambda, $
                             cartesian = cartesian, $ ; project to cartesian coordinates
                             double = double          ; use double precision if available

COMPILE_OPT IDL2

type = (keyword_set(double) && gpudoublecapable()) ? 9 : 6 ; dcomplex or complex

ab = fix(ab_, TYPE = type)      ; Lorenz-Mie scattering coefficients
nc = n_elements(ab[0, *]) - 1   ; number of terms required for convergence
ci = fix(complex(0, 1), TYPE = type)

; convert to spherical coordinates centered on the sphere.
; (r, theta, phi) is the spherical coordinate of the pixel
; at (x,y) in the imaging plane at axial distance z from the
; center of the sphere.

; transfer coordinates to GPU
x = gpuPutArr(x_)
y = gpuPutArr(y_)
x = gpuFix(x, TYPE = type, /OVERWRITE)
y = gpuFix(y, TYPE = type, /OVERWRITE)
npts = x.n_elements
if n_elements(z_) eq 1 then begin
   z = gpuMake_Array(npts, VALUE = z_, TYPE = type)
endif else begin
   z = gpuPutArr(z_)
   z = gpuFix(z, TYPE = type, /OVERWRITE)
endelse

if y.n_elements ne npts or z.n_elements ne npts then $
   message, "dimensions of x, y and z must agree."

; rho   = sqrt(x^2 + y^2)
rho = gpuMult(x, x, /NONBLOCKING)
kr  = gpuMult(y, y)
rho = gpuAdd(rho, kr, LHS = rho, /NONBLOCKING)

; r = sqrt(rho^2 + z^2)
kr = gpuMult(z, z, LHS = kr)
kr = gpuAdd(kr, rho, LHS = kr)

rho = gpuSqrt(rho, LHS = rho, /NONBLOCKING)
kr = gpuSqrt(kr, LHS = kr)

; polar angle
sintheta = gpuDiv(rho, kr, /NONBLOCKING)
costheta = gpuDiv(z, kr, /NONBLOCKING)

; azimuthal angle
cosphi = gpuDiv(x, rho, /NONBLOCKING)
sinphi = gpuDiv(y, rho, /NONBLOCKING)
;;; NOTE: this results in a divide-by-zero error when x = y = 0.
;;; Have not yet found an efficient fix
;;; The following fails because gpuAtan2 does not work for complex
;;; data types.  Addressing this by making some variables DOUBLE
;;; rather than DCOMPLEX slows the routine by a factor of 3 because
;;; of all of the additional type conversions this requires.
;cosphi = gpuAtan2(y, x)
;sinphi = gpuSin(cosphi)
;cosphi = gpuCos(cosphi, LHS = cosphi, /NONBLOCKING)

; scale distances by wavenumber
kr = gpuAdd(2.*!pi/lambda, kr, 0., kr, 0., LHS = kr)

; starting points for recursive function evaluation ...
; ... Riccati-Bessel radial functions, page 478

; sinkr = sin(kr)
; coskr = cos(kr)
; xi_nm2 = dcomplex(coskr, sinkr) ; \xi_{-1}(kr)
; xi_nm1 = dcomplex(sinkr,-coskr) ; \xi_0(kr)
x = gpuCos(kr, LHS = x, /NONBLOCKING)
y = gpuSin(kr, LHS = y)
xi_nm2 = (type eq 6) ? gpuComplex(x, y) : gpuDcomplex(x, y)
xi_nm1 = gpuAdd(-ci, xi_nm2, 0., xi_nm2, 0., /NONBLOCKING)
xi_n   = gpuMake_Array(npts, /NOZERO, TYPE = type)

; ... angular functions (4.47), page 95
;pi_nm1 = 0.d                    ; \pi_0(\cos\theta)
;pi_n   = 1.d                    ; \pi_1(\cos\theta)
pi_nm1 = gpuMake_Array(npts, TYPE = type)
pi_n   = gpuMake_Array(npts, VALUE = 1., TYPE = type)

; storage for vector spherical harmonics: [r,theta,phi]
;Mo1n = dcomplexarr(3,npts)
;Ne1n = dcomplexarr(3,npts)
Mo1n2 = gpuMake_Array(npts, /NOZERO, TYPE = type)
Mo1n3 = gpuMake_Array(npts, /NOZERO, TYPE = type)
Ne1n1 = gpuMake_Array(npts, /NOZERO, TYPE = type)
Ne1n2 = gpuMake_Array(npts, /NOZERO, TYPE = type)
Ne1n3 = gpuMake_Array(npts, /NOZERO, TYPE = type)

; storage for Deirmendjian's derivative
dn = gpuMake_Array(npts, /NOZERO, TYPE = type)

; Storage for accumulating scattered field, Es
;Es = dcomplexarr(3,npts)
Es1 = gpuMake_Array(npts, TYPE = type) ; x/r component
Es2 = gpuMake_Array(npts, TYPE = type) ; y/theta component
Es3 = gpuMake_Array(npts, TYPE = type) ; z/phi component

; variables
swisc = gpuMake_Array(npts, /NOZERO, TYPE = type)
twisc = gpuMake_Array(npts, /NOZERO, TYPE = type)
tau_n = gpuMake_Array(npts, /NOZERO, TYPE = type)

; Compute field by summing multipole contributions
for n = 1., nc do begin

; upward recurrences ...
; ... Legendre factor (4.47)
; Method described by Wiscombe (1980)
;    swisc = pi_n * costheta 
;    twisc = swisc - pi_nm1
;    tau_n = n * twisc - pi_nm1  ; \tau_n(\cos\theta)
   swisc = gpuMult(pi_n, costheta, LHS = swisc)
   twisc = gpuSub(swisc, pi_nm1, LHS = twisc)
   tau_n = gpuSub(1., pi_nm1, n, twisc, 0., LHS = tau_n, /NONBLOCKING)
                                ; NOTE: tau_n appears as -tau_n in the
                                ; scattering formulas -- define accordingly here.

; ... Riccati-Bessel function, page 478
;    xi_n   = (2.d*n - 1.d) * xi_nm1 / kr - xi_nm2    ; \xi_n(kr)
   xi_n = gpuDiv(xi_nm1, kr, LHS = xi_n)
   xi_n = gpuAdd(2.*n-1., xi_n, -1., xi_nm2, 0., LHS = xi_n)

; ... Deirmendjian's derivative
;    dn = (n * xi_n)/kr - xi_nm1
   dn = gpuDiv(xi_n, kr, LHS = dn)
   dn = gpuSub(n, dn, 1., xi_nm1, 0., LHS = dn)
   
; vector spherical harmonics (4.50)
;    Mo1n[0,*] = 0.d             ; no radial component (as initialized)
;    Mo1n[1,*] = pi_n * xi_n     ; ... divided by cosphi/kr
;    Mo1n[2,]* = -tau_n * xi_n   ; ... divided by sinphi/kr
   Mo1n2 = gpuMult(pi_n, xi_n, LHS = Mo1n2, /NONBLOCKING)
   Mo1n3 = gpuMult(tau_n, xi_n, LHS = Mo1n3, /NONBLOCKING)

;    Ne1n[0,*] = n*(n + 1.d) * pi_n * xi_n ; ... divided by cosphi sintheta/kr^2
;    Ne1n[1,*] = -tau_n * dn     ; ... divided by cosphi/kr
;    Ne1n[2,*] = pi_n  * dn      ; ... divided by sinphi/kr
   Ne1n1 = gpuMult(n*(n+1.), pi_n, 1., xi_n, 0., LHS = Ne1n1, /NONBLOCKING)
   Ne1n2 = gpuMult(tau_n, dn, LHS = Ne1n2, /NONBLOCKING)
   Ne1n3 = gpuMult(pi_n, dn, LHS = Ne1n3)

; upward recurrences ...
; ... angular functions (4.47)
; Method described by Wiscombe (1980)
;    pi_nm1 = pi_n
;    pi_n = swisc + (n + 1.d) * twisc / n
   temp = pi_nm1
   pi_nm1 = pi_n
   pi_n = temp
   pi_n = gpuAdd(1., swisc, (n+1.)/n, twisc, 0., LHS = pi_n, /NONBLOCKING)

; ... Riccati-Bessel function
;    xi_nm2 = xi_nm1
;    xi_nm1 = xi_n
   temp = xi_nm2
   xi_nm2 = xi_nm1
   xi_nm1 = xi_n
   xi_n = temp

; Field calculation
; prefactor, page 93
   En = ci^n * (2.*n + 1.)/ n / (n + 1.)
   an = En * ci * ab[0,n]
   bn = En * ab[1,n]
;   if check_math() ne 0 then print, n, ab[*, n], En
; NOTE: can cause floating underflow errors in single precision

; the scattered field in spherical coordinates (4.45)
;    Es += En * (ci * ab[0,n] * Ne1n - ab[1,n] * Mo1n)
   Es1 = gpuAdd(1., Es1, an, Ne1n1, 0., LHS = Es1, /NONBLOCKING)

   Es2 = gpuAdd(1., Es2, an, Ne1n2, 0., LHS = Es2)
   Es2 = gpuSub(1., Es2, bn, Mo1n2, 0., LHS = Es2, /NONBLOCKING)

   Es3 = gpuAdd(1., Es3, an, Ne1n3, 0., LHS = Es3)
   Es3 = gpuSub(1., Es3, bn, Mo1n3, 0., LHS = Es3, /NONBLOCKING)
endfor

; geometric factors were divided out of the vector
; spherical harmonics for accuracy and efficiency ...
; ... put them back at the end.
;Es[0,*] *= cosphi * sintheta / kr^2
;Es[1,*] *= cosphi / kr
;Es[2,*] *= sinphi / kr

x = gpuDiv(cosphi, kr, LHS = x, /NONBLOCKING)
y = gpuDiv(sintheta, kr, LHS = y, /NONBLOCKING)
z = gpuDiv(sinphi, kr, LHS = z)
Es1 = gpuMult(Es1, x, LHS = Es1)
Es1 = gpuMult(Es1, y, LHS = Es1, /NONBLOCKING)
Es2 = gpuMult(Es2, x, LHS = Es2, /NONBLOCKING)
Es3 = gpuMult(Es3, z, LHS = Es3, /NONBLOCKING)

; By default, the scattered wave is returned in spherical
; coordinates.  Project components onto Cartesian coordinates.
; Assumes that the incident wave propagates along z and 
; is linearly polarized along x
if keyword_set(cartesian) then begin
;    Ec = Es
   Ec = xi_n                    ; temporary storage; rename for clarity

;    Ec[0,*] =  Es[0,*] * sintheta * cosphi
;    Ec[0,*] += Es[1,*] * costheta * cosphi
;    Ec[0,*] -= Es[2,*] * sinphi
   x  = gpuMult(cosphi, sintheta, LHS = x, /NONBLOCKING)
   y  = gpuMult(cosphi, costheta, LHS = y)
   Ec = gpuMult(Es1, x, LHS = Ec, /NONBLOCKING)
   Mo1n2  = gpuMult(Es2, y, LHS = Mo1n2)
   Ec = gpuAdd(Ec, Mo1n2, LHS = Ec, /NONBLOCKING)
   Mo1n3  = gpuMult(Es3, sinphi, LHS = Mo1n3)
   Ec = gpuSub(Ec, Mo1n3, LHS = Ec)

   Es1_cpu = gpuGetArr(Ec)

;    Ec[1,*] =  Es[0,*] * sintheta * sinphi
;    Ec[1,*] += Es[1,*] * costheta * sinphi
;    Ec[1,*] += Es[2,*] * cosphi
   x  = gpuMult(sinphi, sintheta, LHS = x, /NONBLOCKING)
   y  = gpuMult(sinphi, costheta, LHS = y)
   Ec = gpuMult(Es1, x, LHS = Ec, /NONBLOCKING)
   Mo1n2  = gpuMult(Es2, y, LHS = Mo1n2)
   Ec = gpuAdd(Ec, z, LHS = Ec, /NONBLOCKING)
   Mo1n3  = gpuMult(Es3, cosphi, LHS = Mo1n3)
   Ec = gpuAdd(Ec, Mo1n3, LHS = Ec)

   Es2_cpu = gpuGetArr(Ec)
   
;    Ec[2,*] =  Es[0,*] * costheta - Es[1,*] * sintheta
   Ec = gpuMult(Es1, costheta, LHS = Ec, /NONBLOCKING)
   Mo1n2 = gpuMult(Es2, sintheta, LHS = Mo1n2)
   Ec = gpuSub(Ec, Mo1n2, LHS = Ec)

   Es3_cpu = gpuGetArr(Ec)
endif else begin
   Es1_cpu = gpuGetArr(Es1)
   Es2_cpu = gpuGetArr(Es2)
   Es3_cpu = gpuGetArr(Es3)
endelse

; GPU objects are cleaned up automatically
; when they go out of scope, starting with GPUlib 1.4.0
;gpuFree, [x, y, z]
;gpuFree, [rho, kr, costheta, sintheta, cosphi, sinphi]
;gpuFree, [swisc, twisc, tau_n]
;gpuFree, [xi_n, xi_nm1, xi_nm2]
;gpuFree, [pi_n, pi_nm1, dn]
;gpuFree, [Mo1n2, Mo1n3]
;gpuFree, [Ne1n1, Ne1n2, Ne1n3]
;gpuFree, [Es1, Es2, Es3]

Es = [[Es1_cpu], [Es2_cpu], [Es3_cpu]]

return, transpose(Es)
end
