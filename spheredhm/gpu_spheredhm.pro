;+
; NAME:
;    gpu_spheredhm
;
; PURPOSE:
;    Calculates the in-line holographic image of a sphere
;    using GPULIB for hardware acceleration
;
; CATEGORY:
;    Holographic microscopy
;
; CALLING SEQUENCE:
;    a = gpu_spheredhm(rp, ap, np, nm, dim)
;
; INPUTS:
;    rp: [xp, yp, zp] position of sphere relative to center of
;       image.
;
;    ap: radius of sphere [micrometer]
;
;    np: refractive index of sphere
;
;    nm: refractive index of medium
;
;    dim: [nx, ny] dimensions of image [pixels]
;
; KEYWORD PARAMETERS:
;    alpha: fraction of incident light scattered by particle.
;
;    lambda: vacuum wavelength of light [micrometers]
;
;    mpp: micrometers per pixel
;
;    precision: relative precision with which fields are calculated.       
;
; OUTPUTS:
;    a: [nx,ny] holographic microscopy image   
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
; MODIFICATION HISTORY:
; Written by David G. Grier, New York University, 5/2007
; 6/9/2007: DGG finally read Section 4.8 in Bohren and Huffman about
;    numerical stability of the recursions used to compute the scattering
;    coefficients.  Feh.  Result is a total rewrite.
; 6/20/2007: DGG Calculate \tau_n(\cos\theta) and \pi_n(\cos\theta)
;    according to recurrence relations in 
;    W. J. Wiscombe, Appl. Opt. 19, 1505-1509 (1980).
;    This is supposed to improve numerical accuracy.
; 2/8/2008: DGG. Replaced single [3,npts] array of input coordinates
;    with two [npts] arrays for x and y, and a separate input for z.
;    Eliminated double() call for coordinates.  Z may have 1 element or
;    npts elements. Small documentation fixes.
; 4/3/2008: Bo Sun (Sephiroth), NYU: Calculate Lorenz-Mie a and b
;    coefficients using continued fractions rather than recursion.
;    Osman Akcakir from Arryx pointed out that the results are
;    more accurate in extreme cases.  Method described in
;    William J. Lentz, "Generating Bessel functions in Mie scattering
;    calculations using continued fractions," Appl. Opt. 15, 668-671
;    (1976).
; 4/4/2008: DGG small code clean-ups and documentation.  Added
;    RECURSIVE keyword for backward compatibility in computing a and b
;    coefficients.
; 4/11/2008: Sephiroth: Corrected small error in jump code for
;    repeated fractions in Mie coefficients.
; 6/25/2008: DGG Don't clobber x input coordinates.
; 9/26/2008: Rewritten for use with cuda via gpulib by DGG and JAG
; 10/9/2008: Bo Sun (Sephiroth) Rewritten to reduce all complex
;            operation into real ones
; 10/11/2008: Bo Sun (Sephiroth) Debug and fixed VRAM leaking
; 10/13/2008: DGG adapted from GPU_SPHEREFIELD to conform to new
;     code layout.  Eliminating temporary GPU variables increases
;     speed by nearly factor of 2 :).
;     Eliminated RECURSIVE keyword.
; 02/05/2009: DGG use low-level gpulib API calls for improved speed.
;     NOTE: Requires cudaThreadSynchronize() calls to prevent random
;     segmentation faults.  Explicitly convert constants to
;     single-precision.
; 04/08/2009: DGG combined functions of GPU_SPHERICALFIELD,
;     SPHEREFIELD, and SPHEREDHM to make GPU acceleration for entire
;     DHM calculation.  NOTE: npts really does have to be cast LONG.
; 04/09/2009: DGG eliminated GPU variables for costheta and z, which
;     are not needed for DHM calculations on plane.
;     Reorganized main loop to eliminate separate storage for
;     SWISC, TWISC and TAU_N.
; 04/10/2009: DGG cleaned up calls to cudaSynchronizeThreads().
;     Reorganized use of temporary variables to eliminate some
;     synchronizations.  Compute pixel coordinates X and Y on GPU.
; 06/18/2010: DGG Added COMPILE_OPT.
;
; Copyright (c) 2007-2010 Bo Sun and David G. Grier.
;-

function gpu_spheredhm, rp, ap, np, nm, dim, $
                        alpha = alpha, $
                        lambda = lambda, $
                        mpp = mpp, $
                        precision = precision

COMPILE_OPT IDL2

if n_elements(lambda) eq 0 then $
   lambda = 0.6328              ; wavelength in vacuum [micrometers]

if n_elements(mpp) eq 0 then $
   mpp = 0.101                  ; length calibration [micrometers/pixel]

; Lorenz-Mie scattering coefficients
ab = sphere_coefficients(ap, np, nm, lambda)
if n_elements(precision) eq 1 then begin
   fac = total(abs(ab[1,*]), 1)
   w = where(fac gt precision*max(fac))
   ab = ab[*,[0,w]]             ; retain first coefficient for bookkeeping
endif
ab = complex(ab)
nc = n_elements(ab[0,*])-1      ; number of terms required for convergence

; scaled wavenumber of light in medium
lambda_m = float(lambda / real_part(nm) / mpp) ; medium wavelength [pixel]
k = 2. * !pi / lambda_m                        ; [pixel^{-1}]

ci = complex(0,1)

err = 0                         ; gpulib error code

;;;;;;
;
; Cartesian coordinates of pixels
;
nx = float(dim[0])
ny = float(dim[1])
npts = long(nx * ny)

x = gpuMake_Array(npts, /INDEX, /FLOAT, ERROR=e) & err OR= e
y = gpuMake_Array(npts, /INDEX, /FLOAT, ERROR=e) & err OR= e
; NOTE: x and y are used here for pixel coordinates and are
;     reused in the main loop as temporary storage

e = gpuFloorFAT(npts, 1., 1./nx, y.handle, 0., 0., y.handle)
e = gpuSubFAT(npts, 1., x.handle, nx, y.handle, rp[0] - 0.5*nx, x.handle)
e = gpuAddFAT(npts, 1., y.handle, 0., y.handle, rp[1] - 0.5*ny, y.handle)
err OR= cudaThreadSynchronize()

z = float(rp[2])                ; scalar: only compute on DHM plane
kz = k*z

;;;;;;
;
; Convert to spherical coordinates centered on the sphere.
; (r, theta, phi) is the spherical coordinate of the pixel
; at (x,y) in the imaging plane at distance z from the
; center of the sphere.

; rho   = sqrt(x^2 + y^2)
rho = gpuMake_Array(npts, /NOZERO, /FLOAT, ERROR=e) & err OR= e
kr  = gpuMake_Array(npts, /NOZERO, /FLOAT, ERROR=e) & err OR= e

e = gpuMultF(npts, x.handle, x.handle, rho.handle)   ; x^2
e = gpuMultF(npts, y.handle, y.handle, kr.handle)    ; y^2
err OR= cudaThreadSynchronize()
e = gpuAddF(npts, rho.handle, kr.handle, rho.handle) ; x^2 + y^2
err OR= cudaThreadSynchronize()

; r = sqrt(rho^2 + z^2)
e = gpuSqrtFAT(npts, 1., 1., rho.handle, z^2, 0., kr.handle)
err OR= cudaThreadSynchronize()

e = gpuSqrtF(npts, rho.handle, rho.handle) ; rho = sqrt(x^2 + y^2) 
err OR= cudaThreadSynchronize()

; polar angle
; NOTE: costheta = z/r -- use gpuDiv instead of defining here
sintheta = gpuMake_Array(npts, /NOZERO, /FLOAT, ERROR=e) & err OR= e
e = gpuDivF(npts, rho.handle, kr.handle, sintheta.handle)

; azimuthal angle
cosphi = gpuMake_Array(npts, /NOZERO, /FLOAT, ERROR=e) & err OR= e
sinphi = gpuMake_Array(npts, /NOZERO, /FLOAT, ERROR=e) & err OR= e
e = gpuDivF(npts, x.handle, rho.handle, cosphi.handle)
e = gpuDivF(npts, y.handle, rho.handle, sinphi.handle)

; scale distances by wavenumber
e = gpuAddFAT(npts, k, kr.handle, 0., kr.handle, 0., kr.handle)
err OR= cudaThreadSynchronize()

;;;;;;
;
; Starting points for recursive function evaluation ...
;
; ... Riccati-Bessel radial functions, page 478
; sinkr = sin(kr)
; coskr = cos(kr)
; xi_nm2 = dcomplex(coskr, sinkr) ; \xi_{-1}(kr)
; xi_nm1 = dcomplex(sinkr,-coskr) ; \xi_0(kr)
xi_nm2_r = gpuMake_Array(npts, /NOZERO, /FLOAT, ERROR=e) & err OR= e
xi_nm2_i = gpuMake_Array(npts, /NOZERO, /FLOAT, ERROR=e) & err OR= e
xi_nm1_r = gpuMake_Array(npts, /NOZERO, /FLOAT, ERROR=e) & err OR= e
xi_nm1_i = gpuMake_Array(npts, /NOZERO, /FLOAT, ERROR=e) & err OR= e
xi_n_r   = gpuMake_Array(npts, /NOZERO, /FLOAT, ERROR=e) & err OR= e
xi_n_i   = gpuMake_Array(npts, /NOZERO, /FLOAT, ERROR=e) & err OR= e

e = gpuCosF(npts, kr.handle, xi_nm2_r.handle)
e = gpuSinF(npts, kr.handle, xi_nm2_i.handle)
e = gpuSinF(npts, kr.handle, xi_nm1_r.handle)
e = gpuCosFAT(npts, -1., 1., kr.handle, 0., 0., xi_nm1_i.handle)

; ... angular functions (4.47), page 95
;pi_nm1 = 0.d                    ; \pi_0(\cos\theta)
;pi_n   = 1.d                    ; \pi_1(\cos\theta)
pi_nm1 = gpuMake_Array(npts, /FLOAT, ERROR=e) & err OR= e
pi_n = gpuMake_Array(npts, VALUE=1., /FLOAT, ERROR=e) & err OR= e

; ... vector spherical harmonics: [r,theta,phi]
;Mo1n = dcomplexarr(3,npts)
;Ne1n = dcomplexarr(3,npts)
Mo1n2_r = gpuMake_Array(npts, /NOZERO, /FLOAT, ERROR=e) & err OR= e
Mo1n2_i = gpuMake_Array(npts, /NOZERO, /FLOAT, ERROR=e) & err OR= e
Mo1n3_r = gpuMake_Array(npts, /NOZERO, /FLOAT, ERROR=e) & err OR= e
Mo1n3_i = gpuMake_Array(npts, /NOZERO, /FLOAT, ERROR=e) & err OR= e
Ne1n1_r = gpuMake_Array(npts, /NOZERO, /FLOAT, ERROR=e) & err OR= e
Ne1n1_i = gpuMake_Array(npts, /NOZERO, /FLOAT, ERROR=e) & err OR= e
Ne1n2_r = gpuMake_Array(npts, /NOZERO, /FLOAT, ERROR=e) & err OR= e
Ne1n2_i = gpuMake_Array(npts, /NOZERO, /FLOAT, ERROR=e) & err OR= e
Ne1n3_r = gpuMake_Array(npts, /NOZERO, /FLOAT, ERROR=e) & err OR= e
Ne1n3_i = gpuMake_Array(npts, /NOZERO, /FLOAT, ERROR=e) & err OR= e

; Storage for Dermidjian's derivative
dn_r = gpuMake_Array(npts, /NOZERO, /FLOAT, ERROR=e) & err OR= e
dn_i = gpuMake_Array(npts, /NOZERO, /FLOAT, ERROR=e) & err OR= e

; Storage for scattered field, Es
;Es = dcomplexarr(3,npts)
Es1_r = gpuMake_Array(npts, /FLOAT, ERROR=e) & err OR= e ; x/r component
Es1_i = gpuMake_Array(npts, /FLOAT, ERROR=e) & err OR= e
Es2_r = gpuMake_Array(npts, /FLOAT, ERROR=e) & err OR= e ; y/theta component
Es2_i = gpuMake_Array(npts, /FLOAT, ERROR=e) & err OR= e
Es3_r = gpuMake_Array(npts, /FLOAT, ERROR=e) & err OR= e ; z/phi component
Es3_i = gpuMake_Array(npts, /FLOAT, ERROR=e) & err OR= e

; reuse GPU variables
swisc = x
twisc = y
tau_n = rho

;;;;;;
;
; Compute field by summing multipole contributions
;
for n = 1., nc do begin

; upward recurrences ...
; ... Legendre factor (4.47):  \tau_n(\cos\theta)
; Method described by Wiscombe (1980)
;    swisc = pi_n * costheta = pi_n * kz/kr
;    twisc = swisc - pi_nm1
;    tau_n = n * twisc - pi_nm1
   e = gpuDivFAT(npts, kz, pi_n.handle, 1., kr.handle, 0., swisc.handle)
   err OR= cudaThreadSynchronize()
   e = gpuSubF(npts, swisc.handle, pi_nm1.handle, twisc.handle)
   err OR= cudaThreadSynchronize()
   e = gpuSubFAT(npts, 1., pi_nm1.handle, n, twisc.handle, 0., tau_n.handle)
   ; NOTE: tau_n appears as -tau_n: redefine accordingly.

; ... Riccati-Bessel function, page 478: \xi_n(kr)
;    xi_n   = (2.d*n - 1.d) * xi_nm1 / kr - xi_nm2
   e = gpuDivF(npts, xi_nm1_r.handle, kr.handle, xi_n_r.handle)
   err OR= cudaThreadSynchronize()
   e = gpuAddFAT(npts, 2.*n-1., xi_n_r.handle, -1., xi_nm2_r.handle, 0., xi_n_r.handle)
   e = gpuDivF(npts, xi_nm1_i.handle, kr.handle, xi_n_i.handle)
   err OR= cudaThreadSynchronize()
   e = gpuAddFAT(npts, 2.*n-1., xi_n_i.handle, -1., xi_nm2_i.handle, 0., xi_n_i.handle)
   err OR= cudaThreadSynchronize()

; vector spherical harmonics (4.50)
;    Mo1n[0,*] = 0.d             ; no radial component (as initialized)

;    Mo1n[1,*] = pi_n * xi_n     ; ... divided by cosphi/kr
   e = gpuMultF(npts, pi_n.handle, xi_n_r.handle, Mo1n2_r.handle)
   e = gpuMultF(npts, pi_n.handle, xi_n_i.handle, Mo1n2_i.handle)

;    Mo1n[2,*] = -tau_n * xi_n   ; ... divided by sinphi/kr
   e = gpuMultF(npts, tau_n.handle, xi_n_r.handle, Mo1n3_r.handle)
   e = gpuMultF(npts, tau_n.handle, xi_n_i.handle, Mo1n3_i.handle)
   
;    dn = (n * xi_n)/kr - xi_nm1
   e = gpuDivF(npts, xi_n_r.handle, kr.handle, dn_r.handle)
   e = gpuDivF(npts, xi_n_i.handle, kr.handle, dn_i.handle)
   err OR= cudaThreadSynchronize()
   e = gpuSubFAT(npts, n, dn_r.handle, 1., xi_nm1_r.handle, 0., dn_r.handle)
   e = gpuSubFAT(npts, n, dn_i.handle, 1., xi_nm1_i.handle, 0., dn_i.handle)
   err OR= cudaThreadSynchronize()

;    Ne1n[0,*] = n*(n + 1.d) * pi_n * xi_n ; ... divided by cosphi sintheta/kr^2
   e = gpuMultFAT(npts, n*(n+1.), pi_n.handle, 1., xi_n_r.handle, 0, Ne1n1_r.handle)
   e = gpuMultFAT(npts, n*(n+1.), pi_n.handle, 1., xi_n_i.handle, 0, Ne1n1_i.handle)

;    Ne1n[1,*] = -tau_n * dn     ; ... divided by cosphi/kr
   e = gpuMultF(npts, tau_n.handle, dn_r.handle, Ne1n2_r.handle)
   e = gpuMultF(npts, tau_n.handle, dn_i.handle, Ne1n2_i.handle)

;    Ne1n[2,*] = pi_n * dn       ; ... divided by sinphi/kr
   e = gpuMultF(npts, pi_n.handle, dn_r.handle, Ne1n3_r.handle)
   e = gpuMultF(npts, pi_n.handle, dn_i.handle, Ne1n3_i.handle)

   err OR= cudaThreadSynchronize()

;;;;;;
;
; Upward recurrences ...
;
; ... angular functions (4.47)
; Method described by Wiscombe (1980)
;    pi_nm1 = pi_n
;    pi_n = swisc + (n + 1.d) * twisc / n
   temp = pi_nm1
   pi_nm1 = pi_n
   pi_n = temp
   e = gpuAddFAT(npts, 1., swisc.handle, (n+1.)/n, twisc.handle, 0., pi_n.handle)
; ... Riccati-Bessel function
;    xi_nm2 = xi_nm1
;    xi_nm1 = xi_n
   temp = xi_nm2_r
   xi_nm2_r = xi_nm1_r
   xi_nm1_r = xi_n_r
   xi_n_r = temp
   temp = xi_nm2_i
   xi_nm2_i = xi_nm1_i
   xi_nm1_i = xi_n_i
   xi_n_i = temp

;;;;;;
;
; Scattered field in spherical coordinates (4.45)
;    Es += En * (ci * ab[0,n] * Ne1n - ab[1,n] * Mo1n)

; prefactor, page 93
   En = ci^n * (2.*n + 1.)/ n / (n + 1.)
   an = En * ci * ab[0,n]
   bn = En * ab[1,n]
   an_r = real_part(an)
   an_i = imaginary(an)
   bn_r = real_part(bn)
   bn_i = imaginary(bn)

;Es1 ; Note that Mo1n1 = 0
   e = gpuAddFAT(npts, an_r, Ne1n1_r.handle, -an_i, Ne1n1_i.handle, 0., x.handle)
   e = gpuAddFAT(npts, an_i, Ne1n1_r.handle,  an_r, Ne1n1_i.handle, 0., y.handle)
   err OR= cudaThreadSynchronize()
   e = gpuAddF(npts, Es1_r.handle, x.handle, Es1_r.handle)
   e = gpuAddF(npts, Es1_i.handle, y.handle, Es1_i.handle)
   err OR= cudaThreadSynchronize()

;Es2
   e = gpuAddFAT(npts,  an_r, Ne1n2_r.handle, -an_i, Ne1n2_i.handle, 0., x.handle)
   e = gpuAddFAT(npts,  an_i, Ne1n2_r.handle,  an_r, Ne1n2_i.handle, 0., y.handle)
   err OR= cudaThreadSynchronize()
   e = gpuAddF(npts, Es2_r.handle, x.handle, Es2_r.handle)
   e = gpuAddF(npts, Es2_i.handle, y.handle, Es2_i.handle)
   e = gpuAddFAT(npts, -bn_r, Mo1n2_r.handle,  bn_i, Mo1n2_i.handle, 0., x.handle)
   e = gpuAddFAT(npts, -bn_i, Mo1n2_r.handle, -bn_r, Mo1n2_i.handle, 0., y.handle)
   err OR= cudaThreadSynchronize()
   e = gpuAddF(npts, Es2_r.handle, x.handle, Es2_r.handle)
   e = gpuAddF(npts, Es2_i.handle, y.handle, Es2_i.handle)
   err OR= cudaThreadSynchronize()

;Es3
   e = gpuAddFAT(npts,  an_r, Ne1n3_r.handle, -an_i, Ne1n3_i.handle, 0., x.handle)
   e = gpuAddFAT(npts,  an_i, Ne1n3_r.handle,  an_r, Ne1n3_i.handle, 0., y.handle)
   err OR= cudaThreadSynchronize()
   e = gpuAddF(npts, Es3_r.handle, x.handle, Es3_r.handle)
   e = gpuAddF(npts, Es3_i.handle, y.handle, Es3_i.handle)
   err OR= cudaThreadSynchronize()
   e = gpuAddFAT(npts, -bn_r, Mo1n3_r.handle,  bn_i, Mo1n3_i.handle, 0., x.handle)
   e = gpuAddFAT(npts, -bn_i, Mo1n3_r.handle, -bn_r, Mo1n3_i.handle, 0., y.handle)
   err OR= cudaThreadSynchronize()
   e = gpuAddF(npts, Es3_r.handle, x.handle, Es3_r.handle)
   e = gpuAddF(npts, Es3_i.handle, y.handle, Es3_i.handle)
   err OR= cudaThreadSynchronize() 
endfor

; geometric factors were divided out of the vector
; spherical harmonics for accuracy and efficiency ...
; ... put them back at the end.
;Es[0,*] *= cosphi * sintheta / kr^2
;Es[1,*] *= cosphi / kr
;Es[2,*] *= sinphi / kr

e = gpuDivF(npts, cosphi.handle, kr.handle, x.handle)
e = gpuDivF(npts, sinphi.handle, kr.handle, y.handle)
err OR= cudaThreadSynchronize()
e = gpuMultF(npts, Es1_r.handle, x.handle, Es1_r.handle)
e = gpuMultF(npts, Es1_i.handle, x.handle, Es1_i.handle)
e = gpuMultF(npts, Es2_r.handle, x.handle, Es2_r.handle)
e = gpuMultF(npts, Es2_i.handle, x.handle, Es2_i.handle)
err OR= cudaThreadSynchronize()

e = gpuMultF(npts, Es1_r.handle, y.handle, Es1_r.handle)
e = gpuMultF(npts, Es1_i.handle, y.handle, Es1_i.handle)
e = gpuMultF(npts, Es3_r.handle, y.handle, Es3_r.handle)
e = gpuMultF(npts, Es3_i.handle, y.handle, Es3_i.handle)
err OR= cudaThreadSynchronize()

;;;;;;
;
; At this point, the scattered wave is expressed in spherical
; coordinates.  Because the reference beam is expected to be
; linearly polarized along the x axis, the interference pattern
; depends on the x component of the scattered field in the focal
; plane.  Project the solution onto the x axis.
Ex_r = xi_nm2_r ; reuse GPU variables
Ex_i = xi_nm2_i                      

;    Ec[0,*] =  Es[0,*] * sintheta * cosphi
;    Ec[0,*] += Es[1,*] * costheta * cosphi
;    Ec[0,*] -= Es[2,*] * sinphi
e = gpuMultF(npts, cosphi.handle, sintheta.handle, x.handle)
e = gpuDivFAT(npts, kz, cosphi.handle, 1., kr.handle, 0., y.handle)
err OR= cudaThreadSynchronize()

e = gpuMultF(npts, Es1_r.handle, x.handle, Ex_r.handle)
e = gpuMultF(npts, Es1_i.handle, x.handle, Ex_i.handle)

e = gpuMultF(npts, Es2_r.handle, y.handle, rho.handle)
e = gpuMultF(npts, Es2_i.handle, y.handle, kr.handle)
err OR= cudaThreadSynchronize()
e = gpuAddF(npts, Ex_r.handle, rho.handle, Ex_r.handle)
e = gpuAddF(npts, Ex_i.handle, kr.handle,  Ex_i.handle)
err OR= cudaThreadSynchronize()
   
e = gpuMultF(npts, Es3_r.handle, sinphi.handle, x.handle)
e = gpuMultF(npts, Es3_i.handle, sinphi.handle, y.handle)
err OR= cudaThreadSynchronize()
e = gpuSubF(npts, Ex_r.handle, x.handle, Ex_r.handle)
e = gpuSubF(npts, Ex_i.handle, y.handle, Ex_i.handle)

;;;;;;
;
; Calculate holographic image
;
dhm = pi_n             ; dhm = Re{ \vec{E}_s(\vec{r}) \cdot \uvec x \exp(ikz) }
e = gpuAddFAT(npts, cos(kz), Ex_r.handle, sin(kz), Ex_i.handle, 0., dhm.handle)

;;;;;;
;
; Include background and self-image
;
if n_elements(alpha) eq 1 then begin
   fsq = pi_nm1                 ; fsq = |E_s|^2

   e = gpuMultF(npts, Es1_r.handle, Es1_r.handle, Es1_r.handle)
   e = gpuMultF(npts, Es1_i.handle, Es1_i.handle, Es1_i.handle)
   err OR= cudaThreadSynchronize()
   e = gpuAddF(npts, Es1_r.handle, Es1_i.handle, fsq.handle) ; |E_r|^2

   e = gpuMultF(npts, Es2_r.handle, Es2_r.handle, Es2_r.handle)
   e = gpuMultF(npts, Es2_i.handle, Es2_i.handle, Es2_i.handle)
   err OR= cudaThreadSynchronize()
   e = gpuAddF(npts, Es2_r.handle, Es2_i.handle, x.handle) ; |E_\theta|^2

   err OR= cudaThreadSynchronize()
   e = gpuAddF(npts, fsq.handle, x.handle, fsq.handle)
   
   e = gpuMultF(npts, Es3_r.handle, Es3_r.handle, Es3_r.handle)
   e = gpuMultF(npts, Es3_i.handle, Es3_i.handle, Es3_i.handle)
   err OR= cudaThreadSynchronize()
   e = gpuAddF(npts, Es3_r.handle, Es3_i.handle, x.handle) ; |E_\phi|^2

   err OR= cudaThreadSynchronize()
   e = gpuAddF(npts, fsq.handle, x.handle, fsq.handle) ; |E_s|^2

   err OR= cudaThreadSynchronize()
   e = gpuAddFAT(npts, alpha^2, fsq.handle, 2.*alpha, dhm.handle, 1., dhm.handle)
endif

;;;;;;
;
; The Answer
;
err OR= cudaThreadSynchronize()
res = gpuGetArr(dhm)

;;;;;;
;
; Free GPU resources
;
gpuFree, [x, y, rho],                                ERROR=e & err OR= e
gpuFree, [kr, sintheta, cosphi, sinphi],             ERROR=e & err OR= e
gpuFree, [xi_n_r, xi_n_i],                           ERROR=e & err OR= e
gpuFree, [xi_nm1_r, xi_nm1_i],                       ERROR=e & err OR= e
gpuFree, [xi_nm2_r, xi_nm2_i],                       ERROR=e & err OR= e
gpuFree, [pi_n, pi_nm1],                             ERROR=e & err OR= e
gpuFree, [Mo1n2_r, Mo1n2_i, Mo1n3_r, Mo1n3_i],       ERROR=e & err OR= e
gpuFree, [dn_r, dn_i],                               ERROR=e & err OR= e
gpuFree, [Ne1n1_r, Ne1n1_i, Ne1n2_r, Ne1n2_i],       ERROR=e & err OR= e
gpuFree, [Ne1n3_r, Ne1n3_i],                         ERROR=e & err OR= e
gpuFree, [Es1_r, Es1_i, Es2_r, Es2_i, Es3_r, Es3_i], ERROR=e & err OR= e

if err ne 0 then $
   message, "gpulib error:" + strtrim(err), /inf

return, reform(res, nx, ny)
end
