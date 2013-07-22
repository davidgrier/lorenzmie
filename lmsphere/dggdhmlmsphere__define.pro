;+
; NAME:
;    DGGdhmLMSphere
;
; PURPOSE:
;    This object uses Lorenz-Mie theory to compute the 
;    hologram of a sphere as would be recorded with 
;    in-line digital video microscopy.
;    The Lorenz-Mie calculation is accelerated using GPULib.
;
; CATEGORY:
;    Holographic video microscopy, object graphics
;
; INHERITS:
;    IDL_Object
;
; PROPERTIES:
;    NOTE: R = Required for initialization
;          G = Get
;          S = Set
;
;    DIM:      (RG ) [nx, ny] dimensions of hologram [pixels].
;    LAMBDA:   (RGS) vacuum wavelength of light [um]
;    MPP:      (RGS) magnification [um/pixel]
;    
;    RP:       ( GS) [XP, YP, ZP] position of the center of the sphere 
;                    relative to the center of the image in the focal
;                    plane. [pixel]
;    XP:       ( GS) x-coordinate of the sphere's center [pixel]
;    YP:       ( GS) y-coordinate of the sphere's center [pixel]
;    ZP:       ( GS) z-coordinate of the sphere's center [pixel]
;    AP:       ( GS) radius of sphere [um]
;    NP:       ( GS) complex refractive index of sphere
;    NM:       ( GS) complex refractive index of medium
;    ALPHA:    ( GS) relative amplitude of illumination
;    DELTA:    ( GS) wavefront distortion [pixel]
; 
;    HOLOGRAM: ( G ) real-valued computed holographic image
;    FIELD:    ( G ) complex-valued scattered field
;    AB:       ( GS) Lorenz-Mie scattering coefficients
;
; KEYWORDS:
;    DEINTERLACE: (Initialization)
;        Compute either the even (DEINTERLACE = 2) or
;        odd (DEINTERLACE = 1) field of an interlaced hologram.
;        Default: Return the full hologram.
;
; METHODS:
;    GetProperty: Get accessible properties
;    SetProperty: Set accessible properties
;        NOTE: GetProperty and SetProperty can be called implicitly
;        for individual properties using the IDL_Object model.
;
;    Compute: Compute hologram using present parameters.
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
; Written by David G. Grier, New York University 09/05/2011.
; 09/10/2011 DGG Added DEINTERLACE keyword.
; 10/02/2011 DGG Fixed divide-by-zero error when center of
;    hologram is aligned with the computational grid.
; 04/26/2012 DGG Fixed center-on-grid code for cases when the
;    center is out of the field of view.  Thanks to Hagay
;    Shpaisman for pointing out the bug.
; 05/03/2012 DGG Updated parameter checking.
; 06/08/2012 DGG & Paige Hasebe (NYUAD) project into spherical
;    coordinates rather than Cartesion for efficiency.
; 06/13/2012 DGG Scattered field is stored internally in spherical
;    coordinates as v.Es1, v.Es2 and v.Es3.  Added FIELD property
;    for access.
; 06/18/2012 DGG Added GetProperty method for AB: scattering
;    coefficients may be provided instead of being computed.
; 10/13/2012 DGG Added DELTA property.  Renamed to DGGdhmLMSphere
; 03/13/2013 DGG origin of coordinate system is moved to lower-left
;    corner, rather than center.
; 05/03/2013 DGG Update for GPULib 1.6.0, which fixes bugs in earlier
;    releases.
; 07/19/2013 DGG Initial version of ComputeCPU.  Refactored
;    deinterlace code.  Corrected divide-by-zero corner case for GPU.
;    Automatically select CPU or GPU calculation.  Reorder CPU arrays
;    for greater efficiency: [npts,3] rather than [3,npts].
; 07/20/2013 DGG Separate CPU geometry calculation.  Preallocate
;    CPU variables associated with geometry.
; 07/21/2013 DGG Reorganizing GPU code for reduced memory footprint,
;    and hopefully for speed after refactoring.
; 07/22/2013 DGG Single precision is not enough for recurrences.
;    Default to CPU if GPU cannot handle double precision.
;
; NOTES:
; Use low-level GPU routines for speed.
; CPU code returns incorrect field (almost never used).
; Permit ap and np to be arrays for core-shell particles
; Integrate sphere_coefficient code?
; Optionally force single precision?
; Allow for indexing points -- calculate only at specified points
; 
; Copyright (c) 2011-2013 David G. Grier
;-    

;;;;
;
; DGGdhmLMSphere::ComputeGeometryCPU
;
; Update geometric factors on CPU
;
pro DGGdhmLMSphere::ComputeGeometryCPU

COMPILE_OPT IDL2, HIDDEN

xc = self.rp[0]
yc = self.rp[1] - (self.deinterlace mod 2)

nx = self.nx
ny = self.ny
npts = nx * ny
stride = (self.deinterlace ne 0) ? 2.d : 1.d

v = self.v

(*v).x = rebin(dindgen(nx) - xc, nx, ny, /sample)
(*v).x = reform((*v).x, npts, /overwrite)
(*v).y = rebin(stride*dindgen(1, ny) - yc, nx, ny, /sample)
(*v).y = reform((*v).y, npts, /overwrite)
z = self.rp[2] > 1.d

; convert to spherical coordinates centered on the sphere.
; (r, theta, phi) is the spherical coordinate of the pixel
; at (x,y) in the imaging plane at distance z from the
; center of the sphere.
(*v).rho   = sqrt((*v).x^2 + (*v).y^2)
(*v).kr    = sqrt((*v).rho^2 + z^2)
(*v).costheta = z/(*v).kr
(*v).sintheta = (*v).rho/(*v).kr
phi   = atan((*v).y, (*v).x)
(*v).cosphi = cos(phi)
(*v).sinphi = sin(phi)

(*v).kr *= self.k                       ; reduced radial coordinate
(*v).sinkr = sin((*v).kr)
(*v).coskr = cos((*v).kr)

end

;;;;
;
; DGGdhmLMSphere::ComputeCPU
;
; Compute hologram on CPU
;
pro DGGdhmLMSphere::ComputeCPU

COMPILE_OPT IDL2, HIDDEN

; turn off math exceptions for floating point underflow errors
;currentexcept = !except
;!except = 0
;void = check_math()

ci = dcomplex(0, 1)             ; imaginary unit

ab = *self.ab                   ; Lorenz-Mie coefficients
nc = n_elements(ab[0,*]) - 1    ; number of terms

self->computeGeometryCPU

npts = self.nx * self.ny        ; number of points in hologram
;v = *(self.v)                   ; preallocated variables
v = self.v

; starting points for recursive function evaluation ...
; ... Riccati-Bessel radial functions, page 478
xi_nm2 = dcomplex((*v).coskr, (*v).sinkr) ; \xi_{-1}(kr)
xi_nm1 = dcomplex((*v).sinkr,-(*v).coskr) ; \xi_0(kr)

; ... angular functions (4.47), page 95
pi_nm1 = 0.d                    ; \pi_0(\cos\theta)
pi_n   = 1.d                    ; \pi_1(\cos\theta)

Mo1n = dcomplexarr(npts, 3, /NOZERO)
Mo1n[*, 0] = dcomplex(0)
Ne1n = dcomplexarr(npts, 3, /NOZERO)
Es = dcomplexarr(npts, 3)

; Compute field by summing multipole contributions
for n = 1.d, nc do begin

; upward recurrences ...
; ... Legendre factor (4.47)
; Method described by Wiscombe (1980)
    swisc = pi_n * (*v).costheta
    twisc = swisc - pi_nm1
    tau_n = pi_nm1 - n * twisc  ; -\tau_n(\cos\theta)

; ... Riccati-Bessel function, page 478
    xi_n = (2.d*n - 1.d) * (xi_nm1 / (*v).kr) - xi_nm2    ; \xi_n(kr)

; vector spherical harmonics (4.50)
;   Mo1n[0, 0] = 0.d             ; no radial component
    Mo1n[0, 1] = pi_n * xi_n     ; ... divided by cosphi/kr
    Mo1n[0, 2] = tau_n * xi_n    ; ... divided by sinphi/kr

    dn = (n * xi_n) / (*v).kr - xi_nm1
    Ne1n[0, 0] = n*(n + 1.d) * pi_n * xi_n ; ... divided by cosphi sintheta/kr^2
    Ne1n[0, 1] = tau_n * dn      ; ... divided by cosphi/kr
    Ne1n[0, 2] = pi_n  * dn      ; ... divided by sinphi/kr

; prefactor, page 93
    En = ci^n * (2.d*n + 1.d) / n / (n + 1.d)

; the scattered field in spherical coordinates (4.45)
    Es += (En * ci * ab[0,n]) * Ne1n - (En * ab[1,n]) * Mo1n

; upward recurrences ...
; ... angular functions (4.47)
; Method described by Wiscombe (1980)
    pi_nm1 = pi_n
    pi_n = swisc + ((n + 1.d) / n) * twisc

; ... Riccati-Bessel function
    xi_nm2 = xi_nm1
    xi_nm1 = xi_n
endfor

; geometric factors were divided out of the vector
; spherical harmonics for accuracy and efficiency ...
; ... put them back at the end.
Es[*, 0] *= (*v).cosphi * (*v).sintheta / (*v).kr^2
Es[*, 1] *= (*v).cosphi / (*v).kr
Es[*, 2] *= (*v).sinphi / (*v).kr

;;; Hologram
;;;
;;; I(\vec{r}) = |\hat{x} + \alpha \exp(-i k zp) \vec{E}_s(\vec{r})|^2
;;;
; NOTE: Project \hat{x} onto spherical coordinates
; \hat{x} = \sin\theta \cos\phi \hat{r} 
;           + \cos\theta \cos\phi \hat{\theta}
;           - \sin\phi \hat{\phi}
;
Es *= self.alpha * exp(dcomplex(0, -self.k*(self.rp[2] + self.delta)))
Es[*, 0] += (*v).cosphi * (*v).sintheta
Es[*, 1] += (*v).cosphi * (*v).costheta
Es[*, 2] -= (*v).sinphi

*self.hologram = reform(total(real_part(Es * conj(Es)), 2), self.nx, self.ny)

;status = check_math()
;!except = currentexcept                ; restore math error checking
end

;;;;
;
; DGGdhmLMSphere::Compute
;
; Compute hologram
;
pro DGGdhmLMSphere::Compute

COMPILE_OPT IDL2, HIDDEN

if ~self.gpu then begin
   self->ComputeCPU
   return
endif

; turn off math exceptions for floating point underflow errors
currentexcept = !except
!except = 0
void = check_math()

ci = dcomplex(0, 1)             ; imaginary unit

ab = *self.ab                   ; Lorenz-Mie coefficients
nc = n_elements(ab[0,*]) - 1    ; number of terms

v = *(self.v)                   ; structure of GPU variables

xc = self.rp[0]
yc = self.rp[1] - (self.deinterlace mod 2)

; handle deinterlacing
nx = self.nx
ny = self.ny
stride = (self.deinterlace ne 0) ? 2.d : 1.d

v.x = gpumake_array(nx, ny, type = self.type, /index)
v.y = gpucopy(v.x, LHS = v.y)
v.y = gpufloor(1.d, 1.d/nx, v.y, 0.d, 0.d, LHS = v.y)
v.x = gpuadd(1.d, v.x, -nx, v.y, -xc, LHS = v.x)
v.y = gpuadd(stride, v.y, 0.d, v.y, -yc, LHS = v.y)
v.z = gpuadd(0.d, v.x, 0.d, v.x, self.rp[2], LHS = v.z, /NONBLOCKING)

; rho = sqrt(x^2 + y^2)
v.rho = gpumult(v.x, v.x, LHS = v.rho, /NONBLOCKING)
v.kr  = gpumult(v.y, v.y, LHS = v.kr)
v.rho = gpuadd(v.rho, v.kr, LHS = v.rho)

; r = sqrt(rho^2 + z^2)
v.kr = gpusqrt(1.d, 1.d, v.rho, self.rp[2]^2, 0.d, LHS = v.kr)
;;; NOTE: Under GPULib 1.4.4, the 5-parameter form of gpusqrt
; incorrectly yields NAN for some input values on some GPUs, whereas
; the 1-parameter form works properly.  This seems to be fixed under
; GPULib 1.6.0.
; v.kr = gpuadd(1.d, v.rho, 0.d, v.rho, self.rp[2]^2, LHS = v.kr)
; v.kr = gpusqrt(v.kr, LHS = v.kr, /NONBLOCKING)
v.rho = gpusqrt(v.rho, LHS = v.rho)

; polar angle
v.sintheta = gpudiv(v.rho, v.kr, LHS = v.sintheta, /NONBLOCKING)
v.costheta = gpudiv(v.z,   v.kr, LHS = v.costheta, /NONBLOCKING)

; azimuthal angle
;;; correct divide-by-zero errors when rho is small
if ((abs(xc) + abs(yc)) mod 1.) lt 1e-3 then begin ; center falls on pixel center
   ;; don't do anything if the center is outside the array
   if (xc ge 0 and xc lt nx and yc ge 0 and yc lt ny) then begin
      rho = gpugetarr(v.rho)
      rho[xc, yc] = 1.
      gpuputarr, rho, v.rho
;      (v.rho)[xc, yc] = 1.
      x = gpugetarr(v.x)
      x[xc, yc] = 1.
      gpuputarr, x, v.x
;      (v.x)[xc, yc] = 1.
   endif
endif
v.sinphi = gpudiv(v.y, v.rho, LHS = v.sinphi, /NONBLOCKING)
v.cosphi = gpudiv(v.x, v.rho, LHS = v.cosphi, /NONBLOCKING)

; scale distances by wavenumber
v.kr = gpuadd(self.k, v.kr, 0.d, v.kr, 0.d, LHS = v.kr)

; starting points for recursive function evaluation ...
; ... Riccati-Bessel radial functions, page 478
; \xi_{0}(kr)  = sin(kr) - i cos(kr)
; \xi_{-1}(kr) = cos(kr) + i sin(kr)
v.x = gpucos(v.kr, LHS = v.x, /NONBLOCKING)
v.y = gpusin(v.kr, LHS = v.y)
v.xi[1] = gpuadd(1.d, v.y, -ci, v.x, 0.d, LHS = v.xi[1], /NONBLOCKING)
v.xi[2] = gpuadd(1.d, v.x,  ci, v.y, 0.d, LHS = v.xi[2], /NONBLOCKING)

; ... angular functions (4.47), page 95
; \pi_0(\cos\theta) = 0
; \pi_1(\cos\theta) = 1
v.pi[1] = gpuadd(0.d, v.z, 0.d, v.z, 0.d, LHS = v.pi[1], /NONBLOCKING)
v.pi[0] = gpuadd(0.d, v.z, 0.d, v.z, 1.d, LHS = v.pi[0], /NONBLOCKING)

; ... scattered field
v.Es1 = gpuadd(0.d, v.z, 0.d, v.z, 0.d, LHS = v.Es1, /NONBLOCKING)
v.Es2 = gpuadd(0.d, v.z, 0.d, v.z, 0.d, LHS = v.Es2, /NONBLOCKING)
v.Es3 = gpuadd(0.d, v.z, 0.d, v.z, 0.d, LHS = v.Es3, /NONBLOCKING)

; Compute field by summing multipole contributions
for n = 1.d, nc do begin

; Field calculation prefactor, page 93
   En = ci^n * (2.d * n + 1.d) / n / (n + 1.d)
   an = En * ci * ab[0, n]
   bn = -En * ab[1, n]

; upward recurrences ...

; ... Legendre factor (4.47)
; Method described by Wiscombe (1980)
;  swisc = \pi_n * \cos\theta
;  twisc = swisc - \pi_{n-1}
;  tau = n * twisc - \pi_{n-1}
   v.x = gpumult(v.pi[0], v.costheta, LHS = v.x, /NONBLOCKING) ; swisc in x

; ... Riccati-Bessel function, page 478
; \xi_n = (2n - 1) xi_{n-1} / kr - xi_{n-2}
   v.xi[0] = gpudiv(v.xi[1], v.kr, LHS = v.xi[0])
   v.xi[0] = gpuadd(2.d*n - 1.d, v.xi[0], -1.d, v.xi[2], 0.d, $
                      LHS = v.xi[0], /NONBLOCKING)

   v.y = gpusub(v.x, v.pi[1], LHS = v.y)                             ; twisc in y
   v.z = gpuadd(1.d, v.pi[1], -n, v.y, 0.d, LHS = v.z, /NONBLOCKING) ; -tau in z

; ... Deirmendjian's derivative
; d_n = (n xi_n)/kr - xi_{n-1}
   v.dn = gpudiv(v.xi[0], v.kr, LHS = v.dn)
   v.dn = gpuadd(n, v.dn, -1.d, v.xi[1], 0.d, LHS = v.dn)

; Vector spherical harmonics (4.50)
;  Mo1n[0,*] = 0.d            ; no radial component
;  Mo1n[1,*] = pi_n * xi_n    ; ... divided by cosphi/kr
;  Mo1n[2,*] = -tau_n * xi_n  ; ... divided by sinphi/kr
   v.Mo1n2 = gpumult(v.pi[0], v.xi[0], LHS = v.Mo1n2, /NONBLOCKING)
   v.Mo1n3 = gpumult(v.z, v.xi[0], LHS = v.Mo1n3, /NONBLOCKING)

;  Ne1n[0,*] = n*(n+1) * pi_n * xi_n ; ... divided by cosphi sintheta/kr^2
;  Ne1n[1,*] = -tau_n * dn    ; ... divided by cosphi/kr
;  Ne1n[2,*] = pi_n * dn      ; ... divided by sinphi/kr
   v.Ne1n1 = gpumult(n*(n+1.d), v.pi[0], 1.d, v.xi[0], 0.d, $
                       LHS = v.Ne1n1, /NONBLOCKING)
   v.Ne1n2 = gpumult(v.z, v.dn, LHS = v.Ne1n2, /NONBLOCKING)
   v.Ne1n3 = gpumult(v.pi[0], v.dn, LHS = v.Ne1n3)
   
; scattered field in spherical coordinates (4.45)
   v.Es1 = gpuadd(1.d, v.Es1, an, v.Ne1n1, 0.d, $
                    LHS = v.Es1, /NONBLOCKING)
   v.Es2 = gpuadd(1.d, v.Es2, an, v.Ne1n2, 0.d, $
                    LHS = v.Es2, /NONBLOCKING)
   v.Es3 = gpuadd(1.d, v.Es3, an, v.Ne1n3, 0.d, LHS = v.Es3)
   v.Es2 = gpuadd(1.d, v.Es2, bn, v.Mo1n2, 0.d, $
                    LHS = v.Es2, /NONBLOCKING)
   v.Es3 = gpuadd(1.d, v.Es3, bn, v.Mo1n3, 0.d, $
                    LHS = v.Es3, /NONBLOCKING)

; upward recurrences ...
; ... Riccati-Bessel function
;  xi_{n-2} = xi_{n-1}
;  xi_{n-1} = xi_n
   v.xi = shift(v.xi, 1)

; ... angular functions (4.47)
; Method described by Wiscombe (1980)
;  pi_{n-1} = pi_n
;  pi_n = swisc + (n + 1) twisc / n
   v.pi = shift(v.pi, 1)
   v.pi[0] = gpuadd(1.d, v.x, (n + 1.d)/n, v.y, 0.d, LHS = v.pi[0], /NONBLOCKING)

endfor

;;; Scattered field 
;;;
;;; \vec{E}_s(\vec{r}) in spherical coordinates
;;;
; Geometric factors were divided out of the vector spherical harmonics
; for accuracy and efficiency.  Put them back in now:
v.x = gpudiv(v.cosphi, v.kr, LHS = v.x, /NONBLOCKING)
v.y = gpudiv(v.sintheta, v.kr, LHS = v.y, /NONBLOCKING)
v.z = gpudiv(v.sinphi, v.kr, LHS = v.z)
v.Es1 = gpumult(v.Es1, v.x, LHS = v.Es1)
v.Es1 = gpumult(v.Es1, v.y, LHS = v.Es1, /NONBLOCKING)
v.Es2 = gpumult(v.Es2, v.x, LHS = v.Es2, /NONBLOCKING)
v.Es3 = gpumult(v.Es3, v.z, LHS = v.Es3)

;;; Hologram
;;;
;;; I(\vec{r}) = |\hat{x} + \alpha \exp(-i k zp) \vec{E}_s(\vec{r})|^2
;;;
; NOTE: Project \hat{x} onto spherical coordinates
; \hat{x} = \sin\theta \cos\phi \hat{r} 
;           + \cos\theta \cos\phi \hat{\theta}
;           - \sin\phi \hat{\phi}
;
fac = self.alpha * exp(dcomplex(0, -self.k*(self.rp[2] + self.delta)))
v.x = gpumult(v.sintheta, v.cosphi, LHS = v.x, /NONBLOCKING)
v.y = gpumult(v.costheta, v.cosphi, LHS = v.y, /NONBLOCKING)
v.xi[2] = gpuadd(-1., v.sinphi, fac, v.Es3, 0., LHS = v.xi[2])
v.xi[1] = gpuadd(1., v.y, fac, v.Es2, 0., LHS = v.xi[1], /NONBLOCKING)
v.xi[0] = gpuadd(1., v.x, fac, v.Es1, 0., LHS = v.xi[0])
;
; Intensity is squared magnitude
;
v.x = gpuconj(v.xi[0], LHS = v.x, /NONBLOCKING)
v.y = gpuconj(v.xi[1], LHS = v.y, /NONBLOCKING)
v.z = gpuconj(v.xi[2], LHS = v.z)
v.xi[0] = gpumult(v.xi[0], v.x, LHS = v.xi[0], /NONBLOCKING)
v.xi[1] = gpumult(v.xi[1], v.y, LHS = v.xi[1], /NONBLOCKING)
v.xi[2] = gpumult(v.xi[2], v.z, LHS = v.xi[2])
v.dn = gpuadd(v.xi[0], v.xi[1], LHS = v.dn)
v.dn = gpuadd(v.dn, v.xi[2], LHS = v.dn)

v.hologram = gpureal(v.dn, LHS = v.hologram)

*self.hologram = gpugetarr(v.hologram, LHS = *self.hologram)

status = check_math()
!except = currentexcept                ; restore math error checking
end

;;;;
;
; DGGdhmLMSphere::SetProperty
;
; Set properties associated with the hologram
;
pro DGGdhmLMSphere::SetProperty, xp = xp, $
                                 yp = yp, $
                                 zp = zp, $
                                 rp = rp, $
                                 ap = ap, $
                                 np = np, $
                                 kp = kp, $
                                 nm = nm, $
                                 km = km, $
                                 alpha = alpha, $
                                 delta = delta, $
                                 lambda = lambda, $
                                 mpp = mpp, $
                                 deinterlace = deinterlace, $
                                 type = type, $
                                 dim = dim

COMPILE_OPT IDL2, HIDDEN

if arg_present(deinterlace) then begin
   message, "The DEINTERLACE flag can only be set at initialization", /inf
   return
endif

if arg_present(type) then begin
   message, "The TYPE keyword can only be set at initialization", /inf
   return
endif

if arg_present(dim) then begin
   message, "The DIMENSION keyword can only be set at initialization", /inf
   return
endif

newcoeffs = 0B ; flag to determine if new Lorenz-Mie coefficients are required

;;; Geometry
if isa(xp, /scalar, /number) then self.rp[0] = double(xp)
if isa(yp, /scalar, /number) then self.rp[1] = double(yp)
if isa(zp, /scalar, /number) then self.rp[2] = double(zp)
if (n_elements(rp) eq 3) then self.rp = double(rp)

if isa(mpp, /scalar, /number) then $
   self.mpp = double(mpp)

;;; Illumination
if isa(alpha, /scalar, /number) then $
   self.alpha = double(alpha)

if isa(delta, /scalar, /number) then $
   self.delta = double(delta)

;;; Properties of particle and medium
if isa(ab, /number, /array) then begin
   sz = size(ab)
   if (sz[0] ne 2) or (sz[1] ne 2)  then begin
      message, 'AB: scattering coefficients must be a [2,n] array', /inf
      return
   endif
   self.ab = ptr_new(dcomplex(ab))
endif

if isa(ap, /scalar, /number) then begin
   self.ap = double(ap)
   newcoeffs = 1B
endif

if isa(np, /scalar, /number) then begin
   self.np = dcomplex(np)
   newcoeffs = 1B
endif

if isa(kp, /scalar, /number) then begin
   self.np = dcomplex(real_part(self.np), kp)
   newcoeffs = 1B
endif

if isa(nm, /scalar, /number) then begin
   self.nm = dcomplex(nm)
   newcoeffs = 1B
endif

if isa(km, /scalar, /number) then begin
   self.nm = dcomplex(real_part(self.nm), km)
   newcoeffs = 1B
endif

if isa(lambda, /scalar, /number) then begin
   self.lambda = double(lambda)
   newcoeffs = 1B
endif

self.k = 2.d * !dpi * real_part(self.nm) * self.mpp / self.lambda

;;; compute new Lorenz-Mie coefficients, if necessary
if (newcoeffs) then begin
   ab = sphere_coefficients(self.ap, self.np, self.nm, self.lambda)
   self.ab = ptr_new(ab, /no_copy)
endif

;;; compute hologram
self->compute

end

;;;;
;
; DGGdhmLMSphere::GetProperty
;
; Retrieve the value of properties associated with the hologram
;
pro DGGdhmLMSphere::GetProperty, hologram = hologram, $
                                 field = field, $
                                 dim = dim, $
                                 ab = ab, $
                                 xp = xp, $
                                 yp = yp, $
                                 zp = zp, $
                                 rp = rp, $
                                 ap = ap, $
                                 np = np, $
                                 kp = kp, $
                                 nm = nm, $
                                 km = km, $
                                 alpha = alpha, $
                                 delta = delta, $
                                 lambda = lambda, $
                                 mpp = mpp, $
                                 deinterlace = deinterlace, $
                                 type = type, $
                                 gpu = gpu

COMPILE_OPT IDL2, HIDDEN

if arg_present(hologram) then $
   hologram = *(self.hologram)

if arg_present(field) then begin
   field = dcomplexarr(3, self.dim[0], self.dim[1], /nozero)
   field[0,*,*] = gpugetarr((*self.v).Es1, LHS = field[0,*,*])
   field[1,*,*] = gpugetarr((*self.v).Es2, LHS = field[1,*,*])
   field[2,*,*] = gpugetarr((*self.v).Es3, LHS = field[2,*,*])
endif

if arg_present(ab) then $
   ab = *(self.ab)

if arg_present(dim) then dim = self.dim
if arg_present(rp) then rp = self.rp
if arg_present(xp) then xp = self.rp[0]
if arg_present(yp) then yp = self.rp[1]
if arg_present(zp) then zp = self.rp[2]
if arg_present(ap) then ap = self.ap
if arg_present(np) then np = self.np
if arg_present(kp) then kp = imaginary(self.np)
if arg_present(nm) then nm = self.nm
if arg_present(km) then km = imaginary(self.nm)
if arg_present(alpha) then alpha = self.alpha
if arg_present(delta) then delta = self.delta
if arg_present(lambda) then lambda = self.lambda
if arg_present(mpp) then mpp = self.mpp
if arg_present(deinterlace) then deinterlace = self.deinterlace
if arg_present(type) then type = self.type
if arg_present(gpu) then gpu = self.gpu

end

;;;;
;
; DGGdhmLMSphere:CPUInit
;
; Allocate CPU memory
;
pro DGGdhmLMSphere::CPUInit

COMPILE_OPT IDL2,  HIDDEN

npts = self.nx * self.ny
v = {SphereDHM_CPU_Variables, $
     x: dblarr(npts), $
     y: dblarr(npts), $
     rho: dblarr(npts), $
     kr: dblarr(npts), $
     costheta: dblarr(npts), $
     sintheta: dblarr(npts), $
     cosphi: dblarr(npts), $
     sinphi: dblarr(npts), $
     coskr: dblarr(npts), $
     sinkr: dblarr(npts) $
    }

self.v = ptr_new(v, /no_copy)
end

;;;;
;
; DGGdhmLMSphere::GPUInit
;
; Initialize GPU and allocate GPU memory
;
function DGGdhmLMSphere::GPUInit

COMPILE_OPT IDL2, HIDDEN

catch, error
if (error ne 0L) then begin
   catch, /cancel
   return, 0B
endif

;;; Initialize GPU
gpuinit, /hardware, /quiet

if ~gpudoublecapable() then $
   return, 0B

self.type = 9
ftype = 4

;;; Allocate GPU memory
nx = self.nx
ny = self.ny

v = {SphereDHM_GPU_Variables, $
     x:        gpumake_array(nx, ny, type = self.type, /NOZERO),  $
     y:        gpumake_array(nx, ny, type = self.type, /NOZERO),  $
     z:        gpumake_array(nx, ny, type = self.type, /NOZERO),  $
     rho:      gpumake_array(nx, ny, type = self.type, /NOZERO),  $
     kr:       gpumake_array(nx, ny, type = self.type, /NOZERO),  $
     sintheta: gpumake_array(nx, ny, type = self.type, /NOZERO),  $
     costheta: gpumake_array(nx, ny, type = self.type, /NOZERO),  $
     sinphi:   gpumake_array(nx, ny, type = self.type, /NOZERO),  $
     cosphi:   gpumake_array(nx, ny, type = self.type, /NOZERO),  $
     xi:      [gpumake_array(nx, ny, type = self.type, /NOZERO),  $
               gpumake_array(nx, ny, type = self.type, /NOZERO),  $
               gpumake_array(nx, ny, type = self.type, /NOZERO)], $
     pi:      [gpumake_array(nx, ny, type = self.type, /NOZERO),  $
               gpumake_array(nx, ny, type = self.type, /NOZERO)], $
     Mo1n2:    gpumake_array(nx, ny, type = self.type, /NOZERO),  $
     Mo1n3:    gpumake_array(nx, ny, type = self.type, /NOZERO),  $
     Ne1n1:    gpumake_array(nx, ny, type = self.type, /NOZERO),  $
     Ne1n2:    gpumake_array(nx, ny, type = self.type, /NOZERO),  $
     Ne1n3:    gpumake_array(nx, ny, type = self.type, /NOZERO),  $
     Es1:      gpumake_array(nx, ny, type = self.type, /NOZERO),  $
     Es2:      gpumake_array(nx, ny, type = self.type, /NOZERO),  $
     Es3:      gpumake_array(nx, ny, type = self.type, /NOZERO),  $
     dn:       gpumake_array(nx, ny, type = self.type, /NOZERO),  $
     hologram: gpumake_array(nx, ny, type = ftype, /NOZERO) $
    }
self.v = ptr_new(v, /no_copy)

return, 1B
end

;;;;
;
; DGGdhmLMSphere::Init
;
; Initialize computational pipeline
;
function DGGdhmLMSphere::Init, dim    = dim,    $ ; dimensions of hologram (R)
                               lambda = lambda, $ ; wavelength [um]        (R)
                               mpp    = mpp,    $ ; mag [um/pixel]         (R)
                               xp = xp,         $ ; sphere position [pixel]
                               yp = yp,         $
                               zp = zp,         $
                               rp = rp,         $ ; 3D position [pixel]
                               ap = ap,         $ ; sphere radius [um]
                               np = np,         $ ; sphere refractive index
                               nm = nm,         $ ; medium refractive index
                               alpha = alpha,   $ ; relative illumination amplitude
                               delta = delta,   $ ; illumination wavefront distortion
                               deinterlace = deinterlace, $
                               gpu = gpu          ; use GPU, if available

COMPILE_OPT IDL2, HIDDEN

;;; Required inputs

if (n_elements(dim) eq 2) then $
   self.dim = long(dim) $
else begin
   message, 'dimensions must be specified with the DIM keyword', /inf
   return, 0B
endelse

if isa(lambda, /scalar, /number) then $
   self.lambda = double(lambda) $
else begin
   message, 'wavelength must be specified with the LAMBDA keyword', /inf
   return, 0B
endelse

if isa(mpp, /scalar, /number) then $
   self.mpp = double(mpp) $
else begin
   message, 'magnification must be specified with the MPP keyword', /inf
   return, 0B
endelse

self.type = 9                   ; default to double-precision

self.nx = self.dim[0]
self.ny = self.dim[1]
if isa(deinterlace, /scalar, /number) then begin
   if (deinterlace gt 0) then begin
      self.deinterlace = 2 - (deinterlace mod 2) ; 1: odd, 2: even
      self.ny = floor(self.ny/2.) + (self.ny mod 2) * (self.deinterlace - 1)
   endif
endif

;;; Initialize GPU
self.gpu = (keyword_set(gpu)) ? self->GPUInit() : 0
if ~self.gpu then self->CPUInit

;;; Optional inputs
if (n_elements(rp) eq 3) then $
   self.rp = double(rp)

if isa(xp, /scalar, /number) then $
   self.rp[0] = double(xp)

if isa(yp, /scalar, /number) then $
   self.rp[1] = double(yp)

if isa(zp, /scalar, /number) then $
   self.rp[2] = double(zp)

if isa(ap, /scalar, /number) then $
   self.ap = double(ap) $
else $
   self.ap = 1.d

if isa(np, /scalar, /number) then $
   self.np = dcomplex(np) $
else $
   self.np = dcomplex(1.4)

if isa(nm, /scalar, /number) then $
   self.nm = dcomplex(nm) $
else $
   self.nm = dcomplex(1.3)

if isa(alpha, /scalar, /number) then $
   self.alpha = double(alpha) $
else $
   self.alpha = 1.d

if isa(delta, /scalar, /number) then $
   self.delta = double(delta) $
else $
   self.delta = 0.d

self.k = 2.d * !dpi * real_part(self.nm) * self.mpp / self.lambda

;;; initial Lorenz-Mie coefficients based on inputs
if isa(ab, /number, /array) then begin
   sz = size(ab)
   if (sz[0] ne 2) or (sz[1] ne 2)  then begin
      message, 'AB: scattering coefficients must be a [2,n] array', /inf
      return, 0B
   endif
endif else $
   ab = sphere_coefficients(self.ap, self.np, self.nm, self.lambda)
self.ab = ptr_new(ab)

;;; allocate storage for the result
a = dblarr(self.nx, self.ny, /nozero)
self.hologram = ptr_new(a, /no_copy)

self->compute

return, 1B
end

;;;;
;
; DGGdhmLMSphere::Cleanup
;
; Free resources used by computational pipeline
;
pro DGGdhmLMSphere::Cleanup

COMPILE_OPT IDL2, HIDDEN

if ptr_valid(self.hologram) then $
   ptr_free, self.hologram

if ptr_valid(self.v) then begin
   if self.gpu then begin
      v = *(self.v)
      gpufree, [v.x, v.y, v.z]
      gpufree, [v.rho, v.kr]
      gpufree, [v.sintheta, v.costheta, v.sinphi, v.cosphi]
      gpufree, v.xi
      gpufree, v.pi
      gpufree, [v.Mo1n2, v.Mo1n3, v.Ne1n1, v.Ne1n2, v.Ne1n3]
      gpufree, [v.Es1, v.Es2, v.Es3]
      gpufree, v.dn
      gpufree, v.hologram
   endif
   ptr_free, self.v
endif

end

;;;;
;
; DGGdhmLMSphere__define
;
; Define the object structure for a DGGdhmLMSphere object
;
pro DGGdhmLMSphere__define

COMPILE_OPT IDL2

struct = {DGGdhmLMSphere,            $
          INHERITS     IDL_OBJECT,   $
          hologram:    ptr_new(),    $ ; computed hologram
          dim:         [0L, 0L],     $ ; dimensions of hologram
          ab:          ptr_new(),    $ ; Lorenz-Mie coefficients
          v:           ptr_new(),    $ ; structure of preallocated variables
          type:        0,            $ ; data type (double if possible)
          rp:          dblarr(3),    $ ; 3D position [pixel]
          ap:          0.D,          $ ; sphere radius [um]
          np:          dcomplex(0.), $ ; sphere refractive index
          nm:          dcomplex(0.), $ ; medium refractive index
          alpha:       0.D,          $ ; relative illumination amplitude
          delta:       0.D,          $ ; illumination wavefront distortion
          mpp:         0.D,          $ ; magnification [um/pixel]
          lambda:      0.D,          $ ; vacuum wavelength [um]
          k:           0.D,          $ ; wavenumber in medium [pixel^-1]
          deinterlace: 0,            $ ; 0: none, 1: odd, 2: even
          nx:          0L,           $ ; width of deinterlaced hologram
          ny:          0L,           $ ; height of deinterlaced hologram
          gpu:         0             $ ; flag: 1 if GPU enabled
         }

end
