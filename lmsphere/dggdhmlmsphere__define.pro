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
;    coordinates as gpu.Es1, gpu.Es2 and gpu.Es3.  Added FIELD property
;    for access.
; 06/18/2012 DGG Added GetProperty method for AB: scattering
;    coefficients may be provided instead of being computed.
; 10/13/2012 DGG Added DELTA property.  Renamed to DGGdhmLMSphere
;
; NOTES:
; Implement CPU calculations for systems without GPULib.
; Permit ap and np to be arrays for core-shell particles
; Integrate sphere_coefficient code?
; 
; Copyright (c) 2011-2012 David G. Grier
;-    
;;;;
;
; DGGdhmLMSphere::ComputeCPU
;
; Compute hologram on CPU
;
pro DGGdhmLMSphere::ComputeCPU

COMPILE_OPT IDL2, HIDDEN

;;; hook for CPU-based computation
;;; not implemented yet

end

;;;;
;
; DGGdhmLMSphere::Compute
;
; Compute hologram
;
pro DGGdhmLMSphere::Compute

COMPILE_OPT IDL2, HIDDEN

; turn off math exceptions for floating point underflow errors
currentexcept = !except
!except = 0
void = check_math()

lambda = self.lambda / real_part(self.nm) / self.mpp ; wavelength [pixels]
k = 2.d * !dpi / lambda                              ; wavenumber [pixel^-1]
ci = dcomplex(0, 1)                                  ; imaginary unit

ab = *self.ab                   ; Lorenz-Mie coefficients
nc = n_elements(ab[0,*]) - 1    ; number of terms

gpu = *(self.gpu)               ; structure of GPU variables

nx = double(self.dim[0])
ny = double(self.dim[1])
xc = (nx - 1.d)/2.d + self.rp[0]
yc = (ny - 1.d)/2.d + self.rp[1]

; handle deinterlacing
fac = 1.d
if (self.deinterlace ne 0) then begin
   ny = floor(ny/2.) + (ny mod 2) * (self.deinterlace - 1)
   yc -= self.deinterlace mod 2
   fac = 2.d
endif

gpu.x = gpumake_array(nx, ny, type = self.type, /index)
gpu.y = gpucopy(gpu.x, LHS = gpu.y)
gpu.y = gpufloor(1.d, 1.d/nx, gpu.y, 0.d, 0.d, LHS = gpu.y)
gpu.x = gpuadd(1.d, gpu.x, -nx, gpu.y, -xc, LHS = gpu.x)
gpu.y = gpuadd(fac, gpu.y, 0.d, gpu.y, -yc, LHS = gpu.y)
gpu.z = gpuadd(0.d, gpu.x, 0.d, gpu.x, self.rp[2], LHS = gpu.z, /NONBLOCKING)

; rho = sqrt(x^2 + y^2)
gpu.rho = gpumult(gpu.x, gpu.x, LHS = gpu.rho, /NONBLOCKING)
gpu.kr  = gpumult(gpu.y, gpu.y, LHS = gpu.kr)
gpu.rho = gpuadd(gpu.rho, gpu.kr, LHS = gpu.rho)

; r = sqrt(rho^2 + z^2)
;;; NOTE: the 5-parameter form of gpusqrt incorrectly yields NAN for
; some input values on some GPUs, whereas the 1-parameter form works properly.
; gpu.kr = gpusqrt(1.d, 1.d, gpu.rho, self.rp[2]^2, 0.d, LHS = gpu.kr)
gpu.kr = gpuadd(1.d, gpu.rho, 0.d, gpu.rho, self.rp[2]^2, LHS = gpu.kr)
gpu.kr = gpusqrt(gpu.kr, LHS = gpu.kr, /NONBLOCKING)
gpu.rho = gpusqrt(gpu.rho, LHS = gpu.rho)

; polar angle
gpu.sintheta = gpudiv(gpu.rho, gpu.kr, LHS = gpu.sintheta, /NONBLOCKING)
gpu.costheta = gpudiv(gpu.z,   gpu.kr, LHS = gpu.costheta, /NONBLOCKING)

; azimuthal angle
gpu.sinphi = gpudiv(gpu.y, gpu.rho, LHS = gpu.sinphi, /NONBLOCKING)
gpu.cosphi = gpudiv(gpu.x, gpu.rho, LHS = gpu.cosphi, /NONBLOCKING)
;;; NOTE: this results in a divide-by-zero error when x = y = 0.
; The following hack fixes the problem under the constraint
; that gpu array subscripting does not yet work.
if ((abs(xc) + abs(yc)) mod 1.) lt 1e-3 then begin
   ; don't do anything if the center is outside the array
   if (xc ge 0 and xc lt nx and yc ge 0 and yc lt ny) then begin
      temp = gpugetarr(gpu.sinphi)
      temp[xc,yc] = 0.
      gpu.sinphi = gpuputarr(temp, LHS = gpu.sinphi)
      temp = gpugetarr(gpu.cosphi, LHS = temp)
      temp[xc,yc] = 1.
      gpu.cosphi = gpuputarr(temp, LHS = gpu.cosphi)
   endif
endif

; scale distances by wavenumber
gpu.kr = gpuadd(k, gpu.kr, 0.d, gpu.kr, 0.d, LHS = gpu.kr)

; starting points for recursive function evaluation ...
; ... Riccati-Bessel radial functions, page 478
; \xi_{0}(kr)  = sin(kr) - i cos(kr)
; \xi_{-1}(kr) = cos(kr) + i sin(kr)
gpu.x = gpucos(gpu.kr, LHS = gpu.x, /NONBLOCKING)
gpu.y = gpusin(gpu.kr, LHS = gpu.y)
gpu.xi[1] = gpuadd(1.d, gpu.y, -ci, gpu.x, 0.d, LHS = gpu.xi[1], /NONBLOCKING)
gpu.xi[2] = gpuadd(1.d, gpu.x,  ci, gpu.y, 0.d, LHS = gpu.xi[2], /NONBLOCKING)

; ... angular functions (4.47), page 95
; \pi_0(\cos\theta) = 0
; \pi_1(\cos\theta) = 1
gpu.pi[1] = gpuadd(0.d, gpu.z, 0.d, gpu.z, 0.d, LHS = gpu.pi[1], /NONBLOCKING)
gpu.pi[0] = gpuadd(0.d, gpu.z, 0.d, gpu.z, 1.d, LHS = gpu.pi[0], /NONBLOCKING)

; ... scattered field
gpu.Es1 = gpuadd(0.d, gpu.z, 0.d, gpu.z, 0.d, LHS = gpu.Es1, /NONBLOCKING)
gpu.Es2 = gpuadd(0.d, gpu.z, 0.d, gpu.z, 0.d, LHS = gpu.Es2, /NONBLOCKING)
gpu.Es3 = gpuadd(0.d, gpu.z, 0.d, gpu.z, 0.d, LHS = gpu.Es3, /NONBLOCKING)

; Compute field by summing multipole contributions
for n = 1.d, nc do begin
; upward recurrences ...
; ... Legendre factor (4.47)
; Method described by Wiscombe (1980)
;  swisc = \pi_n * \cos\theta
;  twisc = swisc - \pi_{n-1}
;  tau = n * twisc - \pi_{n-1}
   gpu.x = gpumult(gpu.pi[0], gpu.costheta, LHS = gpu.x) ; swisc in x
   gpu.y = gpusub(gpu.x, gpu.pi[1], LHS = gpu.y)         ; twisc in y
   gpu.z = gpuadd(1.d, gpu.pi[1], -n, gpu.y, 0.d, $
                  LHS = gpu.z, /NONBLOCKING) ; -tau in z

; ... Riccati-Bessel function, page 478
; \xi_n = (2n - 1) xi_{n-1} / kr - xi_{n-2}
   gpu.xi[0] = gpudiv(gpu.xi[1], gpu.kr, LHS = gpu.xi[0])
   gpu.xi[0] = gpuadd(2.d*n - 1.d, gpu.xi[0], -1.d, gpu.xi[2], 0.d, $
                      LHS = gpu.xi[0])

; ... Deirmendjian's derivative
; d_n = (n xi_n)/kr - xi_{n-1}
   gpu.dn = gpudiv(gpu.xi[0], gpu.kr, LHS = gpu.dn)
   gpu.dn = gpuadd(n, gpu.dn, -1.d, gpu.xi[1], 0.d, $
                   LHS = gpu.dn, /NONBLOCKING)

; Vector spherical harmonics (4.50)
;  Mo1n[0,*] = 0.d            ; no radial component
;  Mo1n[1,*] = pi_n * xi_n    ; ... divided by cosphi/kr
;  Mo1n[2,*] = -tau_n * xi_n  ; ... divided by sinphi/kr
   gpu.Mo1n2 = gpumult(gpu.pi[0], gpu.xi[0], LHS = gpu.Mo1n2, /NONBLOCKING)
   gpu.Mo1n3 = gpumult(gpu.z, gpu.xi[0], LHS = gpu.Mo1n3, /NONBLOCKING)

;  Ne1n[0,*] = n*(n+1) * pi_n * xi_n ; ... divided by cosphi sintheta/kr^2
;  Ne1n[1,*] = -tau_n * dn    ; ... divided by cosphi/kr
;  Ne1n[2,*] = pi_n * dn      ; ... divided by sinphi/kr
   gpu.Ne1n1 = gpumult(n*(n+1.d), gpu.pi[0], 1.d, gpu.xi[0], 0.d, $
                       LHS = gpu.Ne1n1, /NONBLOCKING)
   gpu.Ne1n2 = gpumult(gpu.z, gpu.dn, LHS = gpu.Ne1n2, /NONBLOCKING)
   gpu.Ne1n3 = gpumult(gpu.pi[0], gpu.dn, LHS = gpu.Ne1n3)
   
; upward recurrences ...
; ... angular functions (4.47)
; Method described by Wiscombe (1980)
;  pi_{n-1} = pi_n
;  pi_n = swisc + (n + 1) twisc / n
   gpu.pi = shift(gpu.pi, 1)
   gpu.pi[0] = gpuadd(1.d, gpu.x, (n + 1.d)/n, gpu.y, 0.d, $
                      LHS = gpu.pi[0], /NONBLOCKING)

; ... Riccati-Bessel function
;  xi_{n-2} = xi_{n-1}
;  xi_{n-1} = xi_n
   gpu.xi = shift(gpu.xi, 1)

; Field calculation
; prefactor, page 93
   En = ci^n * (2.d * n + 1.d) / n / (n + 1.d)
   an = En * ci * ab[0, n]
   bn = -En * ab[1, n]
; scattered field in spherical coordinates (4.45)
   gpu.Es1 = gpuadd(1.d, gpu.Es1, an, gpu.Ne1n1, 0.d, $
                    LHS = gpu.Es1, /NONBLOCKING)
   gpu.Es2 = gpuadd(1.d, gpu.Es2, an, gpu.Ne1n2, 0.d, $
                    LHS = gpu.Es2, /NONBLOCKING)
   gpu.Es3 = gpuadd(1.d, gpu.Es3, an, gpu.Ne1n3, 0.d, LHS = gpu.Es3)
   gpu.Es2 = gpuadd(1.d, gpu.Es2, bn, gpu.Mo1n2, 0.d, $
                    LHS = gpu.Es2, /NONBLOCKING)
   gpu.Es3 = gpuadd(1.d, gpu.Es3, bn, gpu.Mo1n3, 0.d, $
                    LHS = gpu.Es3, /NONBLOCKING)
endfor

;;; Scattered field 
;;;
;;; \vec{E}_s(\vec{r}) in spherical coordinates
;;;
; Geometric factors were divided out of the vector spherical harmonics
; for accuracy and efficiency.  Put them back in now:
gpu.x = gpudiv(gpu.cosphi, gpu.kr, LHS = gpu.x, /NONBLOCKING)
gpu.y = gpudiv(gpu.sintheta, gpu.kr, LHS = gpu.y, /NONBLOCKING)
gpu.z = gpudiv(gpu.sinphi, gpu.kr, LHS = gpu.z)
gpu.Es1 = gpumult(gpu.Es1, gpu.x, LHS = gpu.Es1)
gpu.Es1 = gpumult(gpu.Es1, gpu.y, LHS = gpu.Es1, /NONBLOCKING)
gpu.Es2 = gpumult(gpu.Es2, gpu.x, LHS = gpu.Es2, /NONBLOCKING)
gpu.Es3 = gpumult(gpu.Es3, gpu.z, LHS = gpu.Es3)

;;; Hologram
;;;
;;; I(\vec{r}) = |\hat{x} + \alpha \exp(-i k zp) \vec{E}_s(\vec{r})|^2
;;;
; NOTE: Project \hat{x} onto spherical coordinates
; \hat{x} = \sin\theta \cos\phi \hat{r} 
;           + \cos\theta \cos\phi \hat{\theta}
;           - \sin\phi \hat{\phi}
;
fac = self.alpha * exp(dcomplex(0, -k*(self.rp[2] + self.delta)))
gpu.x = gpumult(gpu.sintheta, gpu.cosphi, LHS = gpu.x, /NONBLOCKING)
gpu.y = gpumult(gpu.costheta, gpu.cosphi, LHS = gpu.y, /NONBLOCKING)
gpu.xi[2] = gpuadd(-1., gpu.sinphi, fac, gpu.Es3, 0., LHS = gpu.xi[2])
gpu.xi[1] = gpuadd(1., gpu.y, fac, gpu.Es2, 0., LHS = gpu.xi[1], /NONBLOCKING)
gpu.xi[0] = gpuadd(1., gpu.x, fac, gpu.Es1, 0., LHS = gpu.xi[0])
;
; Intensity is squared magnitude
;
gpu.x = gpuconj(gpu.xi[0], LHS = gpu.x, /NONBLOCKING)
gpu.y = gpuconj(gpu.xi[1], LHS = gpu.y, /NONBLOCKING)
gpu.z = gpuconj(gpu.xi[2], LHS = gpu.z)
gpu.xi[0] = gpumult(gpu.xi[0], gpu.x, LHS = gpu.xi[0], /NONBLOCKING)
gpu.xi[1] = gpumult(gpu.xi[1], gpu.y, LHS = gpu.xi[1], /NONBLOCKING)
gpu.xi[2] = gpumult(gpu.xi[2], gpu.z, LHS = gpu.xi[2])
gpu.dn = gpuadd(gpu.xi[0], gpu.xi[1], LHS = gpu.dn)
gpu.dn = gpuadd(gpu.dn, gpu.xi[2], LHS = gpu.dn)

gpu.hologram = gpureal(gpu.dn, LHS = gpu.hologram)

*self.hologram = gpugetarr(gpu.hologram, LHS = *self.hologram)

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
                                 type = type

COMPILE_OPT IDL2, HIDDEN

if arg_present(hologram) then $
   hologram = *(self.hologram)

if arg_present(field) then begin
   field = dcomplexarr(3, self.dim[0], self.dim[1], /nozero)
   field[0,*,*] = gpugetarr((*self.gpu).Es1, LHS = field[0,*,*])
   field[1,*,*] = gpugetarr((*self.gpu).Es2, LHS = field[1,*,*])
   field[2,*,*] = gpugetarr((*self.gpu).Es3, LHS = field[2,*,*])
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

if (self.type eq 9) then $
   self.type = (gpudoublecapable()) ? 9 : 6 ; dcomplex or complex
ftype = (self.type eq 9) ? 5 : 4            ; double or float

;;; Allocate GPU memory
nx = self.dim[0]
ny = self.dim[1]
if (self.deinterlace ne 0) then begin
   fac = (ny mod 2) * (self.deinterlace - 1)
   ny = floor(ny/2.) + fac
endif

gpu = {SphereDHM_GPU_Variables, $
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
self.gpu = ptr_new(gpu, /no_copy)

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
                               rp = rp,         $  ; 3D position [pixel]
                               ap = ap,         $  ; sphere radius [um]
                               np = np,         $  ; sphere refractive index
                               nm = nm,         $  ; medium refractive index
                               alpha = alpha,   $  ; relative illumination amplitude
                               delta = delta,   $  ; illumination wavefront distortion
                               deinterlace = deinterlace, $
                               single = single

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

self.type = (keyword_set(single)) ? 6 : 9

if isa(deinterlace, /scalar, /number) then $
   self.deinterlace = 2 - (deinterlace mod 2) ; 1: odd, 2: even

;;; Initialize GPU
if ~(self->GPUInit()) then $
   return, 0B

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

;;; allocate storage for the result
ny = dim[1]
if (self.deinterlace ne 0) then begin
   fac = (ny mod 2) * (self.deinterlace - 1)
   ny = floor(ny/2.) + fac
endif
a = dblarr(dim[0], ny, /nozero)
self.hologram = ptr_new(a, /no_copy)

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

;;; safe initial parameters
;res = (self.type eq 6) ? 1e-3 : 1e-7
;if (abs(self.rp[0]) + abs(self.rp[1])) mod 1.d lt res then $
;   self.rp[0] += res

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

if ptr_valid(self.gpu) then begin
   gpu = *(self.gpu)
   gpufree, [gpu.x, gpu.y, gpu.z]
   gpufree, [gpu.rho, gpu.kr]
   gpufree, [gpu.sintheta, gpu.costheta, gpu.sinphi, gpu.cosphi]
   gpufree, gpu.xi
   gpufree, gpu.pi
   gpufree, [gpu.Mo1n2, gpu.Mo1n3, gpu.Ne1n1, gpu.Ne1n2, gpu.Ne1n3]
   gpufree, [gpu.Es1, gpu.Es2, gpu.Es3]
   gpufree, gpu.dn
   gpufree, gpu.hologram
   ptr_free, self.gpu
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
          gpu:         ptr_new(),    $ ; structure of GPU variables
          cpu:         ptr_new(),    $ ; structure of CPU variables
          type:        0,            $ ; data type (double if possible)
          rp:          dblarr(3),    $ ; 3D position [pixel]
          ap:          0.D,          $ ; sphere radius [um]
          np:          dcomplex(0.), $ ; sphere refractive index
          nm:          dcomplex(0.), $ ; medium refractive index
          alpha:       0.D,          $ ; relative illumination amplitude
          delta:       0.D,          $ ; illumination wavefront distortion
          lambda:      0.D,          $ ; vacuum wavelength [um]
          mpp:         0.D,          $ ; magnification [um/pixel]
          deinterlace: 0             $ ; 0: none, 1: odd, 2: even
         }

end
