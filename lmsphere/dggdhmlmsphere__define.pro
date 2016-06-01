;+
; NAME:
;    DGGdhmLMSphere
;
; PURPOSE:
;    This object uses Lorenz-Mie theory to compute the 
;    hologram of a sphere as would be recorded with 
;    in-line digital video microscopy.
;
; CATEGORY:
;    Holographic video microscopy, object graphics
;
; INHERITS:
;    IDL_Object
;
; PROPERTIES:
;    NOTE: R = Required for initialization
;          I = Optionally used for initialization
;          G = Get
;          S = Set
;
;    [RG ] DIM:    [nx, ny] dimensions of hologram [pixels].
;    [IG ] R0:     coordinates of lower-left pixel [pixels]
;                  Default: [0,0]
;    [RGS] LAMBDA: vacuum wavelength of light [um]
;    [RGS] MPP:    magnification [um/pixel]

;    
;    [ GS] RP: [xp, yp, zp] position of the center of the sphere 
;                  relative to the center of the image in the focal
;                  plane. [pixel]
;    [ GS] XP: x-coordinate of the sphere's center [pixel]
;    [ GS] YP: y-coordinate of the sphere's center [pixel]
;    [ GS] ZP: z-coordinate of the sphere's center [pixel]
;    [ GS] AP: radius of sphere [um]
;    [ GS] NP: complex refractive index of sphere
;    [ GS] NM: complex refractive index of medium
;    [ GS] ALPHA: relative amplitude of illumination
;    [ GS] DELTA: wavefront distortion [pixel]
; 
;    [ G ] HOLOGRAM: real-valued computed holographic image
;    [ GS] AB:       Generalized Lorenz-Mie scattering coefficients
;    [ GS] RESOLUTION: Resolution of Lorenz-Mie coefficients.
;                      0: Full resolution
;                      1e-5: No noticable loss of quality.
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
; 5. F. C. Cheong, B. Sun, R. Dreyfus, J. Amato-Grill, K. Xiao,
;    L. Dixon and D. G. Grier,
;    "Flow visualization and flow cytometry with holographic video
;    microscopy," Opt. Express 17, 13071-13079 (2009).
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
;    Optionally limit resolution of Lorenz-Mie coefficients.
; 07/24/2013 DGG Anonymous data structures allow for resizing.
; 07/25/2013 DGG CPU code returns correct fields.
; 02/21/2015 DGG use V pointer property as hook for subclasses
;    to provide alternative geommetries.
; 03/08/2015 DGG remove GPU code.  The idea is to implement
;    all acceleration modes by subclassing the base class.
; 06/01/2016 DGG allow ap and np to be arrays for stratified spheres.
;
; NOTES:
; Integrate sphere_coefficient code?
; Allow for indexing points -- calculate only at specified points
; Use masking to handle deinterlacing.
; 
; Copyright (c) 2011-2016 David G. Grier
;-    

;;;;
;
; DGGdhmLMSphere::UpdateGeometry
;
; Update geometric factors on CPU
;
pro DGGdhmLMSphere::UpdateGeometry

  COMPILE_OPT IDL2, HIDDEN

  v = self.geometry
  coordinates = self.coordinates
  
  (*v).x = (*coordinates).x - self.rp[0]
  (*v).y = (*coordinates).y - self.rp[1]
  
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

  (*v).kr *= self.k             ; reduced radial coordinate
  (*v).sinkr = sin((*v).kr)
  (*v).coskr = cos((*v).kr)
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
;currentexcept = !except
;!except = 0
;void = check_math()

  ci = dcomplex(0, 1)           ; imaginary unit

  ab = *self.ab                 ; Lorenz-Mie coefficients
  nc = n_elements(ab[0,*]) - 1  ; number of terms
  
  v = self.geometry             ; preallocated variables

  ;; starting points for recursive function evaluation ...
  ;; ... Riccati-Bessel radial functions, page 478
  xi_nm2 = dcomplex((*v).coskr, (*v).sinkr) ; \xi_{-1}(kr)
  xi_nm1 = dcomplex((*v).sinkr,-(*v).coskr) ; \xi_0(kr)

  ;; ... angular functions (4.47), page 95
  pi_nm1 = 0.d                  ; \pi_0(\cos\theta)
  pi_n   = 1.d                  ; \pi_1(\cos\theta)

  Mo1n = dcomplexarr((*v).npts, 3, /NOZERO)
  Mo1n[*, 0] = dcomplex(0)
  Ne1n = dcomplexarr((*v).npts, 3, /NOZERO)
  Es = dcomplexarr((*v).npts, 3)

  ;; Compute field by summing multipole contributions
  for n = 1.d, nc do begin

  ;; upward recurrences ...
  ;; ... Legendre factor (4.47)
  ;; Method described by Wiscombe (1980)
     swisc = pi_n * (*v).costheta
     twisc = swisc - pi_nm1
     tau_n = pi_nm1 - n * twisc ; -\tau_n(\cos\theta)

  ;; ... Riccati-Bessel function, page 478
     xi_n = (2.d*n - 1.d) * (xi_nm1 / (*v).kr) - xi_nm2 ; \xi_n(kr)

  ;; ... Deirmendjian's derivative
     dn = (n * xi_n) / (*v).kr - xi_nm1

  ;; vector spherical harmonics (4.50)
  ;;   Mo1n[*, 0] = dcomplex(0.d)   ; no radial component
     Mo1n[0, 1] = pi_n * xi_n   ; ... divided by cosphi/kr
     Mo1n[0, 2] = tau_n * xi_n  ; ... divided by sinphi/kr

     Ne1n[0, 0] = n*(n + 1.d) * pi_n * xi_n ; ... divided by cosphi sintheta/kr^2
     Ne1n[0, 1] = tau_n * dn                ; ... divided by cosphi/kr
     Ne1n[0, 2] = pi_n  * dn                ; ... divided by sinphi/kr

  ;; prefactor, page 93
     En = ci^n * (2.d*n + 1.d) / n / (n + 1.d)

  ;; the scattered field in spherical coordinates (4.45)
     Es += (En * ci * ab[0, n]) * Ne1n
     Es -= (En * ab[1, n]) * Mo1n

  ;; upward recurrences ...
  ;; ... angular functions (4.47)
  ;; Method described by Wiscombe (1980)
     pi_nm1 = pi_n
     pi_n = swisc + ((n + 1.d) / n) * twisc

  ;; ... Riccati-Bessel function
     xi_nm2 = xi_nm1
     xi_nm1 = xi_n
  endfor

  ;; geometric factors were divided out of the vector
  ;; spherical harmonics for accuracy and efficiency ...
  ;; ... put them back at the end.
  Es[*, 0] *= (*v).cosphi * (*v).sintheta / (*v).kr^2
  Es[*, 1] *= (*v).cosphi / (*v).kr
  Es[*, 2] *= (*v).sinphi / (*v).kr

  ;; angular demagnification
  ;; Abbe sine condition
  ;; sinthetap = noil * sintheta/M
  ;; costhetap = sqrt(1 - sinthetap^2)
  
  ;; solid angle projection
  ;; fac = sqrt(costhetap/costheta)
  ;; Es[*, 1] *= fac
  ;; Es[*, 2] *= fac

  ;; k-vector rotation
  ;; costheta = costhetap
  ;; sintheta = sinthetap
  
  ;;; Hologram
  ;;;
  ;;; I(\vec{r}) = |\hat{x} + \alpha \exp(-i k zp) \vec{E}_s(\vec{r})|^2
  ;;;
  ;; NOTE: Project \hat{x} onto spherical coordinates
  ;; \hat{x} = \sin\theta \cos\phi \hat{r} 
  ;;           + \cos\theta \cos\phi \hat{\theta}
  ;;           - \sin\phi \hat{\phi}
  ;;
  Ne1n = (self.alpha * exp(dcomplex(0, -self.k*(self.rp[2] + self.delta)))) * Es
  Ne1n[*, 0] += (*v).cosphi * (*v).sintheta
  Ne1n[*, 1] += (*v).cosphi * (*v).costheta
  Ne1n[*, 2] -= (*v).sinphi

  *self.hologram = reform(total(real_part(Ne1n * conj(Ne1n)), 2), self.nx, self.ny)

;status = check_math()
;!except = currentexcept                ; restore math error checking
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
                                 ap = ap_, $
                                 np = np_, $
                                 kp = kp, $
                                 nm = nm, $
                                 km = km, $
                                 resolution =  resolution, $
                                 alpha = alpha, $
                                 delta = delta, $
                                 lambda = lambda, $
                                 mpp = mpp, $
                                 deinterlace = deinterlace, $
                                 dim = dim

  COMPILE_OPT IDL2, HIDDEN

  if arg_present(deinterlace) then begin
     message, "The DEINTERLACE flag can only be set at initialization", /inf
     return
  endif

  if arg_present(dim) then begin
     message, "The DIMENSION keyword can only be set at initialization", /inf
     return
  endif

  newcoeffs = 0B              ; flag to determine if new Lorenz-Mie coefficients are required

  ;;; Geometry
  if (doupdate = isa(xp, /scalar, /number)) then $
     self.rp[0] = double(xp)
  if (doupdate or= isa(yp, /scalar, /number)) then $
     self.rp[1] = double(yp)
  if (doupdate or= isa(zp, /scalar, /number)) then $
     self.rp[2] = double(zp)
  if (doupdate or= (isa(rp, /number, /array) && (n_elements(rp) eq 3))) then $
     self.rp = double(rp)
  if doupdate then $
     self.UpdateGeometry

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

  if isa(ap_, /number) then begin
     ap = double(ap_)
     self.ap = ptr_new(ap, /no_copy)
     newcoeffs = 1B
  endif

  if isa(np_, /number) then begin
     np = dcomplex(np_)
     self.np = ptr_new(np, /no_copy)
     newcoeffs = 1B
  endif

  if isa(kp, /scalar, /number) then begin
     np = dcomplex(real_part(*self.np), kp)
     self.np = ptr_new(np, /no_copy)
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

  if isa(resolution, /scalar, /number) then begin
     self.resolution = double(resolution)
     newcoeffs = 1B
  endif

  if isa(lambda, /scalar, /number) then begin
     self.lambda = double(lambda)
     newcoeffs = 1B
  endif

  self.k = 2.d * !dpi * real_part(self.nm) * self.mpp / self.lambda

  ;;; compute new Lorenz-Mie coefficients, if necessary
  if (newcoeffs) then begin
     ab = sphere_coefficients(*self.ap, *self.np, self.nm, self.lambda, resolution = self.resolution)
     self.ab = ptr_new(ab, /no_copy)
  endif

  ;;; compute hologram
  self.Compute
end

;;;;
;
; DGGdhmLMSphere::GetProperty
;
; Retrieve the value of properties associated with the hologram
;
pro DGGdhmLMSphere::GetProperty, hologram = hologram, $
                                 dim = dim, $
                                 ab = ab, $
                                 resolution = resolution, $
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
                                 geometry = geometry

  COMPILE_OPT IDL2, HIDDEN

  if arg_present(hologram) then $
     hologram = *(self.hologram)

  if arg_present(ab) then ab = *(self.ab)
  if arg_present(resolution) then resolution = self.resolution

  if arg_present(dim) then dim = self.dim
  if arg_present(rp) then rp = self.rp
  if arg_present(xp) then xp = self.rp[0]
  if arg_present(yp) then yp = self.rp[1]
  if arg_present(zp) then zp = self.rp[2]
  if arg_present(ap) then ap = *self.ap
  if arg_present(np) then np = *self.np
  if arg_present(kp) then kp = imaginary(self.np)
  if arg_present(nm) then nm = self.nm
  if arg_present(km) then km = imaginary(self.nm)
  if arg_present(alpha) then alpha = self.alpha
  if arg_present(delta) then delta = self.delta
  if arg_present(lambda) then lambda = self.lambda
  if arg_present(mpp) then mpp = self.mpp
  if arg_present(deinterlace) then deinterlace = self.deinterlace
  if arg_present(geometry) then geometry = self.geometry
end

;;;;
;
; DGGdhmLMSphere::CreateGeometry
;
pro DGGdhmLMSphere::CreateGeometry

  COMPILE_OPT IDL2,  HIDDEN

  nx = self.nx
  ny = self.ny
  npts = nx * ny
  stride = (self.deinterlace ne 0) ? 2.d : 1.d
  
  x = rebin(dindgen(nx), nx, ny, /sample)
  x = reform(x, npts, /overwrite) - self.r0[0]
  y = rebin(stride*dindgen(1, ny) + (self.deinterlace mod 2), nx, ny, /sample)
  y = reform(y, npts, /overwrite) - self.r0[1]
  
  coordinates = {x: x, $
                 y: y  $
                }
  
  geometry = {x: dblarr(npts),        $
              y: dblarr(npts),        $
              rho: dblarr(npts),      $
              kr: dblarr(npts),       $
              costheta: dblarr(npts), $
              sintheta: dblarr(npts), $
              cosphi: dblarr(npts),   $
              sinphi: dblarr(npts),   $
              coskr: dblarr(npts),    $
              sinkr: dblarr(npts),    $
              npts: npts              $
             }
  
  self.coordinates = ptr_new(coordinates, /no_copy)
  self.geometry = ptr_new(geometry, /no_copy) ; work with locally defined geometry
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
                               r0     = r0,     $ ; coordinates of lower-left pixel
                               xp = xp,         $ ; sphere position [pixel]
                               yp = yp,         $
                               zp = zp,         $
                               rp = rp,         $ ; 3D position [pixel]
                               ap = ap_,        $ ; sphere radius [um]
                               np = np_,        $ ; sphere refractive index
                               nm = nm,         $ ; medium refractive index
                               resolution = resolution, $ ; resolution for Lorenz-Mie coefficients
                               alpha = alpha,   $ ; relative illumination amplitude
                               delta = delta,   $ ; illumination wavefront distortion
                               deinterlace = deinterlace

  COMPILE_OPT IDL2, HIDDEN

  ;;; Required inputs
  if n_elements(dim) ne 2 then begin
     message, 'dimensions must be specified with the DIM keyword', /inf
     return, 0B
  endif
  self.dim = long(dim)

  if ~isa(lambda, /scalar, /number) then begin
     message, 'wavelength must be specified with the LAMBDA keyword', /inf
     return, 0B
  endif
  self.lambda = double(lambda)

  if ~isa(mpp, /scalar, /number) then begin
     message, 'magnification must be specified with the MPP keyword', /inf
     return, 0B
  endif
  self.mpp = double(mpp)

  self.nx = self.dim[0]
  self.ny = self.dim[1]
  if isa(deinterlace, /scalar, /number) then begin
     if (deinterlace gt 0) then begin
        self.deinterlace = 2 - (deinterlace mod 2) ; 1: odd, 2: even
        self.ny = floor(self.ny/2.) + (self.ny mod 2) * (self.deinterlace - 1)
     endif
  endif

  self.r0 = isa(r0, /number, /array) ? double(r0) : [0.d, 0]

  self.CreateGeometry
  
  ;;; Optional inputs
  self.rp = (n_elements(rp) eq 3) ? double(rp) : [dim/2, 100.]

  if isa(xp, /scalar, /number) then $
     self.rp[0] = double(xp)

  if isa(yp, /scalar, /number) then $
     self.rp[1] = double(yp)

  if isa(zp, /scalar, /number) then $
     self.rp[2] = double(zp)

  ap = isa(ap_, /number) ? double(ap_) : 1.d
  self.ap = ptr_new(ap, /no_copy)

  np = isa(np_, /number) ? dcomplex(np_) : dcomplex(1.4)
  if n_elements(np) ne n_elements(*self.ap) then begin
     message, 'ap and np should have same number of elements', /inf
     return, 0B
  endif
  self.np = ptr_new(np, /no_copy)

  self.nm = isa(nm, /scalar, /number) ? dcomplex(nm) : dcomplex(1.3)

  self.alpha = isa(alpha, /scalar, /number) ? double(alpha) : 1.d

  self.delta = isa(delta, /scalar, /number) ? double(delta) : 0.d

  self.k = 2.d * !dpi * real_part(self.nm) * self.mpp / self.lambda

  self.resolution = isa(resolution, /scalar, /number) ? double(resolution) : 0

  ab = sphere_coefficients(*self.ap, *self.np, self.nm, self.lambda, $
                           resolution = self.resolution)
  self.ab = ptr_new(ab)

  ;;; allocate storage for the result
  a = dblarr(self.nx, self.ny, /nozero)
  self.hologram = ptr_new(a, /no_copy)

  self.UpdateGeometry
  self.Compute

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

  if ptr_valid(self.coordinates) then $
     ptr_free, self.coordinates
  
  if ptr_valid(self.geometry) then $
     ptr_free, self.geometry

  if ptr_valid(self.ap) then $
     ptr_free, self.ap

  if ptr_valid(self.np) then $
     ptr_free, self.np
end

;;;;
;
; DGGdhmLMSphere__define
;
; Define the object structure for a DGGdhmLMSphere object
;
pro DGGdhmLMSphere__define

  COMPILE_OPT IDL2, HIDDEN

  struct = {DGGdhmLMSphere,            $
            INHERITS     IDL_OBJECT,   $
            hologram:    ptr_new(),    $ ; computed hologram
            dim:         [0L, 0L],     $ ; dimensions of hologram
            r0:          [0.d, 0],     $ ; coordinates of lower-left corner
            ab:          ptr_new(),    $ ; Lorenz-Mie coefficients
            resolution:  0.D,          $ ; resolution to
            coordinates: ptr_new(),    $ ; Cartesian coordinates of pixels
            geometry:    ptr_new(),    $ ; local copy of preallocated variables
            rp:          dblarr(3),    $ ; 3D position [pixel]
            ap:          ptr_new(),    $ ; sphere radius [um]
            np:          ptr_new(), $ ; sphere refractive index
            nm:          dcomplex(0.), $ ; medium refractive index
            alpha:       0.D,          $ ; relative illumination amplitude
            delta:       0.D,          $ ; illumination wavefront distortion
            mpp:         0.D,          $ ; magnification [um/pixel]
            lambda:      0.D,          $ ; vacuum wavelength [um]
            k:           0.D,          $ ; wavenumber in medium [pixel^-1]
            deinterlace: 0,            $ ; 0: none, 1: odd, 2: even
            nx:          0L,           $ ; width of deinterlaced hologram
            ny:          0L            $ ; height of deinterlaced hologram
         }
end
