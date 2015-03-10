;+
; NAME:
;    generalizedLorenzMie
;
; PURPOSE:
;    This object uses generalized Lorenz-Mie theory to compute the 
;    in-line hologram of a particle with specified Lorenz-Mie scattering
;    coefficients.  The hologram is calculated at specified
;    three-dimensional coordinates under the assumption that the
;    incident illumination is a plane wave linearly polarized along x.
;
; CATEGORY:
;    Holographic video microscopy
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
;    [R  ] XC:     Coordinates at which to compute hologram [pixel]
;    [R  ] YC:
;    [R  ] ZC:
;    [RGS] LAMBDA: Vacuum wavelength of light [um]
;    [RGS] MPP:    Magnification [um/pixel]
;    [RGS] NM:     Complex refractive index of medium
;
;    [IGS] AB:     Generalized Lorenz-Mie scattering coefficients    
;    [IGS] RP:     [xp, yp, zp] Position of the center of the scatterer
;                  in the coordinate system defined by (x,y,z) [pixel]
;    [ GS] XP: x-coordinate of the scatterer's center [pixel]
;    [ GS] YP: y-coordinate of the scatterer's center [pixel]
;    [ GS] ZP: z-coordinate of the scatterer's center [pixel]
;
;    [IGS] ALPHA: relative amplitude of illumination
;    [IGS] DELTA: wavefront distortion [pixel]
; 
;    [ G ] HOLOGRAM: real-valued computed holographic image
;    [ G ] FIELD:    complex-valued scattered field
;
; METHODS:
;    GetProperty: Get accessible properties
;    SetProperty: Set accessible properties
;        NOTE: GetProperty and SetProperty can be called implicitly
;        for individual properties using the IDL_Object model.
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
;
; Copyright (c) 2011-2015 David G. Grier
;-    

;;;;
;
; generalizedLorenzMie::Compute
;
; Compute hologram
;
pro generalizedLorenzMie::Compute

  COMPILE_OPT IDL2, HIDDEN

  if ~ptr_valid(self.ab) then return ; not yet initialized
  ab = *self.ab                      ; Lorenz-Mie coefficients
  nc = n_elements(ab[0,*]) - 1       ; number of terms

  ci = dcomplex(0, 1)           ; imaginary unit

  v = self.v

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

  ;;; Hologram
  ;;;
  ;;; I(\vec{r}) = |\hat{x} + \alpha \exp(-i k zp) \vec{E}_s(\vec{r})|^2
  ;;;
  ;; NOTE: Project \hat{x} onto spherical coordinates
  ;; \hat{x} = \sin\theta \cos\phi \hat{r} 
  ;;           + \cos\theta \cos\phi \hat{\theta}
  ;;           - \sin\phi \hat{\phi}
  ;;
  Es *= self.alpha * exp(-ci*self.k*(self.rp[2] + self.delta))
  Es += (*v).E0

  self.hologram = ptr_new(total(real_part(Es * conj(Es)), 2), /no_copy)
end

;;;;
;
; generalizedLorenzMie::SetProperty
;
; Set properties associated with the hologram
;
pro generalizedLorenzMie::SetProperty, ab = ab, $
                                       rp =  rp, $
                                       xp = xp, $
                                       yp = yp, $
                                       zp = zp, $
                                       lambda = lambda, $
                                       mpp = mpp, $
                                       alpha = alpha, $
                                       delta = delta, $
                                       nm = nm, $
                                       km = km

  COMPILE_OPT IDL2, HIDDEN

  ;;; Properties of scatterer
  if isa(ab, /number, /array) then begin
     sz = size(ab)
     if (sz[0] ne 2) || (sz[1] ne 2)  then begin
        message, 'AB: scattering coefficients must be a [2,n] array', /inf
        return
     endif
     self.ab = ptr_new(dcomplex(ab))
  endif
  
  ;;; Geometry
  orp = self.rp
  if isa(xp, /scalar, /number) then self.rp[0] = double(xp)
  if isa(yp, /scalar, /number) then self.rp[1] = double(yp)
  if isa(zp, /scalar, /number) then self.rp[2] = double(zp)
  if (n_elements(rp) eq 3) then self.rp = double(rp)
  if total(orp ne self.rp) gt 0 then self.UpdateGeometry

  ;;; Optics
  if isa(lambda, /scalar, /number) then self.lambda = double(lambda)
  if isa(mpp,    /scalar, /number) then self.mpp = double(mpp)
  if isa(alpha,  /scalar, /number) then self.alpha = double(alpha)
  if isa(delta,  /scalar, /number) then self.delta = double(delta)

  ;;; Medium
  if isa(nm, /scalar, /number) then self.nm = dcomplex(nm)

  if isa(km, /scalar, /number) then self.nm = dcomplex(real_part(self.nm), km)

  self.k = 2.d * !dpi * real_part(self.nm) * self.mpp / self.lambda
end

;;;;
;
; generalizedLorenzMie::GetProperty
;
; Retrieve the value of properties associated with the hologram
;
pro generalizedLorenzMie::GetProperty, hologram = hologram, $
                                       field = field, $
                                       ab = ab, $
                                       xp = xp, $
                                       yp = yp, $
                                       zp = zp, $
                                       rp = rp, $
                                       nm = nm, $
                                       km = km, $
                                       alpha = alpha, $
                                       delta = delta, $
                                       lambda = lambda, $
                                       mpp = mpp, $
                                       coordinates = coordinates, $
                                       geometry = geometry

  COMPILE_OPT IDL2, HIDDEN

  if arg_present(hologram) then begin
     self.compute
     hologram = *(self.hologram)
  endif
  if arg_present(field) then field = *((*self.v).Es)
  if arg_present(ab) then ab = *(self.ab)
  if arg_present(rp) then rp = self.rp
  if arg_present(xp) then xp = self.rp[0]
  if arg_present(yp) then yp = self.rp[1]
  if arg_present(zp) then zp = self.rp[2]
  if arg_present(nm) then nm = self.nm
  if arg_present(km) then km = imaginary(self.nm)
  if arg_present(alpha) then alpha = self.alpha
  if arg_present(delta) then delta = self.delta
  if arg_present(lambda) then lambda = self.lambda
  if arg_present(mpp) then mpp = self.mpp
  if arg_present(coordinates) then coordinates = self.coordinates
  if arg_present(geometry) then geometry = self.geometry
end

;;;;
;
; generalizedLorenzMie::UpdateGeometry
;
; Update geometric factors on CPU
;
pro generalizedLorenzMie::UpdateGeometry

  COMPILE_OPT IDL2, HIDDEN

  v = self.geometry
  coordinates = self.coordinates
  
  (*v).x = (*coordinates).x - self.rp[0]
  (*v).y = (*coordinates).y - self.rp[1]
  (*v).z = (*coordinates).z + self.rp[2]

; convert to spherical coordinates centered on the sphere.
; (r, theta, phi) is the spherical coordinate of the pixel
; at (x,y) in the imaging plane at distance z from the
; center of the sphere.
  (*v).rho = sqrt((*v).x^2 + (*v).y^2)
  (*v).kr  = sqrt((*v).rho^2 + (*v).z^2)
  (*v).costheta = (*v).z/(*v).kr
  (*v).sintheta = (*v).rho/(*v).kr
  phi = atan((*v).y, (*v).x)
  (*v).cosphi = cos(phi)
  (*v).sinphi = sin(phi)

  (*v).kr *= self.k             ; reduced radial coordinate
  (*v).sinkr = sin((*v).kr)
  (*v).coskr = cos((*v).kr)

  ;; incident field in spherical coordinates
  (*v).E0[*, 0] = (*v).cosphi * (*v).sintheta
  (*v).E0[*, 1] = (*v).cosphi * (*v).costheta
  (*v).E0[*, 2] = -(*v).sinphi

end

;;;;
;
; generalizedLorenzMie::CreateGeometry
;
function generalizedLorenzMie::CreateGeometry, x, y, z

  COMPILE_OPT IDL2,  HIDDEN

  ;; Check coordinate dimensions
  nx = n_elements(x)
  ny = n_elements(y)
  nz = n_elements(z)
  ;; Must be 1 or equal to each other
  npts = nx > ny > nz
  if ((nx ne 1) && (nx ne npts)) || $
     ((ny ne 1) && (ny ne npts)) || $
     ((nz ne 1) && (nz ne npts)) then $
        return, 0B
  
  coordinates = {x: x, $
                 y: y, $
                 z: z  $
                }

  var = dblarr(npts, /NOZERO)
  fvar = dcomplexarr(npts, 3, /NOZERO)
  ;; cvar = dcomplexarr(npts)
  geometry = {x: x, $
              y: y, $
              z: z, $
              rho:      var, $
              kr:       var, $
              costheta: var, $
              sintheta: var, $
              cosphi:   var, $
              sinphi:   var, $
              coskr:    var, $
              sinkr:    var, $
              E0:       fvar, $
  ;; xi_nm2           
  ;; xi_nm1
  ;; xi_n
  ;; Mo1n
  ;; Ne1n
  ;; swisc
  ;; twisc
  ;; tau_n
              npts: npts              $
             }
  
  self.coordinates = ptr_new(coordinates, /no_copy)
  self.geometry = ptr_new(geometry, /no_copy) ; locally defined geometry

  self.v = self.geometry

  return, 1B
end

;;;;
;
; generalizedLorenzMie::Init()
;
; Initialize computational pipeline
;
function generalizedLorenzMie::Init, xc = x, $ ; coordinates of hologram pixels    (R)
                                     yc = y, $
                                     zc = z, $
                                     lambda = lambda, $ ; wavelength [um]         (R)
                                     mpp    = mpp,    $ ; mag [um/pixel]          (R)
                                     nm     = nm,     $ ; medium refractive index (R)
                                     ab = ab,         $ ; Lorenz-Mie coefficients
                                     rp = rp,         $ ; 3D particle position [pixel]
                                     xp = xp,         $ ; particle position [pixel]
                                     yp = yp,         $
                                     zp = zp,         $
                                     alpha = alpha,   $ ; relative illumination amplitude
                                     delta = delta      ; illumination wavefront distortion

  COMPILE_OPT IDL2, HIDDEN

  ;;; Required inputs
  if ~isa(x, /number) || ~isa(y, /number) || ~isa(z, /number) then begin
     message, 'Coordinates must be specified.', /inf
     return, 0B
  endif

  if ~self.CreateGeometry(x, y, z) then begin
     message, 'Coordinates must be scalars or vectors of compatible lengths.', /inf
     return, 0B
  endif
  
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

  self.nm = isa(nm, /scalar, /number) ? dcomplex(nm) : dcomplex(1.3)

  self.k = 2.d * !dpi * real_part(self.nm) * self.mpp / self.lambda
  
  ;;; Optional inputs
  ;; Lorenz-Mie coefficients
  if isa(ab, /number, /array) then begin
     sz = size(ab)
     if (sz[0] ne 2) || (sz[1] ne 2)  then begin
        message, 'AB: scattering coefficients must be a [2,n] array', /inf
        return, 0B
     endif
     self.ab = ptr_new(dcomplex(ab))
  endif

  ;; Particle position
  self.rp = (n_elements(rp) eq 3) ? double(rp) : [mean(x), mean(y), 100.]

  if isa(xp, /scalar, /number) then $
     self.rp[0] = double(xp)

  if isa(yp, /scalar, /number) then $
     self.rp[1] = double(yp)

  if isa(zp, /scalar, /number) then $
     self.rp[2] = double(zp)

  self.UpdateGeometry
  
  ;; Optics
  self.alpha = isa(alpha, /scalar, /number) ? double(alpha) : 1.d
  self.delta = isa(delta, /scalar, /number) ? double(delta) : 0.d

  return, 1B
end

;;;;
;
; generalizedLorenzMie::Cleanup
;
; Free resources used by computational pipeline
;
pro generalizedLorenzMie::Cleanup

  COMPILE_OPT IDL2, HIDDEN

  if ptr_valid(self.coordinates) then $
     ptr_free, self.coordinates
  
  if ptr_valid(self.geometry) then $
     ptr_free, self.geometry

  if ptr_valid(self.hologram) then $
     ptr_free, self.hologram
end

;;;;
;
; generalizedLorenzMie__define
;
; Define the object structure for a generalizedLorenzMie object
;
pro generalizedLorenzMie__define

  COMPILE_OPT IDL2, HIDDEN

  struct = {generalizedLorenzMie,      $
            INHERITS     IDL_OBJECT,   $
            coordinates: ptr_new(),    $ ; Cartesian coordinates of pixels
            ab:          ptr_new(),    $ ; Lorenz-Mie coefficients
            resolution:  0.D,          $ ; resolution limit for LM coefficients
            rp:          dblarr(3),    $ ; 3D particle position [pixel]
            v:           ptr_new(),    $ ; pointer to structure of preallocated variables
            geometry:    ptr_new(),    $ ; local copy of preallocated variables
            lambda:      0.D,          $ ; vacuum wavelength [um]
            nm:          dcomplex(0.), $ ; medium refractive index
            alpha:       0.D,          $ ; relative illumination amplitude
            delta:       0.D,          $ ; illumination wavefront distortion
            mpp:         0.D,          $ ; magnification [um/pixel]
            k:           0.D,          $ ; wavenumber in medium [pixel^-1]
            hologram:    ptr_new()     $ ; computed hologram
         }
end
