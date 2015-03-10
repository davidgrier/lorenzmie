;+
; NAME:
;    LorenzMie
;
; PURPOSE:
;    This object uses Lorenz-Mie theory to compute the in-line
;    hologram of a sphere with specified radius and refractive index.
;    The hologram is calculated at specified three-dimensional
;    coordinates under the assumption that the incident illumination
;    is a plane wave linearly polarized along x.
;
; CATEGORY:
;    Holographic video microscopy
;
; INHERITS:
;    generalizedLorenzMie
;
; PROPERTIES
;    NOTE: R = Required for initialization
;          I = Optionally used for initialization
;          G = Get
;          S = Set
;
;    [IGS] AP:     Radius of sphere [um]
;    [IGS] NP:     Refractive index of sphere
;    [ GS] KP:     Imaginary part of refractive index
;    [IGS] PRECISION: precision at which to keep LM coefficients.
;                  Smaller values (better precision) retain more
;                  terms at the cost of additional computational cost.
;                  Setting to 0 keeps full set of coefficients.
;                  Default: 1e-7
;
; INHERITED PROPERTIES
;    [R  ] XC:     Coordinates at which to compute hologram [pixel]
;    [R  ] YC:
;    [R  ] ZC:
;    [RGS] LAMBDA: Vacuum wavelength of light [um]
;    [RGS] MPP:    Magnification [um/pixel]
;    [RGS] NM:     Complex refractive index of medium
;
;    [ G ] AB:     Generalized Lorenz-Mie scattering coefficients    
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

;;;;;
;
; LorenzMie::SetProperty
;
pro LorenzMie::SetProperty, ap = ap,               $
                            np = np,               $
                            kp = kp,               $
                            precision = precision, $
                            _ref_extra = re
  COMPILE_OPT IDL2, HIDDEN

  updatecoefficients = 0B
  if isa(ap, /number, /scalar) then begin
     self.ap = double(ap)
     updatecoefficients = 1B
  endif

  if isa(np, /number, /scalar) then begin
     if isa(np, 'complex') || isa(np, 'dcomplex') then $
        self.np = dcomplex(np) $
     else $
        self.np = dcomplex(np, imaginary(self.np))
     updatecoefficients = 1B
  endif

  if isa(kp, /number, /scalar) then begin
     self.np = dcomplex(real_part(self.np), kp)
     updatecoefficients = 1B
  endif

  if isa(precision, /number, /scalar) then begin
     self.precision = double(precision) > 0.
     updatecoefficients = 1B
  endif

  if updatecoefficients then begin
     ab = sphere_coefficients(self.ap, self.np, self.nm, self.lambda, $
                              resolution = precision)
     self.ab = ptr_new(ab, /no_copy)
  endif

  self.generalizedLorenzMie::SetProperty, _extra = re
end

;;;;;
;
; LorenzMie::GetProperties
;
pro LorenzMie::GetProperty, ap = ap,               $
                            np = np,               $
                            kp = kp,               $
                            precision = precision, $
                            _ref_extra = re
  COMPILE_OPT IDL2, HIDDEN

  if arg_present(ap) then ap = self.ap
  if arg_present(np) then np = self.np
  if arg_present(kp) then kp = imaginary(self.np)
  if arg_present(precision) then precision = self.precision

  self.generalizedLorenzMie::GetProperty, _extra = re
end

;;;;;
;
; LorenzMie::Init()
;
function LorenzMie::Init, ap = ap,               $
                          np = np,               $
                          precision = precision, $
                          _ref_extra = re
  COMPILE_OPT IDL2, HIDDEN

  if ~self.generalizedLorenzMie::Init(ab = ab, _extra = re) then $
     return, 0B

  self.ap = isa(ap, /number, /scalar) ? double(ap) : 1.D
  self.np = isa(np, /number, /scalar) ? dcomplex(np) : dcomplex(1.4)
  self.precision = isa(precision, /number, /scalar) ? double(precision) > 0d : 1d-7
  ab = sphere_coefficients(self.ap, self.np, self.nm, self.lambda, $
                           resolution = self.precision)
  self.ab = ptr_new(ab, /no_copy)

  return, 1B
end

;;;;;
;
; LorenzMie::Cleanup
;
; Handled by generalizedLorenzMie.
; Nothing to do here.
;

;;;;;
;
; LorenzMie__define
;
; Define the object structure for a LorenzMie object
;
pro LorenzMie__define

  COMPILE_OPT IDL2, HIDDEN

  struct = {LorenzMie, $
            inherits generalizedLorenzMie, $
            ap: 0.D,                       $ ; sphere radius [um]
            np: dcomplex(0),               $ ; sphere refractive index
            precision: 0.D                 $ ; precision limit for LM coefficients
           }
end
