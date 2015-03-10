;+
; NAME:
;    hvmLorenzMie
;
; PURPOSE:
;    This object uses Lorenz-Mie theory to compute the in-line
;    hologram of a sphere with specified radius and refractive index.
;    The hologram is calculated for selected pixels taken from a
;    square array representing the microscope's focal plane, under
;    the assumption that the incident illumination is a plane wave
;    linearly polarized along x.
;
; CATEGORY:
;    Holographic video microscopy
;
; INHERITS:
;    dhmLorenzMie
;
; PROPERTIES
;    NOTE: R = Required for initialization
;          I = Optionally used for initialization
;          G = Get
;          S = Set
;
;    [I  ] FRACTION: Fraction of pixels to consider.
;                  Default: 0.5 unless the DEINTERLACE option is set.
;
;    [I  ] DEINTERLACE: Choose one field of an interlaced video frame.
;                  0 or not set: use all pixels.
;                                FRACTION defaults to 0.5.
;                  EVEN: calculate hologram for even field
;                  ODD:  calculate hologram for odd field
;                  Default: not set, use all pixels
;
;    [ G ] MASK:   Selected pixels within the array whose size
;                  is set by DIMENSIONS.
;
; INHERITED PROPERTIES
;    [RGS] DIMENSIONS: [nx, ny] dimensions of the rectangular
;                  hologram [pixels]
;    [IGS] R0:     [x0, y0] coordinates of the lower-left pixel
;                  of the hologram.  [pixels]
;                  Default: [0,0]
;
;    [RGS] LAMBDA: Vacuum wavelength of light [um]
;    [RGS] MPP:    Magnification [um/pixel]
;    [RGS] NM:     Complex refractive index of medium
;
;    [IGS] AP:     Radius of sphere [um]
;    [IGS] NP:     Refractive index of sphere
;    [ GS] KP:     Imaginary part of refractive index
;    [ G ] AB:     Generalized Lorenz-Mie scattering coefficients
;    [IGS] PRECISION: precision at which to keep LM coefficients.
;                  Smaller values (better precision) retain more
;                  terms at the cost of additional computational cost.
;                  Setting to 0 keeps full set of coefficients.
;                  Default: 1e-7
;
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
; NOTES:
;    Use supplied mask
;    Change deinterlace and fraction on the fly.
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
; hvmLorenzMie::GetProperty
;
pro hvmLorenzMie::GetProperty, hologram = hologram, $
                               deinterlace = deinterlace, $
                               fraction = fraction, $
                               mask = mask, $
                               _ref_extra = re

  COMPILE_OPT IDL2, HIDDEN

  if arg_present(hologram) then $
     self.LorenzMie::GetProperty, hologram = hologram

  if arg_present(deinterlace) then $
     deinterlace = self.deinterlace
  
  if arg_present(fraction) then $
     fraction = self.fraction

  if arg_present(mask) then $
     mask = *self.mask

  self.dhmLorenzMie::GetProperty, _extra = re
end

;;;;;
;
; hvmLorenzMie::MakeCoordinates()
;
function hvmLorenzMie::MakeCoordinates, dimensions, r0

  COMPILE_OPT IDL2, HIDDEN

  if n_params() eq 2 then begin
     grid = self.dhmLorenzMie::MakeCoordinates(dimensions, r0)
     self.grid = grid
  endif else $
     grid = self.grid

  npts = n_elements(grid['x'])
  if self.deinterlace ne 0 then $
     deinterlace = where((grid['y'] mod 2) eq (self.deinterlace mod 2), npts)

  ;;; generate random mask
  if self.fraction ge 1. then $
     mask = deinterlace $
  else begin
     mask = long(npts * randomu(seed, round(self.fraction*npts)))
     mask = mask[sort(mask)]
     u = uniq(mask)
     mask = mask[u]
     if isa(deinterlace) then $
        mask = deinterlace[mask]
  endelse

  ;;; use masked coordinates
  c = hash('x', (grid['x'])[mask], $
           'y', (grid['y'])[mask], $
           'z', grid['z'])

  ;;; save mask
  self.mask = ptr_new(mask, /no_copy)

  return, c
end

;;;;;
;
; hvmLorenzMie::Init()
;
function hvmLorenzMie::Init, deinterlace = deinterlace, $
                             fraction = fraction, $
                             _ref_extra = re

  COMPILE_OPT IDL2, HIDDEN

  self.deinterlace = 0
  if isa(deinterlace, /number, /scalar) && (deinterlace ge 1) then $
        self.deinterlace = 2 - (long(deinterlace) mod 2)
    
  self.fraction = isa(fraction, /number, /scalar) ? $
                  float(fraction) > 0. < 1. : $
                  ((self.deinterlace ge 1) ? 1. : 0.5)

  return, self.dhmLorenzMie::Init(_extra = re)
end

;;;;;
;
; hvmLorenzMie__define
;
pro hvmLorenzMie__define

  COMPILE_OPT IDL2, HIDDEN

  struct = {hvmLorenzMie, $
            inherits dhmLorenzMie, $
            deinterlace: 0, $
            fraction: 0., $
            grid: obj_new(), $
            mask: ptr_new() $
           }
end
