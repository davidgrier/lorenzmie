;+
; NAME:
;    dhmLorenzMie
;
; PURPOSE:
;    This object uses Lorenz-Mie theory to compute the in-line
;    hologram of a sphere with specified radius and refractive index.
;    The hologram is calculated on a square array representing the
;    microscope's focal plane, under the assumption that the incident illumination
;    is a plane wave linearly polarized along x.
;
; CATEGORY:
;    Holographic video microscopy
;
; INHERITS:
;    LorenzMie
;
; PROPERTIES
;    NOTE: R = Required for initialization
;          I = Optionally used for initialization
;          G = Get
;          S = Set
;
; PROPERTIES
;    [RGS] DIMENSIONS: [nx, ny] dimensions of the rectangular
;                  hologram [pixels]
;    [IGS] R0:     [x0, y0] coordinates of the lower-left pixel
;                  of the hologram.  [pixels]
;                  Default: [0,0]
;
; INHERITED PROPERTIES
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
; dhmLorenzMie::GetProperty
;
pro dhmLorenzMie::GetProperty, hologram = hologram, $
                               dimensions = dimensions, $
                               _ref_extra = re
  COMPILE_OPT IDL2, HIDDEN

  if arg_present(hologram) then begin
     self.LorenzMie::GetProperty, hologram = hologram
     hologram = reform(temporary(hologram), self.dimensions)
  endif

  if arg_present(dimensions) then dimensions = self.dimensions
  
  self.LorenzMie::GetProperty, _extra = re
end

;;;;;
;
; dhmLorenzMie::MakeCoordinates
;
; Creates a square-grid coordinate system of specified
; dimensions and with specified coordinate for the lower-left
; corner.  The coordinate system is returned as a hash with
; elements 'x', 'y' and 'z'.
;
; Overload this function to select pixels from this grid.
;
function dhmLorenzMie::MakeCoordinates, dimensions, r0

  COMPILE_OPT IDL2, HIDDEN

  if n_elements(dimensions) ne 2 then begin
     message, 'dimensions must be specified.', /inf
     return, 0B
  endif
  self.dimensions = long(dimensions)

  self.r0 = (n_elements(r0) eq 2) ? double(r0) : [0.d, 0]

  x = rebin(dindgen(self.dimensions[0]), self.dimensions, /sample)
  x = reform(x, n_elements(x), /overwrite) - self.r0[0]
  
  y = rebin(dindgen(1, self.dimensions[1]), self.dimensions, /sample)
  y = reform(y, n_elements(y), /overwrite) - self.r0[1]
  
  z = 0.d                       ; focal plane

  return, hash('x', x, 'y', y, 'z', z)
end

;;;;;
;
; dhmLorenzMie::Init()
;
function dhmLorenzMie::Init, dimensions = dimensions, $
                             r0 = r0,                 $
                             _ref_extra = re
  COMPILE_OPT IDL2, HIDDEN

  v = self.MakeCoordinates(dimensions, r0)

  return, self.LorenzMie::Init(xc = v['x'], yc = v['y'], zc = v['z'], _extra = re)
end

;;;;;
;
; dhmLorenzMie::Cleanup
;
; Handled by LorenzMie.
; Nothing to do here.
;
;;;;;
;
; dhmLorenzMie__define
;
pro dhmLorenzMie__define

  COMPILE_OPT IDL2, HIDDEN

  struct = {dhmLorenzMie,        $
            inherits LorenzMie,  $
            dimensions: [0L, 0], $
            r0: [0.d, 0]         $
           }
end
