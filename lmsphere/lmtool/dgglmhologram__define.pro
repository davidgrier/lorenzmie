;+
; NAME:
;    DGGlmHologram
;
; PURPOSE:
;    Object for analyzing a normalized hologram with Lorenz-Mie theory.
;
; CATEGORY:
;    Digital video microscopy
;
; CALLING SEQUENCE:
;    h = DGGlmHologram(data, lambda, mpp)
;
; INPUTS:
;    data: Normalized hologram
;
;    lambda: Vacuum wavelength of illumination [micrometer]
;
;    mpp: Magnification [micrometer/pixel]
;
; PROPERTIES:
;    data        [RGS] Normalized hologram
;    lambda      [RGS] Vacuum wavelength of illumination [micrometer]
;    mpp         [RGS] Magnification [micrometer/pixel]
;    nm          [IGS] Complex refractive index of medium.
;                      Default: refractive index of water at the given
;                      wavelength and at 24 degrees Celsius.
;    noise       [IGS] Estimate for the pixel noise.
;                      Default: Computed from data with MAD algorithm.
;    deinterlace [IGS] 0, or undefined: No deinterlacing [Default]
;                      1: odd field
;                      2: even field
;
; METHODS:
;    GetProperty
;    SetProperty
;
;    DGGlmHologram::Findfeatures
;        Analyze hologram to find features of interest
;
;    DGGlmHologram::Count()
;        Returns number of features in hologram.
;
;    DGGlmHologram::Get()
;        Retrieve features of type DGGlmFeature from list of
;        discovered features.
;    KEYWORDS:
;        ALL: return array of features
;        POSITION = n: return feature number n from list of
;            features in hologram.
;
;
; NOTES:
;    Remove SMOOTH once CT_ routines have better noise performance.
;
; MODIFICATION HISTORY:
; 01/27/2013 Written by David G. Grier, New York University
; 02/04/2013 DGG Added smooth property as a temporary fix for finding
;   features in noisy images.
;
; Copyright (c) 2013 David G. Grier
;-

;;;;;
;
; DGGlmHologram::Findfeatures
;
pro DGGlmHologram::Findfeatures, smooth = smooth

COMPILE_OPT IDL2, HIDDEN

a = isa(smooth, /number, /scalar) ? smooth(*(self.data), smooth) : *(self.data)
rp = ctfeature(a, noise = self.noise, deinterlace = self.deinterlace, $
               count = nfeatures, /quiet)

if nfeatures gt 0 then begin
   rad = ct_range(a, rp, noise = self.noise, deinterlace = self.deinterlace)
   for n = 0, nfeatures-1 do $
      self.add, dgglmfeature(parent = self, rp = rp[*, n], rad = rad[n])
endif

end

;;;;;
;
; DGGlmHologram::Setproperty
;
pro DGGlmHologram::Setproperty, data = data, $
                                lambda = lambda, $
                                mpp = mpp, $
                                nm = nm, $
                                noise = noise, $
                                deinterlace = deinterlace

COMPILE_OPT IDL2, HIDDEN

if isa(data, /number, /array) then begin
   sz = size(data)
   if sz[0] eq 2 then begin     ; two-dimensional
      med = median(data)
      if med le 1.25 and med ge 0.75 then begin ; normalized
         self.dim = sz[1:2]
         ptr_free, self.data
         self.data = ptr_new(double(data))
         if isa(noise, /scalar, /number) then $
            self.noise = mad(*(self.data))
         self.findfeatures, smooth = self.smooth
      endif
   endif
endif

if isa(lambda, /scalar, /number) then $
   self.lambda = double(lambda)

if isa(mpp, /scalar, /number) then $
   self.mpp = double(mpp)

if isa(nm, /scalar, /number) then $
   self.nm = dcomplex(nm)

if isa(noise, /scalar, /number) then $
   self.noise = noise

if n_elements(deinterlace) eq 1 then $
   self.deinterlace = deinterlace

end
;;;;;
;
; DGGlmHologram::Getproperty
;
pro DGGlmHologram::GetProperty, data = data, $
                                dim = dim, $
                                lambda = lambda, $
                                mpp = mpp, $
                                nm = nm, $
                                noise = noise, $
                                deinterlace = deinterlace

COMPILE_OPT IDL2, HIDDEN

if arg_present(data) then data = *(self.data)
if arg_present(dim) then dim = self.dim
if arg_present(lambda) then lambda = self.lambda
if arg_present(mpp) then mpp = self.mpp
if arg_present(nm) then nm = self.nm
if arg_present(noise) then noise = self.noise
if arg_present(deinterlace) then deinterlace = self.deinterlace
end

;;;;;
;
; DGGlmHologram::Init
;
; Initialize a hologram object
;
function DGGlmHologram::Init, data, $
                              lambda, $
                              mpp, $
                              nm = nm, $
                              noise = noise, $
                              deinterlace = deinterlace

COMPILE_OPT IDL2, HIDDEN

if (self->IDL_Container::Init() ne 1) then $
   return, 0

umsg = 'USAGE: h = DGGlmHologram(data, lambda, mpp)'

if n_params() ne 3 then begin
   message, umsg, /inf
   return, 0
endif

if ~isa(data, /number, /array) then begin
   message, umsg, /inf
   message, 'DATA should be a normalized two-dimensional hologram', /inf
   return, 0
endif
sz = size(data)
if sz[0] ne 2 then begin
   message, umsg, /inf
   message, 'DATA should be a two-dimensional array', /inf
   return, 0
endif
med = median(data)
if med gt 1.25 or med lt 0.75 then begin
   message, umsg, /inf
   message, 'DATA should be normalized', /inf
   return, 0
endif

if ~isa(lambda, /scalar, /number) then begin
   message, umsg, /inf
   message, 'LAMBDA should be a scalar number', /inf
   return, 0
endif

if ~isa(mpp, /scalar, /number) then begin
   message, umsg, /inf
   message, 'MPP should be a scalar number', /inf
   return, 0
endif

self.dim = sz[1:2]
self.data = ptr_new(double(data))
self.lambda = double(lambda)
self.mpp = double(mpp)

if isa(noise, /scalar, /number) then $
   self.noise = double(noise) $
else $
   self.noise = mad(*(self.data))

self.nm = isa(nm, /scalar, /number) ? dcomplex(nm) : refractiveindex(self.lambda, 24.)

if keyword_set(deinterlace) then $
   self.deinterlace = deinterlace

if isa(smooth, /scalar, /number) then $
   self.smooth = smooth

self.findfeatures, smooth = self.smooth

return, 1
end

;;;;;
;
; DGGlmHologram::Cleanup
;
; Free resources claimed by hologram
;
pro DGGlmHologram::Cleanup

COMPILE_OPT IDL2, HIDDEN

self->IDL_Container::Cleanup
ptr_free, self.data

end

;;;;;
;
; DGGlmHologram__define
;
; Define the object structure for a hologram that will be analyzed
; with Lorenz-Mie theory
;
pro DGGlmHologram__define

COMPILE_OPT IDL2

struct = {DGGlmHologram,          $
          inherits IDL_Object,    $
          inherits IDL_Container, $
          data: ptr_new(),        $ ; normalized hologram
          smooth: 0,              $ ; temporary smoothing factor for CT routines
          dim: [0L, 0L],          $ ; dimensions of hologram
          lambda: 0.d,            $ ; wavelength [micrometer]
          mpp: 0.d,               $ ; magnification [micrometer/pixel]
          nm: dcomplex(0),        $ ; refractive index of medium
          noise: 0.,              $ ; estimated pixel noise
          deinterlace: 0          $
         }
end
