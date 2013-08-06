;+
; NAME:
;    DGGlmFeature
;
; PURPOSE:
;    Object describing the hologram of an object obtained with
;    holographic video microscopy, including its analysis with
;    Lorenz-Mie theory.
;
; CATEGORY:
;    Digital video microscopy
;
; CALLING SEQUENCE:
;    f = DGGlmFeature(parent)
;
; INPUTS:
;    parent: DGGlmHologram within which the feature is found.
;
; PROPERTIES:
;    data       [ G ] region of interest cropped from parent
;    rp         [IGS] [xp, yp, zp] position of feature in hologram [pixel]
;    zp         [ GS] position of feature above focal plane [pixel]
;    rad        [IGS] half-width of region of interest [pixel]
;    dim        [ G ] [nx, ny] dimensions of data [pixel].
;                     Computed from rp and rad.
;    r0         [ G ] [x0, y0, 0.] index of lower-left corner of data in
;                     parent.  Computed from rp and rad.
;    ap         [IGS] radius of sphere [micrometers]
;    np         [IGS] complex refractive index of sphere
;    alpha      [IGS] amplitude of illumination at sphere position
;    delta      [IGS] wavefront distortion of illumination [wavelengths]
;    profile    [ G ] radial intensity profile.
;                     Computed from data, rp, and rad.
;    dprofile   [ G ] standard deviation from radial intensity profile
;                     Computed from data, rp and rad.
;    lmprofile  [ G ] Lorenz-Mie computed radial profile
;    lmhologram [ G ] Lorenz-Mie computed hologram
;
; METHODS:
;    GetProperty
;    SetProperty
;
;    DGGlmFeature::Crop
;    DGGlmFeature::Fit
;    DGGlmFeature::FitProfile
;
; NOTES: Use same estimation algorithms as LMFEATURE
;
; MODIFICATION HISTORY:
; 01/27/2013 Written by David G. Grier, New York University
; 07/24/2013 DGG Update call to RS1D.  Use GPU by default.
; 08/05/2013 DGG Support all flags for fitlmsphere.
;
; Copyright (c) 2013 David G. Grier
;-
;;;;;
;
; DGGlmFeature::Fit
;
pro DGGlmFeature::Fit, fixap = fixap, $
                       fixnp = fixnp, $
                       fixkp = fixkp, $
                       fixzp = fixzp, $
                       fixnm = fixnm, $
                       fixkm = fixkm, $
                       fixalpha = fixalpha, $
                       fixdelta = fixdelta, $
                       chisq = chisq, $
                       quiet = quiet

COMPILE_OPT IDL2, HIDDEN

rc = self.rp - self.r0
p0 = [rc, self.ap, real_part(self.np), imaginary(self.np), $
      real_part(self.parent.nm), imaginary(self.parent.nm), self.alpha, self.delta]
aa = aziavg(*(self.data), center = rc, deviates = dev, $
            deinterlace = self.parent.deinterlace)
err = abs(dev)/self.parent.noise > 1.
print, max(err), min(err)

p1 = fitlmsphere(*(self.data), p0, self.parent.lambda, self.parent.mpp, $
                 errors = err, $
                 fixap = fixap, fixnp = fixnp, fixkp = fixkp, $
                 fixzp = fixzp, fixnm = fixnm, fixkm = fixkm, $
                 fixalpha = fixalpha, fixdelta = fixdelta, $
                 deinterlace = self.parent.deinterlace, $
                 chisq = chisq, $
                 /gpu, quiet = quiet)
self.rp = p1[0, 0:1] + self.r0
self.ap = p1[0, 3]
self.np = dcomplex(p1[0, 4], p1[0, 5])
self.alpha = p1[0, 8]
self.delta = p1[0, 9]

end

;;;;;
;
; DGGlmFeature::FitProfile
;
pro DGGlmFeature::FitProfile, fixap = fixap, $
                              fixnp = fixnp, $
                              fixzp = fixzp, $
                              fixalpha = fixalpha, $
                              fixdelta = fixdelta, $
                              chisq = chisq, $
                              quiet = quiet

COMPILE_OPT IDL2, HIDDEN

p0 = [self.rp[2], self.ap, real_part(self.np), imaginary(self.np), $
      real_part(self.parent.nm), imaginary(self.parent.nm), self.alpha, self.delta]
profile = aziavg(*(self.data), center = self.rp - self.r0, rad = self.rad, $
                 deinterlace = self.parent.deinterlace)
p1 = fitlmsphere1d(profile, p0, self.parent.lambda, self.parent.mpp, $
                   chisq = chisq, quiet = quiet, $
                   fixap = fixap, fixnp = fixnp, fixzp = fixzp, $
                   fixalpha = fixalpha, fixdelta = fixdelta)
self.rp[2] = p1[0, 0]
self.ap = p1[0, 1]
self.np = dcomplex(p1[0, 2], p1[0, 3])
self.alpha = p1[0, 6]
self.delta = p1[0, 7]

end

;;;;;
;
; DGGlmFeature::EstimateA
;
; Model interference pattern as Poisson's spot to estimate
; particle radius.
;
; j0n: nth zero of J0()
; rn: radius of nth zero (minimum) of interference pattern
;
; k ap rn / zp = j0n
; ap = zp j0n / (k rn)
;
pro DGGlmFeature::EstimateA

COMPILE_OPT IDL2, HIDDEN

;; Find radii of minima crossings
c = aziavg(*(self.data), center = self.rp - self.r0, rad = self.rad, $
           deinterlace = self.parent.deinterlace) - 1.
rn = float(where(c*c[1:*] lt 0)) ; zero crossings
rn = (rn + rn[1:*])/2.
rn = rn[0:*:2]
lambda = self.parent.lambda
nm = real_part(self.parent.nm)
k = 2.*!pi*nm/lambda
;; Compare radii to those of Bessel function
j0n = [2.4048, 5.5201, 8.6537]                ; zeros of J0(x)
self.ap = self.rp[2]*mean(j0n/rn)/k ; estimated radius [um]

end

;;;;;
;
; DGGlmFeature::EstimateZ
;
; Use Rayleigh-Sommerfeld back-propagation to estimate the axial
; position of the particle
;
pro DGGlmFeature::EstimateZ

COMPILE_OPT IDL2, HIDDEN

z = dindgen(50) * 4.d + 10.d ; FIXME: set range more intelligently
lambda = self.parent.lambda
mpp = self.parent.mpp
nm = real_part(self.parent.nm)
res = rs1d(*(self.data), z, self.rp[0:1] - self.r0, lambda, mpp)
m = max(abs(res), loc)
self.rp[2] = z[loc]
end

;;;;;
;
; DGGlmFeature::Crop
;
; FIXME: update deinterlace for odd/even crop
;
pro DGGlmFeature::Crop

COMPILE_OPT IDL2, HIDDEN

if ~isa(self.parent, 'DGGlmHologram') or self.rad le 0 then return

x0 = (self.rp[0] - self.rad) > 0.d
x1 = (self.rp[0] + self.rad) < (self.parent.dim)[0] - 1.d
y0 = (self.rp[1] - self.rad) > 0.d
y1 = (self.rp[1] + self.rad) < (self.parent.dim)[1] - 1.d

data = (self.parent.data)[x0:x1, y0:y1]
self.data = ptr_new(data, /no_copy)

self.r0 = [x0, y0, 0.d]             ; lower-left corner of cropped image
self.dim = [x1, y1] - self.r0 + 1.d ; dimensions of cropped image

end

;;;;;
;
; DGGlmFeature::Setproperty
;
pro DGGlmFeature::Setproperty, parent = parent, $
                               rp = rp, $
                               rad = rad, $
                               zp = zp, $
                               np = np, $
                               ap = ap, $
                               alpha = alpha, $
                               delta = delta

COMPILE_OPT IDL2, HIDDEN

if isa(parent, 'DGGlmHologram') then $
   self.parent = parent

docrop = 0
if isa(rp, /number, /array) then begin
   case n_elements(rp) of
      2: self.rp[0:1] = rp
      3: self.rp = rp
   endcase
   docrop = 1
endif
if isa(zp, /number, /scalar) then $
   self.rp[2] = zp

if isa(rad, /number) then begin
   self.rad = long(rad[0])
   docrop = 1
endif

if isa(np, /number, /scalar) then $
   self.np = dcomplex(np)

if isa(ap, /number, /scalar) then $
   self.ap = double(ap)

if isa(alpha, /number, /scalar) then $
   self.alpha = double(alpha)

if isa(delta, /number, /scalar) then $
   self.delta = double(delta)

if docrop then self.crop

end

;;;;;
;
; DGGlmFeature::Getproperty
;
pro DGGlmFeature::Getproperty, parent = parent, $
                               data = data, $
                               dim = dim, $
                               r0 = r0, $
                               rp = rp, $
                               rad = rad, $
                               ap = ap, $
                               np = np, $
                               alpha = alpha, $
                               delta = delta, $
                               profile = profile, $
                               dprofile = dprofile, $
                               lmprofile = lmprofile, $
                               lmhologram = lmhologram


COMPILE_OPT IDL2, HIDDEN

if arg_present(parent) then parent = self.parent
if arg_present(data) and ptr_valid(self.data) then data = *self.data
if arg_present(dim) then dim = self.dim
if arg_present(r0) then r0 = self.r0
if arg_present(rp) then rp = self.rp
if arg_present(rad) then rad = self.rad
if arg_present(ap) then ap = self.ap
if arg_present(np) then np = self.np
if arg_present(alpha) then alpha = self.alpha
if arg_present(delta) then delta = self.delta

if arg_present(profile) then $
   profile = aziavg(*(self.data), center = self.rp - self.r0, rad = self.rad, $
                    deinterlace = self.parent.deinterlace)

if arg_present(dprofile) then $
   dprofile = azistd(*(self.data), center = self.rp - self.r0, rad = self.rad, $
                     deinterlace = self.parent.deinterlace)

if arg_present(lmprofile) then $
   lmprofile = lmsphereprofile(findgen(self.rad), $
                               self.rp[2], self.ap, self.np, $
                               self.parent.nm, $
                               self.alpha, self.delta, $
                               self.parent.lambda, self.parent.mpp)

if arg_present(lmhologram) then $
   lmhologram = lmsphere(self.rp - self.r0, self.ap, self.np, $
                         self.parent.nm, $
                         self.parent.lambda, self.parent.mpp, $
                         self.dim, $
                         alpha = self.alpha, delta = self.delta)
end

;;;;;
;
; DGGlmFeature::Init()
;
function DGGlmFeature::Init, parent = parent, $
                             rp = rp, $
                             rad = rad, $
                             ap = ap, $
                             np = np, $
                             zp = zp, $
                             alpha = alpha, $
                             delta = delta

COMPILE_OPT IDL2, HIDDEN

if isa(parent, 'DGGlmHologram') then $
   self.parent = parent
   
if isa(rp, /number, /array) then begin
   case n_elements(rp) of
      2: self.rp[0:1] = rp
      3: self.rp = rp
   endcase
endif
if isa(zp, /number, /scalar) then $
   rp[2] = zp

self.rad = isa(rad, /number, /scalar) ? long(rad) : 100L
self.np = isa(np, /number, /scalar) ? dcomplex(np) : dcomplex(1.5)
self.ap = isa(ap, /number, /scalar) ? double(ap) : 1.d
self.alpha = isa(alpha, /number, /scalar) ? double(alpha) : 1.d
self.delta = isa(delta, /number, /scalar) ? double(delta) : 0.d

self.crop
self.estimatez
self.fitprofile, /fixdelta, /quiet

return, 1
end

;;;;;
;
; DGGlmFeature::Cleanup
;
pro DGGlmFeature::Cleanup

COMPILE_OPT IDL2

ptr_free, data

end

;;;;;
;
; DGGlmFeature__define
;
; Define the object structure for one feature in a hologram,
; as analyzed with Lorenz-Mie theory
;
pro DGGlmFeature__define

COMPILE_OPT IDL2

struct = {DGGlmFeature,        $
          inherits IDL_Object, $
          parent: obj_new(),   $ ; DGGlmHologram
          data: ptr_new(),     $ ; cropped image
          dim: [0, 0],         $ ; dimensions of cropped image
          rp: dblarr(3),       $ ; position of feature in hologram [pixel]
          rad: 0L,             $ ; range around position to crop [pixel]
          ap: 0.d,             $ ; radius of sphere [micrometers]
          np: dcomplex(0),     $ ; refractive index of sphere
          alpha: 1.,           $ ; relative amplitude of illumination
          delta: 0.,           $ ; wavefront distortion of illumination
          r0: dblarr(3)        $ ; origin of coordinate system relative to parent
         }
end
