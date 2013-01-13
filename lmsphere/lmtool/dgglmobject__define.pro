;;;;;
;
; dgglmobject::draw
;
; Update graphic representation of object
;
pro dgglmobject::draw

COMPILE_OPT IDL2, HIDDEN

graphic = *self.graphic
graphic[0, *] += self.rp[0]
graphic[1, *] += self.rp[1]
self -> IDLgrPolyline::SetProperty, data = graphic

end

;;;;;
;
; dgglmobject::setproperty
;
pro dgglmobject::setproperty, rp = rp, $
                              np = np, $
                              nm = nm, $
                              lambda = lambda, $
                              calculate = calculate, $
                              _ref_extra = re

COMPILE_OPT IDL2, HIDDEN

self -> IDLgrPolyline::Setproperty, _extra = re

if isa(rp, /number) then begin
   case n_elements(rp) of
      2: self.rp[0:1] = double(rp)
      3: self.rp = double(rp)
   endcase
endif

self.calculate = isa(calculate, /number, /scalar) ? calculate : 0

if isa(np, /number) then begin
   self.calculate = 1
   self.np = dcomplex(np)
endif

if isa(nm, /number, /scalar) then begin
   self.calculate = 1
   self.nm = dcomplex(nm)
endif

if isa(lambda, /number, /scalar) then begin
   self.calculate = 1
   self.lambda = double(lambda)
endif

end

;;;;;
;
; dgglmobject::getproperty
;
pro dgglmobject::getproperty, rp = rp, $
                              np = np, $
                              nm = nm, $
                              lambda = lambda, $
                              ab = ab, $
                              calculate = calculate, $
                              _ref_extra = re

COMPILE_OPT IDL2, HIDDEN

self -> IDLgrPolyline::GetProperty, _extra = re

if arg_present(rp) then rp = self.rp
if arg_present(np) then np = self.np
if arg_present(nm) then nm = self.nm
if arg_present(lambda) then lambda = self.lambda
if arg_present(ab) then ab = *(self.ab)
if arg_present(calculate) then calculate = self.calculate
end

;;;;;
;
; dgglmobject::init
;
function dgglmobject::init, rp = rp, $
                            np = np, $
                            nm = nm, $
                            lambda = lambda, $
                            graphic = graphic

COMPILE_OPT IDL2, HIDDEN

if isa(rp, /number) and n_elements(rp) eq 3 then $
   self.rp = double(rp)

self.lambda = isa(lambda, /number, /scalar) ? double(lambda) : 0.532d
self.np = isa(np, /number, /scalar) ? dcomplex(np) : 1.5d
self.nm = isa(nm, /number, /scalar) ? dcomplex(nm) : refractiveindex(self.lambda, 0.24)

if ~isa(graphic, /number, /array) then $
   self.graphic = ptr_new(dblarr(2))

return, 1
end

;;;;;
;
; dgglmobject::cleanup
;
; Free resources used by dgglmobject
;
pro dgglmobject::cleanup

COMPILE_OPT IDL2, HIDDEN

ptr_free, self.graphics
ptr_free, self.ab
end

;;;;;
;
; dgglmobject__define
;
; Define an object that represents a Lorenz-Mie scatterer
; 
pro dgglmobject__define

COMPILE_OPT IDL2

struct = {dgglmobject, $
          inherits IDL_Object, $
          inherits IDLgrPolyline, $
          graphic: ptr_new(), $ ; graphical representation
          rp: dblarr(3), $      ; coordinates of center
          np: dcomplex(0.), $   ; complex refractive index of object
          nm: dcomplex(0.), $   ; complex refractive index of medium
          lambda: 0.d, $        ; wavelength [micrometers]
          ab: ptr_new(),  $     ; Lorenz-Mie scattering coefficients
          calcuate: 0b $        ; flag
         }

end
          
