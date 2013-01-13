;;;;;
;
; dgglmsphere::setproperty
;
pro dgglmsphere::setproperty, ap = ap, $
                           _ref_extra = re

COMPILE_OPT IDL2, HIDDEN

self -> dgglmobject::setproperty, _extra = re
self -> dgglmobject::getproperty, calculate = calculate

if isa(ap, /number, /scalar) then begin
   calculate = 1
   self.ap = double(ap)
endif

if calculate then $
   self.ab = ptr_new(sphere_coefficients(self.ap, self.np, self.nm, self.lambda))
end

;;;;;
;
; dgglmsphere::getproperty
;
pro dgglmsphere::getproperty, ap = ap, $
                              _ref_extra = re

COMPILE_OPT IDL2, HIDDEN

self -> dgglmobject::getproperty, _extra = re
if arg_present(ap) then ap = self.ap
end

;;;;;
;
; dgglmsphere::init
;
function dgglmsphere::init, ap = ap, $
                            _ref_extra = re

COMPILE_OPT IDL2, HIDDEN

if ~self -> dgglmobject::init(_extra = re) then $
   return, 0

if isa(ap, /number, /scalar) then begin
   self.ap = double(ap)
   self.ab = ptr_new(sphere_coefficients(self.ap, self.np, self.nm, self.lambda))
endif

; override graphic
npts = 10
theta = 2.*!pi/(npts - 1.) * findgen(1, npts)
self.graphic = ptr_new([sin(theta), cos(theta), replicate(0.1, 1, npts)])

self.name = 'dgglmsphere'

return, 1
end

;;;;;
;
; dgglmsphere__define
;
; Define an object that represents a spherical Lorenz-Mie scatterer
; 
pro dgglmsphere__define

COMPILE_OPT IDL2

struct = {dgglmsphere, $
          inherits dgglmobject, $
          ap: 0.d $             ; radius [micrometers]
         }

end
          
