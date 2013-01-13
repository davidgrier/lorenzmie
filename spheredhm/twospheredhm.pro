function twospheredhm, dhm1, dhm2

COMPILE_OPT IDL2

umsg = 'USAGE: b = twospheredhm(a1, a2)'
if ~isa(dhm1, 'DGGdhmSphereDHM') || ~isa(dhm2, 'DGGdhmSphereDHM') then begin
   message, umsg, /inf
   message, 'a1 and a2 should be objects of type DGGdhmSphereDHM', /inf
   return, -1
endif

; two holograms must be compatible
if (~array_equal(dhm1.dim, dhm2.dim) || $ ; dimensions
    dhm1.lambda ne dhm2.lambda || $       ; wavelength
    dhm1.mpp ne dhm2.mpp || $             ; magnification
    dhm1.nm ne dhm2.nm) then begin        ; medium refractive index
   message, umsg, /inf
   message, 'a1 and a2 must have the same dim, lambda, mpp, and nm', /inf
   return, -1
endif

; holograms of individual spheres
b1 = dhm1.hologram
b2 = dhm2.hologram

; interference between two spheres' scattered fields
ik = dcomplex(0., -2.*!pi*dhm1.nm*dhm1.mpp/dhm1.lambda)
f1 = dhm1.alpha * exp(ik * dhm1.zp) * dhm1.field
f2 = dhm2.alpha * exp(ik * dhm2.zp) * dhm2.field
help,f1,f2
xi = 2. * total(real_part(f1 * conj(f2)), 1)

; hologram of sphere pair
return, b1 + b2 + xi - 1.d

end
