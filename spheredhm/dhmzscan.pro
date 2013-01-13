function dhmzscan, a_, ap, np, nm, lambda=lambda, mpp=mpp

minz = 100.
maxz = 450.
dz = 1.

alpha = 1.

npts = n_elements(a_)
x = findgen(npts)
y = fltarr(npts)

a = a_ - 1.

err = total(a^2)
zp = 0.

for z = minz, maxz, 1 do begin
   field = spherefield(x, y, z, ap, np, k=k, $
                       nm=nm, lambda=lambda, mpp=mpp, /cartesian, $
                       precision=precision, gpu=gpu)
   dhm = 2.d * alpha * field[0,*] * exp(dcomplex(0.,-k*zp))
   dhm += alpha^2 * total(field*conj(field),1)
   dhm = real_part(dhm)
   thiserr = total((dhm - a)^2)
;   print, z, thiserr
   if thiserr lt err then begin
      zp = z
      err = thiserr
   endif
endfor

return, z
end


