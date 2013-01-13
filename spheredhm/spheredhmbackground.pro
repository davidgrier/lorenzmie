function spheredhmbackground, a0, p0, rad = rad

; a0: image
; p0: starting parameters

a = double(a0)
dim = size(a,/dimensions)

ma = mean(a) ; starting guess for background subtraction

; estimate sphere's image
p = fitspheredhm(a/ma, p0, lambda=lambda, /fixnm)
rp = reform(p[0,0:2])
ap = reform(p[0,3])
np = reform(p[0,4])
nm = reform(p[0,5])
alpha = reform(p[0,6])

sphere = spheredhm(rp, ap, np, nm, dim, alpha=alpha)

background = a - ma * sphere + 1.

; better?

return, b
end
