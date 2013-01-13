a = DGGdhmSphereDHM(dim = [201,201], lambda = 0.640, mpp = 0.135d, zp = 100)
a.xp += 20
a.yp += 10
b = DGGdhmSphereDHM(dim = [201,201], lambda = 0.640, mpp = 0.135d, zp = 100.)
objanswer = twospheredhm(a, b)

; Construct fake hologram of two particles, each 10 pixels from the center in the x direction.
dimen = [201,201]

x = rebin(findgen(dimen[0]) - dimen[0]/2., dimen[0], dimen[1])
y = findgen(1, dimen[1]) - dimen[1]/2.
y = rebin(y, dimen[0], dimen[1])
field1 = spherefield(x-a.xp, y-a.yp, a.zp, a.ap, a.np, nm = a.nm, $
                  mpp = a.mpp, lambda = a.lambda, $
                  k = k, /cartesian, gpu=gpu)
                  
field2 = spherefield(x-b.xp, y-b.yp, b.zp, b.ap, b.np, nm = b.nm, $
                  mpp = b.mpp, lambda = b.lambda, $
                  k = k, /cartesian, gpu=gpu)

xhat = 0.*field1
xhat[0,*] = 1.

answer = xhat + field1*exp(dcomplex(0,-k*a.zp)) + field2*exp(dcomplex(0,-k*b.zp))
answer = total(real_part(answer*conj(answer)), 1)
answer = reform(answer, 201, 201)
