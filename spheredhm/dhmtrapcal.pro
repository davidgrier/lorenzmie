function dhmtrapcal, r, drsq

; r: [3,npts] array of particle positions
; drsq: mean-square measurement error

nsteps = n_elements(r[0,*])

; mean (equilibirium) position
rc = total(r,2)/nsteps

for t = 0, nsteps-1 do $
    r[*,t] -= rc

c0 = total(r^2,2)/nsteps
c1 = total(r * shift(r,0,1),2)/nsteps
c2 = total(r * shift(r,0,2),2)/nsteps

w = where(c0 lt drsq, nsmall)
if nsmall gt 0 then $
  message, "measurement error exceeds displacement in" + $
  string(nsmall) + " trajectories",/inf

k = 1.D/c0 
k += drsq * (c0^4 - 2.D*c0^2*c1^2 - 3.D*c1^4 + 2.D*c0*c1^2*c2) / $
     (c0^3 - c0*c1^2)^2

g = -alog(c1/c0)/c0
g += drsq * ((c0^2 - c1^2)*(c0^2 - c1^2 + c0*c2) + $
           (c0^4 - 2.D*c0^2*c1^2 - 3.D*c1^4 + 2.D*c0*c1^2*c2) * alog(c1/c0)) / $
     (c0^3-c0*c1^2)^2

          
dk = sqrt(2.D/nsteps)*sqrt((c0^2 + c1^2)/(c0^2 - c1^2))
dk += (drsq/c0) * 2.D*(c0^6 + c0^4*c1^2 - 9.D*c0^2*c1^4 + 7.D*c1^6 + 6.D*c0^3*c1^2*c2 - 4.D*c0*c1^4*c2) / $
      ((c0^2 - c1^2)^3 * sqrt(nsteps/2.D) * sqrt((c0^2 + c1^2)/(c0^2 - c1^2)))                
dk *= k

dg = sqrt(1.D/nsteps) * sqrt(2.D + $
                             (-c0^2 + c1^2 + 2.D*c1^2*alog(c1/c0))^2 / $
                             (c1^2*(c0^2 - c1^2)*(alog(c1/c0))^2))
dg += (drsq/c0) * ((c0^2 - c1^2)^3*(-c1^2 + c0*(c0 + c2)) + $
                  2.D*c1^2*alog(c1/c0) * $
                  (c0*(c0^2 - c1^2)^2*c2 - 2.D*(c0^2 - c1^2) * $
                   (c0^4 - 3.D*c0^2*c1^2 + 2.D*c1^4 + c0^3*c2)*alog(c1/c0) + $
                   2.D*(c0^6 + c0^4*c1^2 - 9.D*c0^2*c1^4 + 7.D*c1^6 + 6.D*c0^3*c1^2*c2 - $
                        4.D*c0*c1^4*c2)*(alog(c1/c0))^2)) / $
      (c0*c1^2*(c0^2 - c1^2)^3*sqrt(nsteps)*(alog(c1/c0))^3* $
       sqrt(2.D + (-c0^2 + c1^2 + 2.D*c1^2*alog(c1/c0))^2/(c1^2*(c0^2 - c1^2)*(alog(c1/c0))^2)))
dg *= g
           
return,[[k],[dk],[g],[dg]]
end
