function trackspheredhm, filespec, bg, p0, rad=rad, nm=nm, rc=rc

f = file_search(filespec, count=nfiles)

a = read_image(f[0])            ; raw hologram
b = (float(a) - bg)/sqrt(bg)          ; normalized

if n_elements(rad) lt 1 then $
  rad = 50

dim = 2*rad+1

;if n_elements(nm) lt 1 then $
;  nm = 1.33 ; water
nm = p0[6]

sz = size(p0)
if sz[0] eq 1 then begin
    if n_elements(rc) ne 2 then begin
; estimate center of scattering pattern
        plotimage, bytscl(b), /iso
        cursor, x, y, /data
        xc = x
        yc = y
        rc = [xc,yc]
    endif else begin
        xc = rc[0]
        yc = rc[1]
    endelse
    xc = round(xc)
    yc = round(yc)

; establish region of interest
    x0 = xc-rad
    x1 = xc+rad
    y0 = yc-rad
    y1 = yc+rad

    a = a[x0:x1,y0:y1]
    b = bg[x0:x1,y0:y1]
    sb = sqrt(b)

    d = (a-b)/sb
    tvscl,d

; fit the first one
    p = fitspheredhm(d,p0,/quiet,/fixnm)
    dd = p[0,5]*spheredhm(p[0,0:2],p[0,3],p[0,4],nm,[dim,dim])
    tvscl,dd,dim,0
    chi = sqrt(mean((d-dd)^2))
    print, "Chi-squared =", chi^2
    print, transpose(p[0,*])
    p[1,*] *= chi
    dx = p[0,0]
    dy = p[0,1]
    p[0,0] += xc
    p[0,1] += yc

    res = p
    npts = 1
endif else if sz[0] eq 3 then begin
    npts = sz[3]
    res = p0
    p = reform(p0[*,*,npts-1])
    xc = p[0,0]
    yc = p[0,1]
    dx = 0
    dy = 0
endif else $
  message, "incorrect format for initial parameters"
    
for i = npts, nfiles-1 do begin
    xc = round(xc + dx)
    yc = round(yc + dy)
    x0 = xc-rad
    x1 = xc+rad
    y0 = yc-rad
    y1 = yc+rad

    a = read_image(f[i])
    a = a[x0:x1,y0:y1]
    b = bg[x0:x1,y0:y1]
    sb = sqrt(b)

    d = (a-b)/sb
    tvscl,d
    p = fitspheredhm(d,[0.,0.,transpose(p[0,2:*])],/quiet,/fixnm)
    dd = p[0,5]*spheredhm(p[0,0:2],p[0,3],p[0,4],nm,[dim,dim])
    tvscl,dd,dim,0
    print, transpose(p[0,*])
    p[1,*] *= chi
    dx = p[0,0]
    dy = p[0,1]
    p[0,0] += xc
    p[0,1] += yc

    res = [[[res]],[[p]]]
    write_gdf, res, "trackspheredhm.gdf"
endfor

return, res
end
