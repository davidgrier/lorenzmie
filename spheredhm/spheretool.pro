;+
; NAME:
;     spheretool
;
; PURPOSE:
;     Interactively find reasonable fitting parameters for digital
;     holographic microscopy images of colloidal spheres.
;
; CATEGORY:
;     Digital holography; microscopy
;
; CALLING SEQUENCE:
;     IDL> spheretool, a
;
; INPUTS:
;     a: DHM image of a sphere normalized by a background image.
;
; KEYWORD PARAMETERS:
;     lambda: vacuum wavelength of light [micrometers]
;          Default: 0.6328 (He-Ne)
;
;     mpp: length scale calibration [micrometers/pixel]
;          Default: 0.135
;
;     nm: complex refractive index of medium.
;          Default: 1.3326 (water)
;
; KEYWORD FLAGS:
;     gpu: If set, and if GPU-accelerated procedures are
;         available and working, use GPU accelerated calculations.
;
; SIDE EFFECTS:
;     Opens an interactive widget-based application.
;
; PROCEDURE:
;     Computes DHM image using results for the Mie scattering of
;     coherent light by uniform spheres.  
;     Calls SPHEREDHMPROFILE and FITSPHEREDHM, 
;     which in turn rely on SPHEREFIELD.
;
; REFERENCES:
;  1. S. Lee, Y. Roichman, G. Yi, S. Kim, S. Yang, A. van Blaaderen,
;     P. van Oostrum and D. G. Grier,
;     "Chararacterizing and tracking single colloidal particles 
;      with video holographic microscopy,"
;     Optics Express 15, 18275-18282 (2007).
;
;  2. F. C. Cheong, B. Sun, R. Dreyfus, J. Amato-Grill, K. Xiao
;     and D. G. Grier,
;     "Flow visualization and flow cytometry with
;      holographic video microscopy,"
;     Optics Express 17, 13071-13079 (2009).
;
; INSTRUCTIONS:
;  Click on center of sphere's image in IMAGE window.
;  Click on PROFILE tab to examine the hologram's radial profile.
;  Adjust the radius, refractive index, axial displacement and
;  amplitude to get a good fit.  Click on plot to run
;  least-squares fit.
;
; MODIFICATION HISTORY:
; Written 06/11/2007 by David G. Grier, New York University.
; 06/13/2007 DGG: Save parameters, minimal idiot proofing, minimal
;     documentation.
; 06/18/2007 DGG: Added COHERENCE slider.
; 06/19/2007 DGG: Added standard deviations to the radial plot
; 11/03/2007 DGG: Major additions and internal code cleanup:
;     IMAGE tab: added XC, YC, and LAMBDA.
;     PROFILE tab: Reorganized widgets.
;         Implemented fitting with FITSPHEREDHM.  Fit starts with
;         click on profile plot.  Added checkboxes to fix parameters
;         during fitting.  Widgets updated with fit values.
;     Corrected vacuum wavelength for HeNe laser.  Removed unused
;     internal variables.
; 02/08/2008 DGG: Added "Save Procedure" to File menu.
; 04/16/2008 DGG: Added MPP keyword.  Updated documentation.
; 10/08/2008 DGG: Corrected GUI action of XC, YC and LAMBDA widgets.
; 01/15/2009 DGG: Revised cropping code to allow for spheres closer
;     than RAD to the edge of the image.
; 03/18/2009 DGG: Added MPP widget to image tab.
; 03/18/2009 DGG: Code reorganization.  Separate into component
; files. 
; 03/18/2009 DGG: Support for complex refractive indexes
; 03/20/2009 DGG: Added automatic GPU support.  Added COMPILE_OPT.
;    Pass pointer to state variable rather than the entire structure.
; 06/10/2010 DGG: Compiled routines into a single source file.
; 06/29/2010 DGG: Added GPU flag.
; 12/22/2010 DGG: Fixed bug when nm specified on command line.
; 06/02/2012 DGG Streamlined deinterlacing
;
; Copyright (c) 2007-2012 David G. Grier.
;-

;;;
;;; SPHERETOOL_WIDGETS
;;; 
;;; Creates the SPHERETOOL widget hierarchy
;;;
function spheretool_widgets, p, width, height

COMPILE_OPT IDL2, HIDDEN

; The base widget: Contains the controls and the tabbed windows.
base = WIDGET_BASE(TITLE = 'Sphere Tool', MBAR = bar, $
                   UVALUE = 'base', /COLUMN)

file_menu = WIDGET_BUTTON(bar, VALUE = 'File', /MENU)
file_bn = WIDGET_BUTTON(file_menu, VALUE = '&Save Values...', $
                        UVALUE = 'SAVEVALUES', ACCELERATOR = "Alt+S")
file_bn = WIDGET_BUTTON(file_menu, VALUE = 'Save &Procedure...', $
                        UVALUE = 'SAVEPROGRAM', ACCELERATOR = "Alt+P")
file_bn = WIDGET_BUTTON(file_menu, VALUE = '&Quit', $
                        UVALUE = 'QUIT', ACCELERATOR = "Alt+Q")

; tabbed windows showing results
wtab = WIDGET_TAB(base)

; hologram tab: Show image
wimagetab = WIDGET_BASE(wtab, TITLE = 'IMAGE', /COLUMN)
wimage = WIDGET_DRAW(wimagetab, $
                     /BUTTON_EVENTS, UVALUE = 'IMAGE', $
                     xsize = width, ysize = height, RETAIN = 2)

void = WIDGET_LABEL(wimagetab, $
                    VALUE = 'Particle Position [pixels]', $
                    /ALIGN_LEFT)
wparticle = WIDGET_BASE(wimagetab, COLUMN = 3, /FRAME, /GRID_LAYOUT)
wxc = CW_FIELD(wparticle, TITLE = 'xc', UVALUE = 'XC', /RETURN_EVENTS, $
               /FLOATING, VALUE = p.rc[0])
wyc = CW_FIELD(wparticle, TITLE = 'yc', UVALUE = 'YC', /RETURN_EVENTS, $
               /FLOATING, VALUE = p.rc[1])
wrad = CW_FIELD(wparticle, TITLE = 'radius', UVALUE = 'RAD', /RETURN_EVENTS, $
                /INTEGER, VALUE = p.rad)

void = WIDGET_LABEL(wimagetab, $
                    VALUE = 'Instrumental Parameters [micrometers]', $
                    /ALIGN_LEFT)
winstrument = WIDGET_BASE(wimagetab, COLUMN = 2, /FRAME, /GRID_LAYOUT)
wlambda = CW_FIELD(winstrument, TITLE = 'lambda', UVALUE = 'LAMBDA', $
                   /RETURN_EVENTS, /FLOATING, VALUE = p.lambda)
wmpp = CW_FIELD(winstrument, TITLE = 'mpp', UVALUE = 'MPP', $
                /RETURN_EVENTS, /FLOATING, VALUE = p.mpp)

; profile tab: Interact with hologram's radial profile
wprofiletab = WIDGET_BASE(wtab, TITLE = 'PROFILE', /COLUMN)
wprofile = WIDGET_DRAW(wprofiletab, $
                       /BUTTON_EVENTS, UVALUE = 'IMAGEEVENT', $
                       xsize = width, ysize = height, RETAIN = 2)

bnames = ["ap", "zp", "np", "nm", "alpha", "kp", "deinterlace"]
wbuttons = CW_BGROUP(wprofiletab, bnames, COLUMN = n_elements(bnames), $
                     /NONEXCLUSIVE, UVALUE = 'BUTTONS', $
                     LABEL_LEFT = 'Fixed parameters', $
                     SET_VALUE = p.fixed)

wcontrols = WIDGET_BASE(wprofiletab, ROW = 3)
wap = CW_FSLIDER(wcontrols, TITLE = "Particle Radius",$
                 FORMAT = '("ap:", G13.6, " [um]")', $
                 MIN = 0.1, MAX = 10., /DRAG, /EDIT, $
                 UVALUE = 'AP', VALUE = p.ap, XSIZE = width/2)

wzp = CW_FSLIDER(wcontrols, TITLE = "Axial Position", $
                 FORMAT = '("zp:", G13.6, " [pixels]")', $
                 MIN = 10., MAX = 500., /DRAG, /EDIT, $
                 UVALUE = 'ZP', VALUE = p.zp, XSIZE = width/2)

wnp = CW_FSLIDER(wcontrols, TITLE = "Particle Refractive Index", $
                 FORMAT = '("np:", G13.6)', $
                 MIN = 1., MAX = 3., /DRAG, /EDIT, $
                 UVALUE = 'NP', VALUE = p.np, XSIZE = width/2)

wkp = CW_FSLIDER(wcontrols, TITLE = "Particle Extinction Coefficient", $
                 FORMAT = '("kp:", G13.6)', $
                 MIN = 0., MAX = 5., /DRAG, /EDIT, $
                 UVALUE = 'KP', VALUE = p.kp, XSIZE = width/2)

wnm = CW_FSLIDER(wcontrols, TITLE = "Medium Refractive Index", $
                 FORMAT = '("nm:", G13.6)', $
                 MIN = 1., MAX = 5., /DRAG, /EDIT, $
                 UVALUE = 'NM', VALUE = p.nm, XSIZE = width/2)

walpha = CW_FSLIDER(wcontrols, TITLE = "Amplitude", $
                    FORMAT = '("alpha:", G13.6)', $
                    MIN = 0., MAX = 3., /DRAG, /EDIT, $
                    UVALUE = 'ALPHA', VALUE = p.alpha, XSIZE = width/2)

; structure containing identifiers for all widgets
w = {base:     base, $          ; the base widget
     wtab:     wtab, $          ; top level widget containing tabbed widgets
     wimage:   wimage, $        ; holographic image
     wprofile: wprofile, $      ; radial profile with controls
     wap:      wap, $           ; particle radius control
     wzp:      wzp, $           ; particle's axial position control
     wnp:      wnp, $           ; particle index control
     wkp:      wkp, $           ; particle extinction coefficient
     wnm:      wnm, $           ; medium index control
     walpha:   walpha, $        ; uniformity
     wxc:      wxc, $           ; particle center, x
     wyc:      wyc, $           ; particle center, y
     wrad:     wrad, $          ; radius of region of interest
     wlambda:  wlambda, $       ; wavelength of light
     wmpp:     wmpp}            ; micrometers per pixel

return, w
end

;;;
;;; SPHERETOOL_GPU_DETECT
;;;
;;; Checks for availability of GPU accelerated computation
;;;
function spheretool_gpu_detect

COMPILE_OPT IDL2, HIDDEN

catch, error
if (error ne 0L) then begin
   catch, /cancel
   return, 0
endif

gpuinit, /HARDWARE

return, 1
end


;;;
;;; SPHERETOOL_SAVEROUTINE
;;;
;;; Metaprogramming routine that writes a program
;;; to analyze holograms based on current settings
;;;
pro spheretool_saveroutine, s

COMPILE_OPT IDL2, HIDDEN

fn = DIALOG_PICKFILE(/WRITE, /OVERWRITE_PROMPT, $
                     DEFAULT_EXTENSION = 'pro', FILTER=['*.pro'], $
                     TITLE = "Save Procedure")
name = file_basename(fn, '.pro', /FOLD_CASE)
if strlen(fn) lt 1 then return
sz = size((*s).p.a,/DIMENSIONS)
rc = floor((*s).p.rc)
rad = min([(*s).p.rad, rc, sz-rc])
deinterlace = (*s).p.fixed[6]
openw, lun, fn, /get_lun
printf, lun, "pro " + name
printf, lun, ";;; routine generated automatically by SPHERETOOL.PRO"
printf, lun, ";;; " + systime()
printf, lun, ";;; Specify the background image"
printf, lun, "backgroundfile = XXX"
printf, lun, ";;; Specify the image files to be processed:"
printf, lun, "filespec = YYY"
printf, lun, ";;; Specify name of output data file"
printf, lun, 'filename = "'+name+'.gdf"'
printf, lun, ";;; nothing should need to be changed below this line"
printf, lun, "bg = double(read_image(backgroundfile)) > 1.d"
printf, lun, "sz = size(bg,/dimensions)"
printf, lun, "f = file_search(filespec, count=nfiles)"
printf, lun, ";;; Region of interest"
printf, lun, format='("rc = [", I4, ",", I4, "]")', rc
printf, lun, "rad =", rad
printf, lun, ";;; Flags"
printf, lun, "fixap =", (*s).p.fixed[0]
printf, lun, "fixzp =", (*s).p.fixed[1]
printf, lun, "fixnp =", (*s).p.fixed[2]
printf, lun, "fixnm =", (*s).p.fixed[3]
printf, lun, "fixalpha =", (*s).p.fixed[4]
printf, lun, ";;; Starting parameters for fits"
printf, lun, "zp =", (*s).p.zp
printf, lun, "ap =", (*s).p.ap
printf, lun, "np =", (*s).p.np
printf, lun, "alpha =", (*s).p.alpha
printf, lun, "nm =", (*s).p.nm
printf, lun, "lambda =", (*s).p.lambda
printf, lun, "mpp =", (*s).p.mpp
printf, lun, "p = dblarr(2,7,nfiles)"
printf, lun, "pin = double([0, 0, zp, ap, np, alpha, nm])"
printf, lun, "for n = 0, nfiles-1 do begin"
printf, lun, "   x0 = round(rc[0]-rad) > 0"
printf, lun, "   x1 = round(rc[0]+rad) < sz[0]-1"
printf, lun, "   y0 = round(rc[1]-rad) > 0"
printf, lun, "   y1 = round(rc[1]+rad) < sz[1]-1"
printf, lun, "   nx = x1 - x0 + 1"
printf, lun, "   ny = y1 - y0 + 1"
printf, lun, "   xp = rc[0] - x0 - nx/2."
printf, lun, "   yp = rc[1] - y0 - ny/2."
printf, lun, "   a = double((read_image(f[n]))[x0:x1,y0:y1])"
printf, lun, "   a /= bg[x0:x1,y0:y1]"         
printf, lun, "   pin[0:1] += [xp,yp]"
printf, lun, "   pout = fitspheredhm(a, pin, lambda=lambda, mpp=mpp, $"
if deinterlace then $
   printf, lun, "     deinterlace = y0, $"
printf, lun, "     fixap=fixap,fixzp=fixzp,fixnm=fixnm,fixalpha=fixalpha)"
printf, lun, "   pout[0,0:1] += [x0+nx/2., y0+ny/2]"
printf, lun, "   rc = pout[0,0:1]"
printf, lun, "   pin[2:*] = pout[0,2:*]"
printf, lun, "   p[*,*,n] = pout"
printf, lun, "   print, n, pout"
printf, lun, "   write_gdf, p, filename"
printf, lun, "endfor"
printf, lun, "end"
close, lun
free_lun, lun
end


;;;
;;; SPHERETOOL_CLEANUP
;;;
;;; Free resources at shutdown
;;;
pro spheretool_cleanup, w

COMPILE_OPT IDL2, HIDDEN

WIDGET_CONTROL, w, GET_UVALUE = s
PTR_FREE, s

end

;;;
;;; SPHERETOOL_EVENT
;;;
;;; Process events generated by the SPHERETOOL GUI.
;;;
pro spheretool_event, ev

COMPILE_OPT IDL2, HIDDEN

; The base widget's uvalue is the structure containing all the
; information about the program's state.
WIDGET_CONTROL, ev.TOP, GET_UVALUE = s

; Clicking on a tab is handled automatically.  Nothing to do, so
; move along.
if (ev.ID eq (*s).w.wtab) then return

; Other widgets can give us something to do. 
WIDGET_CONTROL, ev.ID, GET_UVALUE = uval

sz = size((*s).p.a, /DIMENSIONS)

CASE uval of
    'IMAGE' : begin
        if (ev.press) then begin
            !x.s = (*s).geom.xs
            !y.s = (*s).geom.ys
            val = convert_coord(ev.x, ev.y, /device, /to_data)
            (*s).p.rc = val[0:1]
            if (*s).p.fixed[6] then begin
               (*s).p.aa = aziavg((*s).p.a, center = (*s).p.rc, /deinterlace)
               (*s).p.daa = azistd((*s).p.a, center = (*s).p.rc, /deinterlace)
            endif else begin
               (*s).p.aa = aziavg((*s).p.a, center = (*s).p.rc)
               (*s).p.daa = azistd((*s).p.a, center = (*s).p.rc)
            endelse
            WIDGET_CONTROL, (*s).w.wxc, SET_VALUE = val[0]
            WIDGET_CONTROL, (*s).w.wyc, SET_VALUE = val[1]
        endif
    end

    'BUTTONS' : begin
       (*s).p.fixed[ev.value] = ev.select
    end

    'AP' : begin
        WIDGET_CONTROL, (*s).w.wap, GET_VALUE = val
        (*s).p.ap = val
    end

    'ZP' : begin
       WIDGET_CONTROL, (*s).w.wzp, GET_VALUE = val
       (*s).p.zp = val
    end

    'NP' : begin
       WIDGET_CONTROL, (*s).w.wnp, GET_VALUE = val
       (*s).p.np = val
    end

    'KP': begin
       WIDGET_CONTROL, (*s).w.wkp, GET_VALUE = val
       (*s).p.kp = val
    end

    'NM' : begin
       WIDGET_CONTROL, (*s).w.wnm, GET_VALUE = val
       (*s).p.nm = val
    end

    'ALPHA': begin
        WIDGET_CONTROL, (*s).w.walpha, GET_VALUE = val
        (*s).p.alpha = val
    end

    'XC': begin
       WIDGET_CONTROL, (*s).w.wxc, GET_VALUE = val
       (*s).p.rc[0] = val
       if (*s).p.fixed[6] then begin
          (*s).p.aa = aziavg((*s).p.a, center = (*s).p.rc, /deinterlace)
          (*s).p.daa = azistd((*s).p.a, center = (*s).p.rc, /deinterlace)
       endif else begin
          (*s).p.aa = aziavg((*s).p.a, center = (*s).p.rc)
          (*s).p.daa = azistd((*s).p.a, center = (*s).p.rc)
       endelse
    end

    'YC': begin
       WIDGET_CONTROL, (*s).w.wyc, GET_VALUE = val
       (*s).p.rc[1] = val
       if (*s).p.fixed[6] then begin
          (*s).p.aa = aziavg((*s).p.a, center = (*s).p.rc, /deinterlace)
          (*s).p.daa = azistd((*s).p.a, center = (*s).p.rc, /deinterlace)
       endif else begin
          (*s).p.aa = aziavg((*s).p.a, center = (*s).p.rc)
          (*s).p.daa = azistd((*s).p.a, center = (*s).p.rc)
       endelse
    end

    'RAD': begin
       WIDGET_CONTROL, (*s).w.wrad, GET_VALUE = val
       (*s).p.rad = val
     end

    'LAMBDA': begin
       WIDGET_CONTROL, (*s).w.wlambda, GET_VALUE = val
       (*s).p.lambda = val
    end

    'MPP': begin
       WIDGET_CONTROL, (*s).w.wmpp, GET_VALUE = val
       (*s).p.mpp = val
    end

    'IMAGEEVENT': begin
       if ev.press then begin
          print, "Fitting to Mie scattering theory ..."
          WIDGET_CONTROL, /HOURGLASS
          rc = (*s).p.rc
          rad = (*s).p.rad
          x0 = round(rc[0] - rad) > 0
          x1 = round(rc[0] + rad) < sz[0]-1
          y0 = round(rc[1] - rad) > 0
          y1 = round(rc[1] + rad) < sz[1]-1
          ac = (*s).p.a[x0:x1, y0:y1]
          nx = x1 - x0 + 1
          ny = y1 - y0 + 1
          xp = rc[0] - x0 - nx/2.
          yp = rc[1] - y0 - ny/2.
          pin = [xp, yp, (*s).p.zp, (*s).p.ap, (*s).p.np, (*s).p.kp, (*s).p.nm, (*s).p.km, $
                 (*s).p.alpha, (*s).p.delta]
          if (*s).p.fixed[6] then begin
             plotimage, bytscl(deinterlace(ac, /odd)), /iso
             pout = fitspheredhm(ac, pin, $
                                 lambda = (*s).p.lambda, $
                                 mpp = (*s).p.mpp, $
                                 deinterlace = y0, $
                                 fixap = (*s).p.fixed[0], $
                                 fixzp = (*s).p.fixed[1], $
                                 fixnp = (*s).p.fixed[2], $
                                 fixnm = (*s).p.fixed[3], $
                                 fixalpha = (*s).p.fixed[4], $
                                 gpu = (*s).p.gpu)
          endif else begin
             plotimage, bytscl(ac), /iso
             pout = fitspheredhm(ac, pin, $
                                 lambda = (*s).p.lambda, $
                                 mpp = (*s).p.mpp, $
                                 fixap = (*s).p.fixed[0], $
                                 fixzp = (*s).p.fixed[1], $
                                 fixnp = (*s).p.fixed[2], $
                                 fixnm = (*s).p.fixed[3], $
                                 fixalpha = (*s).p.fixed[4], $
                                 gpu = (*s).p.gpu)
          endelse
          print, pout
          (*s).p.rc = [x0+nx/2., y0+ny/2.] + pout[0, 0:1] 
          WIDGET_CONTROL, (*s).w.wxc, SET_VALUE = (*s).p.rc[0]
          WIDGET_CONTROL, (*s).w.wyc, SET_VALUE = (*s).p.rc[1]
          (*s).p.zp = pout[0, 2]
          WIDGET_CONTROL, (*s).w.wzp, SET_VALUE = (*s).p.zp
          (*s).p.ap = pout[0, 3]
          WIDGET_CONTROL, (*s).w.wap, SET_VALUE = (*s).p.ap
          (*s).p.np = pout[0, 4]
          WIDGET_CONTROL, (*s).w.wnp, SET_VALUE = (*s).p.np
          (*s).p.kp = pout[0, 5]
          (*s).p.nm = pout[0, 6]
          WIDGET_CONTROL, (*s).w.wnm, SET_VALUE = (*s).p.nm
          (*s).p.km = pout[0, 7]
          (*s).p.alpha = pout[0, 8]
          WIDGET_CONTROL, (*s).w.walpha, SET_VALUE = (*s).p.alpha
          (*s).p.delta = pout[0, 9]

          (*s).p.aa = aziavg((*s).p.a, center = (*s).p.rc)
          (*s).p.daa = azistd((*s).p.a, center = (*s).p.rc)
          WIDGET_CONTROL, HOURGLASS = 0
       endif
    end

    'SAVEPROGRAM': begin
       spheretool_saveroutine, s
    end

    'SAVEVALUES': begin
        fn = DIALOG_PICKFILE(/WRITE, /OVERWRITE_PROMPT, $
                             TITLE = "Save Values", FILTER = "*.gdf", $
                             DEFAULT_EXTENSION = "gdf")
        if strlen(fn) ge 1 then begin
           if (*s).p.fit[1,0] gt 0. then $
              write_gdf, (*s).p.fit, fn, /ASCII $
           else begin
              write_gdf, [(*s).p.rc, $
                          (*s).p.zp, $
                          (*s).p.ap, $
                          (*s).p.np, $
                          (*s).p.kp, $
                          (*s).p.nm, $
                          (*s).p.km, $
                          (*s).p.alpha, $
                          (*s).p.delta], fn, /ASCII
           endelse
        endif
    end

    'QUIT' : begin
        WIDGET_CONTROL, ev.TOP, /DESTROY
        return
    end

    ELSE: print, 'Event type not yet implemented'
ENDCASE

; Process the present image according to the present settings.

; Show updated region of interest
WIDGET_CONTROL, (*s).w.wimage, GET_VALUE = ndx
wset, ndx
dev = !d.name
set_plot, 'z'
plotimage, bytscl((*s).p.a, top=253), /iso
rc = (*s).p.rc
rad = (*s).p.rad
x0 = round(rc[0] - rad) > 0
x1 = round(rc[0] + rad) < sz[0]-1
y0 = round(rc[1] - rad) > 0
y1 = round(rc[1] + rad) < sz[1]-1
plots, rc[0], rc[1], psym=circ()
plots, [x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0], linestyle = 3, color = 254
image = tvrd()
set_plot, dev
tv, image

; show update DHM profile
WIDGET_CONTROL, (*s).w.wprofile, GET_VALUE = ndx
wset, ndx
rho = findgen(n_elements((*s).p.aa))
b = spheredhmprofile(rho, (*s).p.zp, (*s).p.ap, dcomplex((*s).p.np, (*s).p.kp), $
                     dcomplex((*s).p.nm, (*s).p.km), (*s).p.alpha, (*s).p.delta, $
                     (*s).p.lambda, (*s).p.mpp)
set_plot, 'z'
plot, (*s).p.aa, thick=2, xtitle='r [pixel]', ytitle='B(r)', $
  xrange = [0,(*s).p.rad], /xstyle
oplot, (*s).p.aa + (*s).p.daa, linestyle=2
oplot, (*s).p.aa - (*s).p.daa, linestyle=2
oplot, rho, b, color = 254
image = tvrd()
set_plot, dev
tv, image

WIDGET_CONTROL, ev.TOP, SET_UVALUE = s
end

;;;
;;; SPHERETOOL
;;;
;;; The main routine
;;;
pro spheretool, a, nm = nm, lambda = lambda, mpp = mpp, gpu = gpu

COMPILE_OPT IDL2

width = 640
height = 480

if n_params() ne 1 then begin
    message, "USAGE: spheretool, A", /inf
    message, "A must be a normalized holographic microscope image.", /inf
    return
endif

sz = size(a)
if sz[0] ne 2 then begin
    message, "USAGE: IDL> spheretool, A", /inf
    message, "A must be a two-dimensional normalized hologram", /inf
    return
endif

aa = aziavg(a)
daa = azistd(a)
xc = sz[1]/2.
yc = sz[2]/2.
rad = floor(xc < yc)
fit = dblarr(2,7)

; p contains all the parameters used by SPHEREDHM
; Storing them here avoids the need for common blocks,
; and thus allows for multiple instances of SPHERETOOL.
p = {a:      a, $                  ; the image
     aa:     aa, $                 ; azimuthally averaged image
     daa:    daa, $                ; azimuthal standard deviation
     fit:    fit, $                ; results of fits
     rc:     [xc,yc], $            ; in-plane center [pixels]
     rad:    rad, $                ; radius of ROI
     lambda: 0.632816d, $          ; wavelength of light [micrometers] (HeNe)
     mpp:    0.135d, $             ; length calibration [micrometers/pixel]
     fixed:  [0,0,0,1,1,1,0], $    ; flags for fixed parameters
     zp:     100.d, $              ; axial position [pixels]
     ap:     1.d, $                ; particle radius [micrometers]
     np:     1.5d, $               ; particle refractive index
     kp:     0.d, $                ; particle extinction coefficient
     nm:     1.3326d, $            ; medium refractive index (water)
     km:     0.d, $                ; medium extinction coefficient
     alpha:  1.d, $                ; illumination amplitude
     delta:  0.d, $                ; wavefront distortion [pixels]
     gpu:    0}                    ; gpu acceleration flag

if keyword_set(gpu) then $
   p.gpu = spheretool_gpu_detect()

if n_elements(nm) eq 1 then begin
   p.nm = real_part(nm)
   p.km = imaginary(nm)
endif

if n_elements(lambda) eq 1 then $
   p.lambda = double(lambda)

if n_elements(mpp) eq 1 then $
   p.mpp = double(mpp)

red   = indgen(255)
green = red
blue  = red
red[254] = 0
blue[254] = 0
tvlct, red, green, blue

w = spheretool_widgets(p, width, height)

WIDGET_CONTROL, w.base, /REALIZE

; show image with region of interest
WIDGET_CONTROL, w.wimage, GET_VALUE = ndx
wset, ndx
plotimage, bytscl(p.a, top=253), /iso
plots, p.rc[0], p.rc[1], psym=circ()
x0 = xc - rad > 0
x1 = xc + rad < sz[1] - 1
y0 = yc - rad > 0
y1 = yc + rad < sz[2] - 1
plots, [x0,x1,x1,x0,x0], [y0,y0,y1,y1,y0], linestyle=3, color=254
geom = {xs:!x.s, ys:!y.s}

; show radial profile
WIDGET_CONTROL, w.wprofile, GET_VALUE = ndx
wset, ndx
rho = findgen(100)
b = spheredhmprofile(rho, p.zp, p.ap, dcomplex(p.np, p.kp), $
                     dcomplex(p.nm, p.km), p.alpha, p.delta, $
                     p.lambda, p.mpp)
plot, p.aa, thick = 2, $
      xtitle = 'r [pixel]', ytitle = 'B(r)', $
      xrange=[0, p.rad], /xstyle
oplot, rho, b, color = 254

; program state variables
s = {w:w, geom:geom, p:p}
ps = PTR_NEW(s, /NO_COPY)

; transfer state variables to callback routine to process events
WIDGET_CONTROL, w.base, SET_UVALUE = ps, /NO_COPY

; process events
XMANAGER, 'spheretool', w.base, /NO_BLOCK, CLEANUP='spheretool_cleanup'

end
