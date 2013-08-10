;+
; NAME:
;    lmtool
;
; PURPOSE:
;    Interactively guided Lorenz-Mie analysis of inline holograms
;
; CATEGORY:
;    Holographic video microscopy
;
; CALLING SEQUENCE:
;    lmtool, hologram, lambda, mpp
;
; INPUTS:
;    hologram: normalized hologram
;
;    lambda: wavelength of illumination [micrometers]
;
;    mpp: magnification [micrometer/pixel]
;
; KEYWORD PARAMETERS:
;    nm: Complex refractive index of medium
;
; KEYWORD FLAG:
;    deinterlace: 0: no deinterlacing
;                 1: analyze odd field
;                 2: analyze even field
;
; SIDE EFFECTS:
;    Creates GUI
;
; PROCEDURE:
;    Uses lmtool procedures and functions to analyze in-line hologram.
;
; NOTES:
;    View residuals
;    Implement and test GPU
;
; MODIFICATION HISTORY:
; 01/29/2013 Written by David G. Grier, New York University
; 02/03/2013 DGG Initial versions of save procedures.  Increase range of
;   ZP to 300 pixels.
; 07/24/2013 DGG Use GPU acceleration.
; 08/05/2013 DGG Set radius with fringe number.
; 08/10/2013 DGG Added RESIDUALS tab.
;
; Copyright (c) 2013 David G. Grier
;-

;;;;;
;
; LMT_FIT
;
pro lmt_fit, ev

COMPILE_OPT IDL2, HIDDEN

widget_control, ev.top, get_uvalue = s
widget_control, /hourglass
h = (*s).p.h                      ; the hologram
f = h -> get(position = (*s).p.n) ; selected feature
fix = (*s).p.fix                  ; fixed parameters
f.fit, fixap = fix.ap, fixnp = fix.np, fixkp = fix.kp, $
       fixzp = fix.zp, fixnm = fix.nm, fixkm = fix.km, $
       fixalpha = fix.alpha, fixdelta = fix.delta
lmt_update_fit_plot, s
lmt_update_res_plot, s
lmt_update_parameter_widgets, s
lmt_update_profile_plot, s
widget_control, hourglass = 0
end

;;;;;
;
; LMT_SELECT_FEATURE()
;
function lmt_select_feature, win, x, y, button, keymods, clicks

COMPILE_OPT IDL2, HIDDEN

s = win.uvalue
objlist = win.hittest(x, y)
foreach obj, objlist do begin
   if isa(obj, 'polygon') then begin
      flist = (*s).w.f
      flist[(*s).p.n].update
      (*s).p.n = obj.uvalue
      flist[(*s).p.n].update, /selected
      lmt_update_parameter_widgets, s
      lmt_update_profile_plot, s
      lmt_update_fit_plot, s
      lmt_update_res_plot, s
      break
   endif
endforeach

return, 0                       ; do not call the default event handler
end

;;;;;
;
; LMT_UPDATE_PARAMETER_WIDGETS
;
pro lmt_update_parameter_widgets, s

COMPILE_OPT IDL2, HIDDEN

f = (*s).p.h -> get(position = (*s).p.n)
widget_control, (*s).w.ap, set_value = f.ap
widget_control, (*s).w.zp, set_value = (f.rp)[2]
widget_control, (*s).w.np, set_value = real_part(f.np)
widget_control, (*s).w.kp, set_value = imaginary(f.np)
widget_control, (*s).w.alpha, set_value = f.alpha
widget_control, (*s).w.delta, set_value = f.delta
widget_control, (*s).w.rad, set_value = f.rad

end

;;;;;
;
; LMT_UPDATE_FIT_PLOT
;
pro lmt_update_fit_plot, s, compare = compare

COMPILE_OPT IDL2, HIDDEN

h = (*s).p.h
f = h.get(position = (*s).p.n) ; selected feature
if keyword_set(compare) then begin
   (*s).w.hologram -> setdata, deinterlace(f.data, h.deinterlace)-f.lmhologram
endif else begin
   (*s).w.imhologram -> setdata, deinterlace(f.data, h.deinterlace)
   (*s).w.imfit -> setdata, f.lmhologram
endelse
end

;;;;;
;
; LMT_UPDATE_RES_PLOT
;
pro lmt_update_res_plot, s

COMPILE_OPT IDL2, HIDDEN

h = (*s).p.h
f = h.get(position = (*s).p.n) ; selected feature
(*s).w.imres -> setdata, deinterlace(f.data, h.deinterlace)-f.lmhologram
end

;;;;;
;
; LMT_UPDATE_PROFILE_PLOT
;
pro lmt_update_profile_plot, s

COMPILE_OPT IDL2, HIDDEN

h = (*s).p.h
f = h.get(position = (*s).p.n)  ; selected feature
profile = f.profile
dprofile = f.dprofile
lmprofile = f.lmprofile
pplot = (*s).w.pplot
pplot.xrange = [0., f.rad]
pplot -> setdata, profile
x = findgen(f.rad+1)
ppoly = (*s).w.ppoly
ppoly -> setdata, [x, reverse(x)], [profile+dprofile, reverse(profile-dprofile)]
pfit = (*s).w.pfit
pfit -> setdata, lmprofile

end

;;;;;
;
; LMT_DEINTERLACE
;
; Handle changes in deinterlace choice
;
pro lmt_deinterlace, s, deinterlace

COMPILE_OPT IDL2, HIDDEN

h = (*s).p.h
h.deinterlace = deinterlace
(*s).w.im -> putdata, deinterlace(h.data, h.deinterlace)
lmt_update_profile_plot, s
lmt_update_fit_plot, s
lmt_update_res_plot, s

end

;;;;;
;
; LMT_UPDATE_FIX_FLAGS
;
; Handle changes in fix-parameter flags
;
pro lmt_update_fix_flags, s, flags

COMPILE_OPT IDL2, HIDDEN

(*s).p.fix.ap = flags[0]
(*s).p.fix.zp = flags[1]
(*s).p.fix.np = flags[2]
(*s).p.fix.kp = flags[3]
(*s).p.fix.alpha = flags[4]
(*s).p.fix.delta = flags[5]
(*s).p.fix.nm = flags[6]
(*s).p.fix.km = flags[7]
end

;;;;;
;
; LMT_UPDATE_PARAMETER
;
pro lmt_update_parameter, ev

COMPILE_OPT IDL2, HIDDEN

widget_control, ev.id, get_uvalue = uval, get_value = val
widget_control, ev.top, get_uvalue = s
h = (*s).p.h
f = h -> get(position = (*s).p.n)
case uval of
   'ZP': f.zp = val
   'AP': f.ap = val
   'NP': f.np = dcomplex(val, imaginary(f.np))
   'KP': f.np = dcomplex(real_part(f.np), val)
   'ALPHA': f.alpha = val
   'DELTA': f.delta = val
   'NM': h.nm = dcomplex(val, imaginary(h.nm))
   'KM': h.nm = dcomplex(real_part(h.nm), val)
   'LAMBDA': h.lambda = val
   'MPP' : h.mpp = val
   'DEINTERLACE' : lmt_deinterlace, s, val
   'FIX' : lmt_update_fix_flags, s, val
   'NFRINGES' : begin
      h.nfringes = val
      foreach feature, (*s).w.f do $
         feature.update
      ((*s).w.f)[(*s).p.n].update, /selected
      lmt_update_profile_plot, s
      lmt_update_parameter_widgets, s
   end
   'RAD' : begin
      f.rad = val
      ((*s).w.f)[(*s).p.n].update, /selected
      lmt_update_profile_plot, s
   end
   'ESTIMATE': begin
      fix = (*s).p.fix
      f.fitprofile, fixap = fix.ap, fixnp = fix.np, fixzp = fix.zp, $
                            fixalpha = fix.alpha, fixdelta = fix.delta
      lmt_update_parameter_widgets, s
   end
   else: message, uval, /inf
endcase

(*s).w.pfit -> setdata, f.lmprofile

end

;;;;;
;
; LMT_EVENT
;
pro lmt_event, ev

COMPILE_OPT IDL2, HIDDEN

widget_control, ev.top, get_uvalue = s

case tag_names(ev, /structure_name) of
   'WIDGET_BUTTON': begin
      widget_control, ev.id, get_uvalue = uval
      case uval of
         'QUIT': begin
            widget_control, ev.top, /destroy
            return
         end

         else: help, ev
      endcase
   end

   'WIDGET_BASE': begin         ; handle base resize events
      widget_control, ev.id, tlb_get_size = newsize
      wdraw = widget_info(ev.top, find_by_uname = 'IMAGEDRAW')
      xy = newsize - (*s).w.pad
      widget_control, wdraw, $
                      draw_xsize = xy[0], draw_ysize = xy[1], $
                      scr_xsize = xy[0], scr_ysize = xy[1]
      wdraw = widget_info(ev.top, find_by_uname = 'PROFILEDRAW')
      xy = newsize - (*s).w.pad
      widget_control, wdraw, $
                      draw_xsize = xy[0], draw_ysize = xy[1], $
                      scr_xsize = xy[0], scr_ysize = xy[1]
      wdraw = widget_info(ev.top, find_by_uname = 'FITDRAW')
      xy = newsize - (*s).w.pad
      widget_control, wdraw, $
                      draw_xsize = xy[0], draw_ysize = xy[1], $
                      scr_xsize = xy[0], scr_ysize = xy[1]
   end

   'WIDGET_TAB': begin
      case ev.tab of 
         2: lmt_update_fit_plot, s
         3: lmt_update_res_plot, s
         else:
      endcase
   end
      
   else: help, ev                      ; do nothing
endcase

end

;;;;;
;
; LMT_WIDGETS()
;
; Create widget hierarchy
;
function lmt_widgets, p

COMPILE_OPT IDL2, HIDDEN

f = p.h -> get(position = p.n)  ; selected feature

base = widget_base(title = 'Lorenz-Mie Analysis Tool', $
                   mbar = mbar, /column, /tlb_size_events)

file_menu = widget_button(mbar, value = 'File', /menu)
void = widget_button(file_menu, value = '&Save Values...', $
                     event_pro = 'lmt_savevalues')
void = widget_button(file_menu, value = 'Save &Routine...', $
                     event_pro = 'lmt_saveroutine')
void = widget_button(file_menu, value = '&Quit', $
                     uvalue = 'QUIT')

; tabbed windows
wtab = widget_tab(base)

; hologram tab: show image
wimagetab = widget_base(wtab, title = 'Hologram', /column)
wimagedraw = widget_window(wimagetab, uname = 'IMAGEDRAW', $
                           mouse_down_handler = 'lmt_select_feature')
void = widget_label(wimagetab, value = 'Imaging Parameters', /align_left)
wimaging = widget_base(wimagetab, row = 2, /frame, /grid_layout, $
                       event_pro = 'lmt_update_parameter')
void = cw_field(wimaging, title = 'lambda [um]', uvalue = 'LAMBDA', $
                value = p.h.lambda, $
                /column, /return_events, /floating, /tab_mode)
void = cw_field(wimaging, title = 'mpp [um/pixel]', uvalue = 'MPP', $
                value = p.h.mpp, $
                /column, /return_events, /floating, /tab_mode)
void = cw_bgroup(wimaging, ['none', 'odd', 'even'], label_top = 'deinterlace', $
                 /frame, /row, /exclusive, set_value = p.h.deinterlace, $
                 uvalue = 'DEINTERLACE', /no_release)
void = cw_field(wimaging, title = 'nm', uvalue = 'NM', value = real_part(p.h.nm), $
                /column, /return_events, /floating, /tab_mode)
void = cw_field(wimaging, title = 'km', uvalue = 'KM', value = imaginary(p.h.nm), $
                /column, /return_events, /floating, /tab_mode)
void = cw_field(wimaging, title = 'nfringes', uvalue = 'NFRINGES', value = p.h.nfringes, $
                /column, /return_events, /floating, /tab_mode)

; profile tab
wprofiletab = widget_base(wtab, title = 'Profile', /column)
wprofiledraw = widget_window(wprofiletab, uname = 'PROFILEDRAW')

void = widget_label(wprofiletab, value = 'Particle Parameters', /align_left)
wprofile = widget_base(wprofiletab, row = 2, /frame, /grid_layout, $
                       event_pro = 'lmt_update_parameter')
wzp = cw_fslider(wprofile, title = 'zp [pixels]', uvalue = 'ZP', value = (f.rp)[2], $
                 min = 0., max = 300., /double, /edit, /drag, /tab_mode)
wnp = cw_fslider(wprofile, title = 'np', uvalue = 'NP', value = real_part(f.np), $
                 min = 1.33, max = 3.d, /double, /edit, /drag, /tab_mode)
walpha = cw_fslider(wprofile, title = 'alpha', uvalue = 'ALPHA', value = f.alpha, $
                    min = 0.4, max = 2., /double, /edit, /drag, /tab_mode)
wrad = cw_field(wprofile, title = 'rad [pixel]', uvalue = 'RAD', value = fix(f.rad), $
                /column, /return_events, /integer, /tab_mode)
wap = cw_fslider(wprofile, title = 'ap [um]', uvalue = 'AP', value = f.ap, $
                 min = 0.05d, max = 4.d, /double, /edit, /drag, /tab_mode)
wkp = cw_fslider(wprofile, title = 'kp', uvalue = 'KP', value = imaginary(f.np), $
                 min = 0.d, max = 10.d, /double, /edit, /drag, /tab_mode)
wdelta = cw_fslider(wprofile, title = 'delta [lambda]', uvalue = 'DELTA', value = f.delta, $
                    min = -2., max = 2, /double, /edit, /drag, /tab_mode)
void = widget_button(wprofile, value = 'Estimate', uvalue = 'ESTIMATE', /no_release)

; fitting tab
wfittab = widget_base(wtab, title = 'LM Fit', /column)
wfitdraw = widget_window(wfittab, uname = 'FITDRAW')
void = widget_label(wfittab, value = 'Fitting Parameters', /align_left)
wfitting = widget_base(wfittab, row = 1, /frame, $
                       event_pro = 'lmt_update_parameter')
void = cw_bgroup(wfitting, ['ap', 'zp', 'np', 'kp', 'alpha', 'delta', 'nm', 'km'], $
                 label_top = 'fixed parameters', /frame, /row, /nonexclusive, $
                 uvalue = 'FIX', $
                 set_value = [p.fix.ap, p.fix.zp, p.fix.np, p.fix.kp, $
                              p.fix.alpha, p.fix.delta, p.fix.nm, p.fix.km])
void = widget_button(wfitting, value = '  Fit  ', event_pro = 'lmt_fit', /no_release)

; residual tab
wrestab = widget_base(wtab, title = 'Residuals', /column)
wresdraw = widget_window(wrestab, uname = 'RESDRAW')

;;; object references to draw objects are available only after
;;; the base widget is realized

widget_control, base, /realize

;;; Hologram
widget_control, wimagedraw, get_value = wimage
wimage.select
im = image(deinterlace(p.h.data, p.h.deinterlace), $
           axis_style = 2, margin = [0.05, 0.05, 0.01, 0.01], /current)
flist = list()
features = p.h -> get(/all)
foreach feature, features, ndx do $
   flist.add, lmtfeature(im, feature, ndx)
flist[p.n].update, /selected

;;; Profile
widget_control, wprofiledraw, get_value = wprofile
wprofile.select
xrange = [0, f.rad]
profile = f.profile
dprofile = f.dprofile
pplot = plot(profile, thick = 2, $
             xrange = xrange, /xstyle, $
             xtitle = 'r [pixel]', ytitle = 'B(r)', /current)
x = findgen(f.rad+1)
ppoly = polygon([x, reverse(x)], [profile+dprofile, reverse(profile-dprofile)], $
                linestyle = 2, /data, $
                /fill_background, fill_color = 'light yellow')
pfit = plot(f.lmprofile, color = 'red', thick = 2, /over)
baseline = plot(xrange, [1, 1], linestyle = 2, /over)
pplot.order, /send_to_back
baseline.order, /send_to_back
ppoly.order, /send_to_back

;;; Fitting
widget_control, wfitdraw, get_value = wfit
wfit.select
imhologram = image(f.data, layout = [2, 1, 1], $
                   axis_style = 2, margin = replicate(0.01, 4), /current)
ax = imhologram.axes
ax[1].hide = 1
imfit = image(f.lmhologram, layout = [2, 1, 2], $
              axis_style = 2, margin = replicate(0.01, 4), /current)
ax = imfit.axes
ax[1].hide = 1

;;; Residuals
widget_control, wresdraw, get_value = wres
wres.select
imres = image(f.data-f.lmhologram, axis_style = 2, /current)

;;; cache padding between the base and the draw
widget_control, base, tlb_get_size = basesize
pad = basesize[0:1] - [640, 512]

;;; widget hierarchy
w = {lmt_widgets, $
     pad: pad, $
     im: im, $
     imhologram: imhologram, $
     imfit: imfit, $
     imres: imres, $
     f: flist, $
     pplot: pplot, $
     ppoly: ppoly, $
     pfit: pfit, $
     zp: wzp, $
     ap: wap, $
     np: wnp, $
     kp: wkp, $
     alpha: walpha, $
     delta: wdelta, $
     rad: wrad $
}
 
;;; program state
s = {lmt_state, $
     w: w, $
     p: p $
    }
ps = ptr_new(s, /no_copy)
wimage.uvalue = ps
widget_control, base, set_uvalue = ps;, /no_copy

return, base
end

;;;;;
;
; LMT_CLEANUP
;
pro lmt_cleanup, base

COMPILE_OPT IDL2, HIDDEN

widget_control, base, get_uvalue = s
ptr_free, s
end

;;;;;
;
; LMTOOL
;
pro lmtool, a, lambda, mpp, nm = nm, deinterlace = deinterlace, smooth = smooth

COMPILE_OPT IDL2

umsg = 'USAGE: dgglmtool, a, lambda, mpp'

if n_params() ne 3 then begin
   message, umsg, /inf
   return
endif

h = dgglmhologram(a, lambda, mpp, nm = nm, $
                  deinterlace = deinterlace, smooth = smooth)
if ~isa(h, 'DGGlmHologram') then begin
   message, umsg, /inf
   message, 'could not initialize', /inf
   return
endif

if h.count() le 0 then begin
   message, 'no features found in hologram', /inf
   return
endif

fix = {lmt_fix, $
       ap: 0, $
       zp: 0, $
       np: 0, $
       kp: 1, $
       alpha: 0, $
       delta: 1, $
       nm: 1, $
       km: 1 $
      }

p = {lmt_parameters, $
     h: h, $                    ; hologram object
     f: list(), $               ; list of features
     n: 0, $                    ; index of selected feature
     fix: fix $
    }

xmanager, 'lmt', lmt_widgets(p), /no_block, cleanup = 'lmt_cleanup'

end
