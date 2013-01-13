;;;;
;
; LMT_FIT_UPDATE
;
; Update values and plots in fitting tab
;
pro lmt_fit_update, s

COMPILE_OPT IDL2, HIDDEN

if (*s).p.count le 0 then return

; estimate radius for feature
rp = [(*s).p.xp, (*s).p.yp]
sz = size(*((*s).p.a), /dimensions)
hm = ct_hitormiss(*((*s).p.a), [(*s).p.xp, (*s), p.yp], (*s).p.noise, /coordinates)
rad = round(100.*sqrt(total(hm[2, *])/n_elements(hm[2, *])))
; crop image to region of interest
x0 = round((*s).p.xp - rad) > 0L
x1 = round((*s).p.xp + rad) < sz[0] - 1L
y0 = round((*s).p.yp - rad) > 0L
x1 = round((*s).p.yp + rad) < sz[1] - 1L
aa =  a[x0:x1, y0:y1]
r0 = double([x0, y0])
r1 = double([x1, y1])
rc = rp - r0
origin = (r0 + r1)/2.d
drc =  rp - origin

b = aziavg(aa, rc = rc, rad = rad)
(*s).p.ap -> putdata, findgen(rad), b

end

;;;;
;
; LMT_FIT_EVENT
;
; Handle events for the fitting tab
;
pro lmt_fit_event, ev

COMPILE_OPT IDL2, HIDDEN

widget_control, ev.id, get_uvalue = uval, get_value = val
widget_control, ev.top, get_uvalue = s
case uval of
   'XP': (*s).p.xp = val
   'YP': (*s).p.yp = val
   'ZP': (*s).p.zp = val
   'AP': (*s).p.ap = val
   'NP': (*s).p.np = dcomplex(val, imaginary((*s).p.np))
   'KP': (*s).p.np = dcomplex(real_part((*s).p.np), val)
endcase

end

;;;;
;
; LMT_FD_UPDATE
;
; Update values and plots affected by the feature detection tab
;
pro lmt_fd_update, s

COMPILE_OPT IDL2, HIDDEN

features = ctfeature(*((*s).p.a), range = (*s).p.range, noise = (*s).p.noise, $
                     smooth = (*s).p.smooth, threshold = (*s).p.threshold, $
                     count = count, ct = ct, /quiet)
(*s).w.fd -> putdata, bytscl(ct)
widget_control, (*s).w.fdmax, set_value = max(ct)
widget_control, (*s).w.fdcount, set_value = count
if count le 0 then begin
   features = [-1, -1]
   (*s).p.ndx = 0
endif

(*s).w.fp -> putdata, features[0, *], features[1, *]

end

;;;;
;
; LMT_FD_EVENT
;
; Handle events for the feature detection tab
;
pro lmt_fd_event, ev

COMPILE_OPT IDL2, HIDDEN

widget_control, ev.id, get_uvalue = uval, get_value = val
widget_control, ev.top, get_uvalue = s
case uval of
   'RANGE'     : (*s).p.range     = val
   'NOISE'     : (*s).p.noise     = val
   'SMOOTH'    : (*s).p.smooth    = val
   'THRESHOLD' : (*s).p.threshold = val
   'FEATURES'  : (*s).p.count     = val
endcase

lmt_fd_update, s

end

;;;;
;
; LMT_IMAGE_EVENT
;
; Handle events for the image tab
;
pro lmt_image_event, ev

COMPILE_OPT IDL2, HIDDEN

widget_control, ev.id, get_uvalue = uval, get_value = val
widget_control, ev.top, get_uvalue = s
case uval of
   'LAMBDA' : (*s).p.lambda = val
   'MPP'    : (*s).p.mpp = val
   'NM'     : (*s).p.nm = dcomplex(val, imaginary((*s).p.nm))
   'KM'     : (*s).p.nm = dcomplex(real_part((*s).p.nm), val)
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

base = widget_base(title = 'Lorenz-Mie Analysis Tool', mbar = mbar, /column, $
                   /tlb_size_events)

file_menu = widget_button(mbar, value = 'File', /menu)
void = widget_button(file_menu, value = '&Save Values...', $
                     event_pro = 'lmt_savevalues', accelerator = 'Alt+S')
void = widget_button(file_menu, value = 'Save &Routine...', $
                     event_pro = 'lmt_saveroutine', accelerator = 'Alt+R')
void = widget_button(file_menu, value = '&Quit', $
                     uvalue = 'QUIT', accelerator = 'Alt+Q')

; tabbed windows
wtab = widget_tab(base)

; hologram tab: show image
wimagetab = widget_base(wtab, title = 'Hologram', /column)
wimagedraw = widget_window(wimagetab, uvalue = 'imagedraw', uname = 'IMAGEDRAW')
void = widget_label(wimagetab, value = 'Imaging Parameters', /align_left)
wimaging = widget_base(wimagetab, row = 2, /frame, /grid_layout, event_pro = 'lmt_image_event')
void = cw_field(wimaging, title = 'lambda', uvalue = 'LAMBDA', value = p.lambda, $
                /return_events, /floating)
void = cw_field(wimaging, title = '   mpp', uvalue = 'MPP', value = p.mpp, $
                /return_events, /floating)
void = cw_field(wimaging, title = '    nm', uvalue = 'NM', value = real_part(p.nm), $
                /return_events, /floating)
void = cw_field(wimaging, title = '    km', uvalue = 'KM', value = imaginary(p.nm), $
                /return_events, /floating)

; feature detection tab
wfeaturetab =  widget_base(wtab,  title = 'Feature Detection', /column)
wfeaturedraw = widget_window(wfeaturetab, uvalue = 'featuredraw', uname = 'FEATUREDRAW')
void = widget_label(wfeaturetab, value = 'Feature Detection Parameters', /align_left)
wfeature = widget_base(wfeaturetab, row = 2, /frame, /grid_layout, event_pro = 'lmt_fd_event')
void    = cw_field(wfeature, title = '    range', uvalue = 'RANGE', value = p.range, $
                   /return_events, /integer)
void    = cw_field(wfeature, title = '    noise', uvalue = 'NOISE', value = p.noise, $
                   /return_events, /floating)
void    = cw_field(wfeature, title = '   smooth', uvalue = 'SMOOTH', value = p.smooth, $
                   /return_events, /integer)
void    = cw_field(wfeature, title = 'threshold', uvalue = 'THRESHOLD', value = p.threshold, $
                   /return_events, /integer)
fdmax   = cw_field(wfeature, title = '      max', value = 0, /integer)
fdcount = cw_field(wfeature, title = 'nfeatures', value = 0, /integer)

; fitting tab: fit to a specific feature
wfittab = widget_base(wtab, title = 'Fitting', /column)
wfitdraw = widget_window(wfittab, uvalue = 'fitdraw', uname = 'FITDRAW')
void = widget_label(wfittab, value = 'Fitting Parameters', /align_left)
wfitting = widget_base(wfittab, row = 2, /frame, /grid_layout, event_pro = 'lmt_fit_event')
void = cw_field(wfitting, title = 'xp', uvalue = 'XP', value = p.xp, $
                /return_events, /floating)
void = cw_field(wfitting, title = 'yp', uvalue = 'YP', value = p.yp, $
                /return_events, /floating)
void = cw_fslider(wfitting, title = 'zp', uvalue = 'ZP', value = p.zp, $
                  min = 0. , max = 200. , /double, /drag)
void = cw_fslider(wfitting, title = 'ap', uvalue = 'AP', value = p.ap, $
                  min = 0.05d, max = 4.d, /double, /drag)
void = cw_fslider(wfitting, title = 'np', uvalue = 'NP', value = real_part(p.np), $
                  min = 1.001d * real_part(p.nm), max = 3.d, /double, /drag)
void = cw_fslider(wfitting, title = 'kp', uvalue = 'KP', value = imaginary(p.np), $
                  min = 0.d, max = 10.d, /double, /drag)

widget_control, base, /realize

;;; object references to draw objects are only available 
;;; after the widget is realized
;;; Image
widget_control, wimagedraw, get_value = wimage
wimage.select
im = image(*(p.a), axis_style = 2, margin = [0.05, 0.05, 0.01, 0.01], /current)
fp = plot([-1], [-1], symbol = 'o', color = 'red', linestyle = 6, /over)

;;; Feature Detection
widget_control, wfeaturedraw, get_value = wfeature
wfeature.select
fd = image(*(p.a), axis_style = 2, margin = [0.05, 0.05, 0.01, 0.01], /current)

;;; Fitting
widget_control, wfitdraw, get_value = wfit
wfit.select
ap = plot([0, 100], [1, 1], xtitle = 'r [pixel]', ytitle = 'B(r)', /current)

;;; cache padding between the base and the draw
widget_control, base, tlb_get_size = basesize
pad = basesize[0:1] - [640, 512]

;;; widget hierarchy
w = {lmt_widgets, $
     pad: pad, $
     im: im, $
     fp: fp, $
     fd: fd, $
     fdmax: fdmax, $
     fdcount: fdcount, $
     ap: ap $
    }

;;; program state
s = {lmt_state, $
     w: w,      $
     p: p       $
    }
ps = ptr_new(s, /no_copy)
lmt_fd_update, ps
widget_control, base, set_uvalue = ps, /no_copy
return, base
end

;;;;;
;
; LMT_WIDGETVALUE()
;
; Return value of widget
;
function lmt_widgetvalue, wid

COMPILE_OPT IDL2, HIDDEN

widget_control, wid, get_value = val
return, val
end

;;;;;
;
; LMT_EVENT
;
; Event loop handler
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

   'WIDGET_BASE': begin
      ;; Handle base resize events
      widget_control, ev.id, tlb_get_size = newsize
      wdraw = widget_info(ev.top, find_by_uname = 'IMAGEDRAW')
      xy = newsize - (*s).w.pad
      widget_control, wdraw, $
                      draw_xsize = xy[0], draw_ysize = xy[1], $
                      scr_xsize = xy[0], scr_ysize = xy[1]
      wdraw = widget_info(ev.top, find_by_uname = 'FEATUREDRAW')
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

   else:                        ; do nothing
endcase

end

;;;;;
;
; LMT_CLEANUP
;
; Free resources
;
pro lmt_cleanup, base

COMPILE_OPT IDL2, HIDDEN

widget_control, base, get_uvalue = s
ptr_free = s
end

;;;;;
;
; LMTOOL
;
; Main routine
;
pro lmtool, a, lambda, mpp, nm = nm, dimensions = dimensions

COMPILE_OPT IDL2

umsg = 'USAGE: lmtool, a, lambda, mpp'

if n_params() ne 3 then begin
   message, umsg, /inf
   return
endif

if ~isa(a, /number, /array) then begin
   message, umsg, /inf
   message, 'A must be a numerical array', /inf
   return
endif

if size(a, /n_dimensions) ne 2 then begin
   message, umsg, /inf
   message, 'A must be a two dimensional array', /inf
   return
endif

if round(mean(a)) ne 1 then begin
   message, umsg, /inf
   message, 'A must be a normalized hologram', /inf
   return
endif

if ~isa(lambda, /number, /scalar) then begin
   message, umsg, /inf
   message, 'LAMBDA must be a numerical scalar', /inf
   return
endif

if ~isa(mpp, /number, /scalar) then begin
   message, umsg, /inf
   message, 'MPP must be a numerical scalar', /inf
   return
endif

if ~isa(nm, /number, /scalar) then $
   nm = refractiveindex(lambda, 24.)

if n_elements(dimensions) ne 2 then $
   dimensions = [640, 480]

features = ctfeature(a, range = range, noise = noise, $
                     smooth = smooth, threshold = threshold, $
                     count = count, /quiet)

pa = ptr_new(a)

p = {lmt_parameters, $
     a: pa, $                     ; normalized hologram
     lambda: double(lambda), $    ; wavelength [micrometers]
     mpp: double(mpp), $          ; micrometers per pixel
     nm: dcomplex(nm), $          ; refractive index of medium
     xp: 0.d, $                   ; coordinate of particle [pixel]
     yp: 0.d, $                   ;
     zp: 0.d, $                   ;
     ap: 0.d, $                   ; radius of particle [micrometers]
     np: dcomplex(0.), $          ; refractive index of particle
     range: fix(range), $         ; range for circletransform
     noise: float(noise), $       ; noise level for feature identification
     smooth: fix(smooth), $       ; smoothing factor for circletransform
     threshold: fix(threshold), $ ; count threshold for feature identification
     count: fix(count), $         ; specified number of features to find
     ndx: 0 $                     ; index of the feature to fit
    }

xmanager, 'lmt', lmt_widgets(p), /no_block, cleanup = 'lmt_cleanup'

end
