;;;;;
;
; LMT_SAVEVALUES
;
; Component of LMTool
;
; Save values describing the currently selected feature to a file
;
; MODIFICATION HISTORY
; 02/01/2013 Written by David G. Grier, New York University
;
; Copyright (c) 2013 David G. Grier
;
pro lmt_savevalues, ev

COMPILE_OPT IDL2, HIDDEN

fn = dialog_pickfile(/write, /overwrite_prompt, $
                     title = 'Save Values', filter = '*.gdf', $
                     default_extension = 'gdf')

if strlen(fn) lt 1 then return

widget_control, ev.top, get_uvalue = s
h = (*s).p.h
f = h.get(position = (*s).p.n)
write_gdf, [f.rp, f.ap, $
            real_part(f.np), imaginary(f.np), $
            real_part(h.nm), imaginary(h.nm), $
            f.alpha, f.delta], fn, /ASCII

end
