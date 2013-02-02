;;;;;
;
; LMT_SAVEROUTINE
;
; Component of LMTool
;
; Create a batch file that processes holographic movies.
;
; MODIFICATION HISTORY
; 02/02/2013 Written by David G. Grier, New York University
;
; Copyright (c) 2013 David G. Grier
;
pro lmt_saveroutine, ev

COMPILE_OPT IDL2, HIDDEN

fn = dialog_pickfile(/write, /overwrite_prompt, $
                     title = 'Save Routine', filter = '*.pro', $
                     default_extension = 'pro')

if strlen(fn) lt 1 then return

widget_control, ev.top, get_uvalue = s
h = (*s).p.h

openw, lun, fn, /get_lun
printf, lun, ';;; routine generated automatically by LMTOOL'
printf, lun, ';;; '+systime()
printf, lun, ';;; Specify the fully-qualified name of the VOB file'
printf, lun, 'fn = XXX'
printf, lun, ';;; Specify the dimensions of a frame'
printf, lun, 'dim = [640, 480]'
printf, lun, 'vob = DGGgrMPlayer(fn, /gray, dim = dim)'
printf, lun, ';;; Compute background estimate'
printf, lun, 'bg = bgestimate(vob) > 1.'
printf, lun, 'vob.rewind'
printf, lun, ';;; Vacuum wavelength of illumination [micrometers]'
printf, lun, 'lambda = '+strtrim(h.lambda)
printf, lun, ';;; Magnification [micrometers/pixel]'
printf, lun, 'mpp = '+strtrim(h.mpp)
printf, lun, ';;; Extract features from VOB'
printf, lun, 'features = list()'
printf, lun, 'while ~vob.eof do begin          & $'
printf, lun, '   a = float(vob.next)/bg        & $'
printf, lun, '   p = lmfeature(a, lambda, mpp) & $'
printf, lun, '   features.add, p               & $'
printf, lun, 'endwhile'
printf, lun, ';;; Save features'
printf, lun, ';;; file name:
printf, lun, 'filename = XXX'
printf, lun, 'save, features, filename = filename, $'
printf, lun, '   description = "Features detected in "+fn+" on "+systime()'
printf, lun, 'end'
close, lun
free_lun, lun

end
