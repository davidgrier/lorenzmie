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
; 08/04/2013 DGG updated generated code for VOB analysis.
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
printf, lun, ';;; Fully-qualified name of the VOB file (without ".VOB")'
printf, lun, 'run = VTS_01_1'
printf, lun, ';;; Required parameters:'
printf, lun, ';;; Vacuum wavelength of illumination [micrometers]'
printf, lun, 'lambda = '+strtrim(h.lambda)
printf, lun, ';;; Magnification [micrometers/pixel]'
printf, lun, 'mpp = '+strtrim(h.mpp)
printf, lun, ';;; Setregion of interest'
printf, lun, 'x0 = 10 & x1 = 645'
printf, lun, 'y0 =  0 & y1 = 479'
printf, lun, ';;; File names:'
printf, lun, 'filename = run+".VOB"'
printf, lun, 'background = run+".bg.gdf"'
printf, lun, 'filename = run+"_"'+dgtimestamp(/date)+'.sav'
printf, lun, ';;; Open video file:'
printf, lun, 'vob = DGGgrMPlayer(fn, /gray, dim = [656,480])'
printf, lun, ';;; Compute background estimate'
printf, lun, 'if file_test(background, /read) then     $'
printf, lun, '  bg = read_gdf(background)              $'
printf, lun, 'else begin                               $'
printf, lun, '  bg = bgestimate(vob)                 & $'
printf, lun, '  write_gdf, bg, bgname                & $'
printf, lun, '  vob.rewind                           & $'
printf, lun, 'end'
printf, lun, 'bg = bg[x0:x1, y0:y1]'
printf, lun, ';;; Extract features from VOB'
printf, lun, 'features = list()'
printf, lun, 'while ~vob.eof do begin                & $'
printf, lun, '  print, vob.framenumber               & $'
printf, lun, '  a = float(vob.next)/bg               & $'
printf, lun, '  p = lmfeature(a, lambda, mpp, /gpu)  & $'
printf, lun, '  features.add, p                      & $'
printf, lun, '  if (vob.framenumber mod 10) eq 0 then  $'
printf, lun, '    save, po, pe, filename = filename  & $'
printf, lun, 'endwhile'
printf, lun, ';;; Save features'
printf, lun, 'save, features, filename = filename, $'
printf, lun, '   description = "Features detected in "+fn+" on "+systime()'
printf, lun, 'end'
close, lun
free_lun, lun

end
