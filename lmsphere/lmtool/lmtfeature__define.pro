;+
; NAME:
;    LMTFeature
;
; PURPOSE:
;    Object describing the graphical representation of a feature in
;    the interactive Lorenz-Mie hologram analysis tool.
;
; CATEGORY:
;    Holographic video microscopy
;
; CALLING SEQUENCE:
;    f = LMTFeature(parent, feature, index)
;
; INPUTS:
;    parent: Graphic object into which the feature is to be introduced
;
;    feature: DGGlmFeature object describing the feature
;
;    index: index of the feature in the list of features represented
;        by the parent.
;
; METHODS:
;    LMTFeature::Update [, /selected]
;        Update the graphical representation of the feature.  If the
;        SELECTED flag is set, graphical representation reflects this
;        status.
;
; MODIFICATION HISTORY:
; 01/27/2013 Written by David G. Grier, New York University
;
; Copyright (c) 2013 David G. Grier
;-

;;;;;
;
; LMTFEATURE::Update
;
pro LMTFeature::Update, selected = selected

feature = self.feature
r0 = feature.r0
rad = feature.rad
r1 = r0 + 2*rad + 1
self.poly -> setdata, [r0[0], r1[0], r1[0], r0[0], r0[0]], $
                      [r0[1], r0[1], r1[1], r1[1], r0[1]]
self.poly.color =  keyword_set(selected) ? 'light green' : 'green'

self.point -> setdata, [(feature.rp)[0]], [(feature.rp)[1]]
end

;;;;;
;
; LMTFEATURE::Init()
;
function LMTFeature::Init, parent, feature, ndx

COMPILE_OPT IDL2, HIDDEN

umsg = 'USAGE: LMTFeature(parent, feature, index)'
if n_params() ne 3 then begin
   message, umsg, /inf
   return, 0
endif

if ~isa(parent, 'graphic') then begin
   message, umsg, /inf
   message, 'PARENT must be a graphic object', /inf
   return, 0
endif

if ~isa(feature, 'DGGlmFeature') then begin
   message, umsg, /inf
   message, 'FEATURE must be a DGGlmFeature', /inf
   return, 0
endif

if ~isa(ndx, /scalar, /number) then begin
   message, umsg, /inf
   message, 'INDEX should be the integer index of the FEATURE', /inf
   return, 0
endif

self.feature = feature

r0 = feature.r0
rad = feature.rad
r1 = r0 + 2*rad + 1
self.poly = polygon([[r0[0], r1[0], r1[0], r0[0], r0[0]], $
                     [r0[1], r0[1], r1[1], r1[1], r0[1]]], $
                    target = parent, /data, fill_transparency = 90, $
                    linestyle = 2, color = 'green', thick = 2, uvalue = ndx)
self.point = plot([(feature.rp)[0]], [(feature.rp)[1]], $
                  color = 'red', thick = 2, symbol = 'o', $
                  overplot = parent)

return, 1
end

;;;;;
;
; LMTFEATURE::Cleanup
;
pro LMTFeature::Cleanup

COMPILE_OPT IDL2, HIDDEN

end

;;;;;
;
; LMTFEATURE__DEFINE
;
pro lmtfeature__define

COMPILE_OPT IDL2

struct = {LMTFeature, $
          inherits IDL_Object, $
          feature: obj_new(), $
          point: obj_new(), $
          poly: obj_new() $
         }
end
