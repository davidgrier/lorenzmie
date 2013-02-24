;+
; NAME:
;    lmf_range
;
; PURPOSE:
;    Compute range from a specified feature over which to fit using lmfeature
;
; CATEGORY:
;    Holographic video microscopy
;
; CALLING SEQUENCE:
;    range = lmf_range(hologram, rc)
;
; INPUTS:
;    HOLOGRAM: normalized hologram
;
;    RC: (x,y) coordinate of feature
;
; KEYWORD PARAMETERS:
;    deinterlace: analyze odd field if set to an odd number,
;        even field if set to an even number.
;        Default: Analyze entire image
;
; OUTPUTS:
;    range: range in pixels over which to analyze hologram.
;
; PROCEDURE:
;    calls AZISTD to compute radial intensity profile around specificed
;    center and its standard deviation.  Range is set when signal is
;    weaker than uncertainty.
;
; MODIFICATION HISTORY:
; 02/24/2013 Written by David G. Grier, New York University
;
; Copyright (c) 2013 David G. Grier
;-

function lmf_range, a, rc, deinterlace = deinterlace

COMPILE_OPT IDL2

sigma = azistd(a, avg, center = rc, deinterlace = deinterlace)
w = where(abs(avg - 1.) ge sigma, ngood)

return, (ngood gt 0) ? max(w) > 30 : 30
end
