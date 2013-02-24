;+
; NAME:
;    lmf_report
;
; PURPOSE:
;    Report parameters obtained with lmfeature
;
; CATEGORY:
;    Holographic video microscopy
;
; MODIFICATION HISTORY:
; 02/24/2013 Written by David G. Grier, New York University
;
; Copyright (c) 2013 David G. Grier
;-

pro lmf_report, msg, p

COMPILE_OPT IDL2, HIDDEN

message, msg, /inf
message, string(p[0:2], $
                format = '("  rp = (",F0.2,", ",F0.2,", ",F0.2,")")'), /inf
message, string(p[3:4], $
                format = '("  ap = ",F0.3," um, np = ",F0.3)'), /inf
message, string(p[8:9], $
                format = '("  alpha = ",F0.3,", delta = ",F0.3)'), /inf
message, /inf

return
end
