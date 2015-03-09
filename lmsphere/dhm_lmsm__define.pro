;;;;;
;
; DHM_LMsm::ComputeGeometry
;
pro DHM_LMsm::ComputeGeometryCPU

  COMPILE_OPT IDL2, HIDDEN

  self.DGGdhmLMsphere::ComputeGeometryCPU

  mask = *self.mask
  v = {x:   self.geometry.x[mask], $
       y:   self.geometry.y[mask], $
       rho: self.geometry.rho[mask], $
       kr:  self.geometry.kr[mask], $
       costheta: self.geometry.costheta[mask], $
       sintheta: self.geometry.sintheta[mask], $
       cosphi: self.geometry.cosphi[mask], $
       sinphi: self.geometry.sinphi[mask], $
       coskr: self.geometry.coskr[mask], $
       sinkr: self.geometry.sinkr[mask] $
      }
  self.mgeometry = ptr_new(v, /no_copy)
end
  
;;;;;
;
; DHM_LMsm::Init()
;
function DHM_LMsm::Init, mask, $
                         _ref_extra = re

  COMPILE_OPT IDL2, HIDDEN

  if ~isa(mask, /number, /array) then begin
     message, 'USAGE: DHM_LMsm(mask)', /inf
     return, 0B
  endif
  
  if ~self.DGGdhmLMsphere::Init, _extra = re then $
     return, 0B

  self.mask = ptr_new(mask)
  self.ComputeGeometryCPU
  self.v = self.mgeometry
  
end

;;;;;
;
; DHM_LMsm::Cleanup
;
pro DHM_LMsm::Cleanup

  COMPILE_OPT IDL2, HIDDEN

  self.DGGdhmLMsphere::Cleanup

  if ptr_valid(mask) then $
     ptr_free, mask

  if ptr_valid(mgeometry) then $
     ptr_free, mgeometry
end

;;;;;
;
; DHM_LMsm__define
;
pro dhm_lmsm__define

  COMPILE_OPT IDL2, HIDDEN

  struct = {DHM_LMsm, $
            inherits DGGdhmLMsphere, $
            mask: ptr_new(),  $
            mgeometry: ptr_new() $
           }
end

