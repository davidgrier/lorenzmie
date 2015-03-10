;;;;;
;
; hvmLorenzMie::GetProperty
;
pro hvmLorenzMie::GetProperty, hologram = hologram, $
                               fraction = fraction, $
                               mask = mask, $
                               _ref_extra = re

  COMPILE_OPT IDL2, HIDDEN

  if arg_present(hologram) then $
     self.LorenzMie::GetProperty, hologram = hologram

  if arg_present(fraction) then $
     fraction = self.fraction

  if arg_present(mask) then $
     mask = *self.mask

  self.dhmLorenzMie::GetProperty, _extra = re
end

;;;;;
;
; hvmLorenzMie::MakeCoordinates()
;
function hvmLorenzMie::MakeCoordinates, dimensions, r0

  COMPILE_OPT IDL2, HIDDEN

  grid = self.dhmLorenzMie::MakeCoordinates(dimensions, r0)
  
  ;;; generate random mask
  npts = n_elements(grid['x'])
  mask = long(npts * randomu(seed, round(self.fraction*npts)))
  mask = mask[sort(mask)]
  u = uniq(mask)
  mask = mask[u]

  ;;; save original coordinates
  self.grid = grid
  
  ;;; use masked coordinates
  c = hash('x', (grid['x'])[mask], $
           'y', (grid['y'])[mask], $
           'z', (grid['z'])[mask])

  ;;; save mask
  self.mask = ptr_new(mask, /no_copy)

  return, c
end

;;;;;
;
; hvmLorenzMie::Init()
;
function hvmLorenzMie::Init, deinterlace = deinterlace, $
                             fraction = fraction, $
                             _ref_extra = re

  COMPILE_OPT IDL2, HIDDEN

  self.deinterlace = keyword_set(deinterlace)
  
  self.fraction = isa(fraction, /number, /scalar) ? float(fraction) > 0 < 1 : 0.5

  return, self.dhmLorenzMie::Init(_extra = re)
end
;;;;;
;
; hvmLorenzMie__define
;
pro hvmLorenzMie__define

  COMPILE_OPT IDL2, HIDDEN

  struct = {hvmLorenzMie, $
            inherits dhmLorenzMie, $
            deinterlace: 0, $
            fraction: 0., $
            grid: obj_new(), $
            mask: ptr_new() $
           }
end
