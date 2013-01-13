function dhmbackground, a, input_bg, $
                        lambda = lambda, $
                        mpp = mpp, $
                        threshold = threshold, $
                        pickn = pickn, $
                        beta = beta, $
                        graphics = graphics, $
                        gpu = gpu, $
                        quiet = quiet, $
                        foreground = foreground

COMPILE_OPT IDL2
      
umsg = 'USAGE: bg = dhmbackground(a, [input_bg])'
                  
gpu = keyword_set(gpu)
quiet = keyword_set(quiet)

if ~isa(beta, /number, /scalar) then $
    beta = 0.5

convergence = 1

nbg = median(a)
ndelta = max(a)
if n_params() eq 2 then begin
   if ~array_equal(size(a,/dimensions), size(input_bg, /dimensions)) then begin
      message, umsg, /inf
      message, 'a and input_bg must have the same dimensions', /inf
      return, nbg
   endif
   nbg = float(input_bg)
   ndelta = mean(nbg)
endif

repeat begin
   bg = nbg 
   delta = ndelta
   b = float(a)/bg
   p = dhmfeature(b, lambda = lambda, mpp = mpp, $
                  threshold = threshold, $
                  pickn = pickn, $
                  graphics = graphics, $
                  gpu = gpu, $
                  quiet = quiet)
   fg = feature2dhm(p, size(a, /dimensions), lambda = lambda, mpp = mpp, gpu = gpu)
   nbg = (1.-beta)*bg + beta*float(a)/fg
   ndelta = max(abs(bg-nbg))
endrep until delta - ndelta lt convergence

if arg_present(foreground) then $
   foreground = fg

return, nbg
end
