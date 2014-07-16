;+
; NAME:
;    SVRfeature
;
; PURPOSE:
;    Lorenz-Mie analysis of in-line holograms of colloidal spheres
;    using Support Vector Regression
;
; CATEGORY:
;    Holographic video microscopy
;
; EXAMPLE:
; Train SVRfeature to analyze holograms:
;    svr = svrfeature()
;    svr.lambda = 0.447
;    svr.mpp = 0.1351
;    svr.nm = refractiveindex(svr.lambda, 24)
;    svr.zrange = [50., 300.]
;    svr.arange = [0.4, 1.5]
;    svr.nrange = [1.36, 1.8]
;    svr.rad = 100
;    svr.prefix = 'svr_HCS_n'
;    svr.train
;    svr.save
;
; Use trained SVR to analyse the even field of an interlaced hologram:
;    svr = svrfeature(prefix = 'svr_HCS_n')
;    p = svr.feature(a, deinterlace=2)
;
; INHERITS:
;    IDL_Object
;
; PROPERTIES:
;    NOTE:
;        I: Initialization
;        G: Get
;        S: Set
;
;    rad     [ GS] radius for azimuthal average
;    lambda  [ GS] vacuum wavelength of illumination [micrometers]
;    mpp     [ GS] magnification [micrometers/pixel]
;    nm      [ GS] refractive index of medium
;    zrange  [ GS] range of axial positions [pixels]
;    arange  [ GS] range of sphere radius [micrometers]
;    nrange  [ GS] range of refractive index values
;    asvr    [ GS] python SVR object for computing radius of sphere
;    nsvr    [ GS] python SVR object for computing refractive index
;    zsvr    [ GS] python SVR object for computing axial position
;    path    [IGS] directory containing data files
;    prefix  [IGS] filename for SVR data
;
; METHODS:
;    p = SVRFeature::Feature(a)
;    DESCRIPTION:
;        Analyzes in-line hologram using support vector regression.
;    INPUTS:
;        a: 2-dimensional normalized hologram
;    KEYWORD PARAMETERS:
;        pickn: Maximum number of features to seek.
;            Default: All, as determined by CTFEATURE
;        count: On output, number of features actually found
;        deinterlace:
;            0 or absent: Analyze entire hologram
;            even value: Analyze even field
;            odd value : Analyze odd field
;    OUTPUT:
;        p: [5,count] array of feature parameters.
;            p[0,*]: x coordinate [pixels]
;            p[1,*]: y coordinate [pixels]
;            p[2,*]: z coordinate [pixels]
;            p[3,*]: sphere radius [micrometers]
;            p[4,*]: sphere refractive index
;    PROCEDURE:
;        Uses CTFEATURE to identify candidate features in a hologram.
;        Uses AZIMEDIAN to obtain the radial profile of each candidate
;            feature.
;        Applies support vector regression on each profile to obtain
;            estimates for the sphere's axial position, radius, and
;            refractive index.  SVRfeature must be trained with 
;            SVRfeature::Train using appropriate instrumental parameters
;            for Feature() to yield meaningful results.
;
;    SVRfeature::Train
;    DESCRIPTION:
;        Trains Support Vector Machines to recognize and analyze 
;        holograms of spheres created with in-line holographic video
;        microscopy.
;    PROCEDURE:
;        Calls LMSPHEREPROFILE to generate the radial profiles
;            of spheres' holograms over the range of parameters.
;        Uses the SVR class from scikit learn to train support vector
;            machines to fit to the training data, and then to attempt
;            to predict independently computed cross-correlation data.
;            Training proceeds by searching for optimal values of
;            gamma and C for the training set until the R^2 of the
;            cross-correlation data reaches its goal for each
;            parameter.  Training terminates when the goal is met,
;            or when limits are reached in the size of the training
;            set or on the time required for training.
;    NOTE:
;        Training can be _very_ time-consuming!
;
;    SVRfeature::Validate()
;    DESCRIPTION:
;        Computes the R^2 statistic for simulated data over the
;        range of values in the SVR training set.
;    KEYWORD PARAMETERS:
;        nsets: Number of data sets to use for validation.
;            Default: 1000L
;    OUTPUT:
;        rsq: [rsq_z, rsq_a, rsq_n]: R^2 statistics for
;            axial position (z), radius (a) and refractive index (n)
;
;    SVRfeature::Save
;    DESCRIPTION:
;        Saves current SVR system as an HDF5 file
;    KEYWORD PARAMETERS:
;        prefix: Name under which to save SVR system.
;            Default: Present prefix
;    SIDE EFFECTS:
;        Writes HDF5 file.
;    
;    result = SVRfeature::Load()
;    DESCRIPTION:
;        Reads SVR data from specified HDF5 file
;    KEYWORD PARAMETERS:
;        prefix: Name under which SVR system has been saved.
;            Default: Present prefix
;        filename: Complete filename for HDF5 file to be loaded.
;    OUTPUT:
;        result: 1 for success; 0 for failure
;
; REFERENCES:
; 1. A. Yevick, M. Hannel and D. G. Grier, "Machine-learning approach
;    to holographic particle characterization," Opt. Express,
;    submitted for publication (2014).
;
; 2. F. Pedregosa, et al., "Scikit-learn: Machine learning in Python,"
;    J. Mach. Learn. Res. 12, 2825-2830 (2011).
;
; 3. C.-C. Chang and C.-J. Lin, "LIBSVM: A library for Support Vector
;    Machines," ACM Trans. Intel. Sys. Techn. 2, 27 (2011).
;
; 4. R. Kling, "Using Python from IDL," (Kling Software, 2014)
; 
; DEPENDENCIES:
; Requires a working installation of Python, including scikit-learn.
; Requires a working installation of Slither from Jacquette
;    Consulting: http://www.slither4idl.com/
;
; MODIFICATION HISTORY:
; 07/08/2014 Written by David G. Grier, New York University
; 07/13/2014 DGG Scale target values into range [0,1] before fitting.
;
; Copyright (c) 2014 David G. Grier
;-
;;;;;
;
; SVRfeature::Feature()
;
function SVRfeature::Feature, a, $
                              pickn = pickn, $
                              count = count, $
                              deinterlace = deinterlace

COMPILE_OPT IDL2, HIDDEN

if ~self.trained then begin
   message, 'SVR system has not been trained', /inf
   return, 0
endif

rc = ctfeature(a, pickn = pickn, count = count)
if (count lt 1) then $
   return, rc

res = fltarr(5, count)
res[0:1, *] = rc[0:1, *]

zsvr = self.data['zsvr']
asvr = self.data['asvr']
nsvr = self.data['nsvr']
for i = 0, count-1 do begin
   b = azimedian(a, center = rc[0:1, i], rad = self.rad, deinterlace = deinterlace)
   res[2:4, i] = [zsvr.predict(b), asvr.predict(b), nsvr.predict(b)]
endfor
res[2, *] = res[2, *] * (self.zrange[1] - self.zrange[0]) + self.zrange[0]
res[3, *] = res[3, *] * (self.arange[1] - self.arange[0]) + self.arange[0]
res[4, *] = res[4, *] * (self.nrange[1] - self.nrange[0]) + self.nrange[0]

return, res
   
end

;;;;;
;
; SVRfeature::CreateTrainingData
;
pro svrfeature::GenerateData, zrange = zrange, $
                              arange = arange, $
                              nrange = nrange, $
                              rad = rad, $
                              nsets = nsets, $
                              quiet = quiet

COMPILE_OPT IDL2, HIDDEN

if isa(zrange, /number, /array) then begin
   self.zrange = float(zrange[0:1])
   zp = randomu(seed, nsets) * (self.zrange[1] - self.zrange[0]) + self.zrange[0]
endif else begin
   message, 'zrange must be specified', /inf
   return
endelse

if isa(arange, /number, /array) then begin
   self.arange = float(arange[0:1])
   ap = randomu(seed, nsets) * (self.arange[1] - self.arange[0]) + self.arange[0]
endif else begin
   message, 'arange must be specified', /inf
   return
endelse

if isa(nrange, /number, /array) then begin
   self.nrange = float(nrange[0:1])
   np = randomu(seed, nsets) * (self.nrange[1] - self.nrange[0]) + self.nrange[0]
endif else begin
   message, 'nrange must be specified', /inf
   return
endelse

if isa(rad, /number, /scalar) then $
   self.rad = rad

nsets = isa(nsets, /number, /scalar) ? long(nsets) : 10000L

noprint = keyword_set(quiet)
message, 'Computing profiles ...', /inf, noprint = noprint
if ~noprint then tic
profiles = fltarr(self.rad+1, nsets, /nozero)
for i = 0, nsets-1 do $
   profiles[*, i] = lmsphereprofile(findgen(self.rad+1), zp[i], ap[i], np[i], $
                                    self.nm, 1., 0.,  self.lambda, self.mpp)
if ~noprint then toc

self.data['zp'] = (zp - self.zrange[0])/(self.zrange[1] - self.zrange[0])
self.data['ap'] = (ap - self.arange[0])/(self.arange[1] - self.arange[0])
self.data['np'] = (np - self.nrange[0])/(self.nrange[1] - self.nrange[0])
self.data['profiles'] = profiles
end

;;;;;
;
; SVRfeature::Validate()
;
function svrfeature::Validate, nsets = nsets

COMPILE_OPT IDL2, HIDDEN

if ~isa(nsets, /scalar, /number) then nsets = 1000L

self -> GenerateData, zrange = self.zrange, $
                      arange = self.arange, $
                      nrange = self.nrange, $
                      rad = self.rad, $
                      nsets = nsets

profiles = self.data['profiles']
zsvr = self.data['zsvr']
asvr = self.data['asvr']
nsvr = self.data['nsvr']
zscore = zsvr.score(profiles, self.data['zp'])
ascore = asvr.score(profiles, self.data['ap'])
nscore = nsvr.score(profiles, self.data['np'])

return, [zscore, ascore, nscore]
end

;;;;;
;
; SVRfeature::Validate
;
pro svrfeature::Validate, nsets = nsets

COMPILE_OPT IDL2, HIDDEN

print, self.Validate(nsets = nsets)
end

;;;;;
;
; SVRfeature::TrainProperty, name
;
pro svrfeature::TrainProperty, name, $
                               goal = goal, $
                               quiet = quiet

COMPILE_OPT IDL2, HIDDEN

message, 'Training '+name+' ...', /inf

nsets = 2000L

if ~isa(goal, /number, /scalar) then $
   goal = 0.95

svr = self.data[name+'svr']
svr.tol = 1e-6
svr.epsilon = 0.001

self -> GenerateData, zrange = self.zrange, $
                      arange = self.arange, $
                      nrange = self.nrange, $
                      rad = self.rad, $
                      nsets = nsets, $
                      quiet = quiet
profiles = self.data['profiles']
tprofiles = profiles[*, 0:*:2]
vprofiles = profiles[*, 1:*:2]
prop = self.data[name+'p']
tprop = prop[0:*:2]
vprop = prop[1:*:2]

message, 'optimizing gamma ...', /inf
svr.C = 1.
bestgamma = 0
rsq = 0
for loggamma = -2., 1, 0.05 do begin
   gamma = 10.^loggamma
   svr.gamma = gamma
   res = svr.fit(tprofiles, tprop)
   nrsq = svr.score(vprofiles, vprop)
   if nrsq lt (rsq + 0.001) then break
   rsq = nrsq
   bestgamma = gamma 
   print, gamma, rsq
endfor
svr.gamma = bestgamma
print, '... done: gamma =', bestgamma, ' R^2 =', rsq

message, 'optimizing C ...', /inf
rsq =  0.
bestC = 0.01
for logC = -2., 3., 0.05 do begin
   C = 10.^logC
   svr.C = C
   res = svr.fit(tprofiles, tprop)
   nrsq = svr.score(vprofiles, vprop)
   if nrsq lt rsq then break
   rsq = nrsq
   bestC = C
endfor
svr.C = bestC
print, '... done: C =', bestC, ' R^2 =', rsq

message, 'building support vector machine ...', /inf
while rsq lt goal do begin
   nsets *= 1.5
   if nsets gt 200000L then break
   t0 = systime(1)
   self -> GenerateData, zrange = self.zrange, $
                         arange = self.arange, $
                         nrange = self.nrange, $
                         rad = self.rad, $
                         nsets = nsets, $
                         quiet = quiet

   profiles = self.data['profiles']
   tprofiles = profiles[*, 0:*:2]
   vprofiles = profiles[*, 1:*:2]
   prop = self.data[name+'p']
   tprop = prop[0:*:2]
   vprop = prop[1:*:2]
   message, 'optimizing gamma ...', /inf
   rsq = 0.
   bestgamma = svr.gamma
   for gamma = svr.gamma, 3.*svr.gamma, 0.05*svr.gamma do begin
      svr.gamma = gamma
      res = svr.fit(tprofiles, tprop)
      nrsq = svr.score(vprofiles, vprop)
      if nrsq lt (rsq + 0.001) then break
      rsq = nrsq
      bestgamma = gamma 
      print, gamma, rsq
   endfor
   svr.gamma = bestgamma
   print, '... done: gamma =', bestgamma, ' R^2 =', rsq
   message, 'optimizing C ...', /inf
   bestC = svr.C
   for C = 1.05*svr.C, 3.*svr.C, 0.05*svr.C do begin
      svr.C = C
      res = svr.fit(tprofiles, tprop)
      nrsq = svr.score(vprofiles, vprop)
      if nrsq lt (rsq + 0.001) then break
      rsq = nrsq
      bestC = C 
      print, C, rsq
   endfor
   svr.C = bestC
   print, '... done: C =', bestC, ' R^2 =', rsq
;   res = svr.fit(tprofiles, tprop)
;   rsq = svr.score(vprofiles, vprop)
   dt = systime(1) - t0
   print, 'Support vectors:', nsets, ' Elapsed time:', dt
   if dt ge 1800 then break
endwhile

end

;;;;;
;
; SVRfeature::Train
;
pro svrfeature::Train, goal = goal

COMPILE_OPT IDL2, HIDDEN

self.TrainProperty, 'z', goal = 0.99, /quiet
self.TrainProperty, 'a', goal = 0.99, /quiet
self.TrainProperty, 'n', goal = 0.99, /quiet

end

;;;;;
;
; SVRfeature::Save
;
pro svrfeature::Save, prefix = prefix

COMPILE_OPT IDL2, HIDDEN

if isa(prefix, 'string') then $
   self.prefix = prefix

fname = self.path + self.prefix + '.h5'
fid = h5f_create(fname)

;;; save SVRs
;; SVR group
svr_gid = h5g_create(fid, 'SVR')
;; SVR data
zsvr = self.pickle.dumps(self.data['zsvr'])
asvr = self.pickle.dumps(self.data['asvr'])
nsvr = self.pickle.dumps(self.data['nsvr'])
;; data types
zsvr_tid = h5t_idl_create(zsvr)
asvr_tid = h5t_idl_create(asvr)
nsvr_tid = h5t_idl_create(nsvr)
;; data spaces
zsvr_sid = h5s_create_simple(1)
asvr_sid = h5s_create_simple(1)
nsvr_sid = h5s_create_simple(1)
;; data sets
zsvr_id = h5d_create(svr_gid, 'zsvr', zsvr_tid, zsvr_sid)
asvr_id = h5d_create(svr_gid, 'asvr', asvr_tid, asvr_sid)
nsvr_id = h5d_create(svr_gid, 'nsvr', nsvr_tid, nsvr_sid)
;; write to HDF5 file
h5d_write, zsvr_id, zsvr
h5d_write, asvr_id, asvr
h5d_write, nsvr_id, nsvr
;; close SVR group
h5g_close, svr_gid

;;; Supporting data
data_gid = h5g_create(fid, 'data')

tid = h5t_idl_create(self.lambda)
sid = h5s_create_simple(1)
id = h5d_create(data_gid, 'lambda', tid, sid)
h5d_write, id, self.lambda
h5d_close, id

tid = h5t_idl_create(self.mpp)
sid = h5s_create_simple(1)
id = h5d_create(data_gid, 'mpp', tid, sid)
h5d_write, id, self.mpp
h5d_close, id

tid = h5t_idl_create(self.nm)
sid = h5s_create_simple(1)
id = h5d_create(data_gid, 'nm', tid, sid)
h5d_write, id, self.nm

h5g_close, data_gid

;;; Training ranges
range_gid = h5g_create(fid, 'ranges')

tid = h5t_idl_create(self.zrange)
sid = h5s_create_simple(2)
id = h5d_create(range_gid, 'zrange', tid, sid)
h5d_write, id, self.zrange
h5d_close, id

tid = h5t_idl_create(self.arange)
sid = h5s_create_simple(2)
id = h5d_create(range_gid, 'arange', tid, sid)
h5d_write, id, self.arange
h5d_close, id

tid = h5t_idl_create(self.nrange)
sid = h5s_create_simple(2)
id = h5d_create(range_gid, 'nrange', tid, sid)
h5d_write, id, self.nrange

h5g_close, range_gid

;;; Close the HDF5 database
h5f_close, fid

end

;;;;;
;
; SVRfeature::Load()
;
function svrfeature::Load, prefix = prefix, $
                           filename = filename
  
COMPILE_OPT IDL2, HIDDEN

if isa(filename, 'string') then $
   fn = filename

if isa(prefix, 'string') then $
   self.prefix = prefix

if ~isa(fn, 'string') then $
   fn = self.path + self.prefix + '.h5'

fn = file_search(fn, /test_regular, /test_read, count = count)
if count le 0 then begin
   message, 'Could not read ' + fn, /inf
   return, 0B
endif
fn = fn[0]

if ~h5f_is_hdf5(fn) then begin
   message, fn + ' is not an HDF5 file', /inf
   return, 0B
endif

self.lambda = h5_getdata(fn, 'data/lambda')
self.mpp = h5_getdata(fn, 'data/mpp')
self.nm = h5_getdata(fn, 'data/nm')
self.zrange = h5_getdata(fn, 'ranges/zrange')
self.arange = h5_getdata(fn, 'ranges/arange')
self.nrange = h5_getdata(fn, 'ranges/nrange')
self.data['zsvr'] = self.pickle.loads((h5_getdata(fn, 'SVR/zsvr'))[0])
self.data['asvr'] = self.pickle.loads((h5_getdata(fn, 'SVR/asvr'))[0])
self.data['nsvr'] = self.pickle.loads((h5_getdata(fn, 'SVR/nsvr'))[0])

zsvr = self.data['zsvr']
sv = zsvr.support_vectors_
self.rad = n_elements(sv[*, 0]) - 1L

self.trained = 1B

return, 1B
end

;;;;;
;
; SVRfeature::SetProperty
;
pro svrfeature::SetProperty, path = path, $
                             prefix = prefix, $
                             lambda = lambda, $
                             mpp = mpp, $
                             nm = nm, $
                             zrange = zrange, $
                             arange = arange, $
                             nrange = nrange, $
                             rad = rad
                             
COMPILE_OPT IDL2, HIDDEN

if isa(path, 'string') then begin
   mypath = file_search(path, /TEST_DIRECTORY, /MARK_DIRECTORY, count = count)
   if count gt 0 then $
      self.path = mypath[0]
endif

if isa(prefix, 'string') then $
   self.prefix = prefix

if isa(lambda, /scalar, /number) then begin
   self.trained = lambda eq self.lambda
   self.lambda = float(lambda)
endif

if isa(mpp, /scalar, /number) then begin
   self.trained = mpp eq self.mpp
   self.mpp = float(mpp)
endif

if isa(nm, /scalar, /number) then begin
   self.trained = nm eq self.nm
   self.nm = float(nm)
endif

if isa(zrange, /number, /array) then begin
   self.trained = total(zrange - self.zrange) eq 0
   self.zrange = zrange
endif

if isa(arange, /number, /array) then begin
   self.trained = total(arange - self.arange) eq 0
   self.arange = arange
endif

if isa(nrange, /number, /array) then begin
   self.trained = total(nrange - self.nrange) eq 0
   self.nrange = nrange
endif

if isa(rad, /scalar, /number) then begin
   self.trained = rad eq self.rad
   self.rad = rad
endif

end

;;;;;
;
; SVRfeature::GetProperty
;
pro svrfeature::GetProperty, path = path, $
                             prefix = prefix, $
                             zsvr = zsvr, $
                             nsvr = nsvr, $
                             asvr = asvr, $
                             lambda = lambda, $
                             mpp = mpp, $
                             nm = nm, $
                             zrange = zrange, $
                             arange = arange, $
                             nrange = nrange

COMPILE_OPT IDL2, HIDDEN

if arg_present(path) then $
   path = self.path

if arg_present(prefix) then $
   prefix = self.prefix

if arg_present(zsvr) then $
   zsvr = self.data['zsvr']

if arg_present(nsvr) then $
   nsvr = self.data['nsvr']

if arg_present(asvr) then $
   asvr = self.data['asvr']

if arg_present(lambda) then $
   lambda = self.lambda

if arg_present(mpp) then $
   mpp = self.mpp

if arg_present(nm) then $
   nm = self.nm

if arg_present(zrange) then $
   zrange = self.zrange

if arg_present(arange) then $
   arange = self.arange

if arg_present(nrange) then $
   nrange = self.nrange

end

;;;;;
;
; SVRfeature::InitPython()
;
function svrfeature::InitPython

COMPILE_OPT IDL2, HIDDEN

catch, error
if (error ne 0L) then begin
   catch, /cancel
   return, 0B
endif

self.pickle = pyimport('pickle')
self.svm = pyimport('sklearn.svm')

ok = isa(self.pickle, 'pythonobject') && isa(self.svm, 'pythonobject')

return, ok
end

;;;;;
;
; SVRfeature::Init()
;
function svrfeature::Init, path = path, $
                           prefix = prefix, $
                           lambda = lambda, $
                           mpp = mpp, $
                           nm = nm

COMPILE_OPT IDL2, HIDDEN

;;; path to data files
if ~isa(path, 'string') then begin
   fn = file_search(strsplit(!PATH, path_sep(/SEARCH_PATH), /EXTRACT) + $
                    '/svrfeature__define.pro')
   path = file_dirname(fn[0], /MARK_DIRECTORY) + 'svr/'
endif
mypath = file_search(path, /TEST_DIRECTORY, /MARK_DIRECTORY, count = count)
if count gt 0 then $
   self.path = mypath[0] $
else $
   return, 0B

self.prefix = isa(prefix, 'string') ? prefix : 'svr'

;;; Use Slither to import python modules
if ~self.InitPython() then $
   return, 0B

self.data = hash()

;;; Load default SVR data
if ~self.Load() then begin
   self.trained = 0B
   self.data['zsvr'] = self.svm.SVR()
   self.data['asvr'] = self.svm.SVR()
   self.data['nsvr'] = self.svm.SVR()
endif

return, 1B
end

;;;;;
;
; SVRfeature__define
;
pro svrfeature__define

COMPILE_OPT IDL2

struct = {SVRfeature, $
          INHERITS IDL_OBJECT, $
          svm:     obj_new(),  $ ; python SVM class
          pickle:  obj_new(),  $ ; python pickle class
          rad:     0L,         $ ; radius for azimuthal average
          lambda:  0.,         $ ; wavelength [micrometers]
          mpp:     0.,         $ ; magnification [micrometers/pixel]
          nm:      0.,         $ ; refractive index of medium
          data:    hash(),     $ ; training data
          arange:  [0., 0],    $ ; range of sphere radius [micrometers]
          nrange:  [0., 0],    $ ; range of refractive index values
          zrange:  [0., 0],    $ ; range of axial positions [pixels]
          trained: 0B,         $ ; flag
          path:    '',         $ ; directory containing data files
          prefix:  ''          $ ; filename for SVR data
         }

end
