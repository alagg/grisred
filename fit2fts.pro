;+
; NAME
;
;    FIT2FTS
;
; EXPLANATION
;
;    The function returns a continuum correction curve. This curve is a
;    polynomial and is retrieved from fitting the FTS spectrum to a GRIS
;    flat field profile.
;    The fitting is based on the genetic algorithm PIKAIA and has the
;    following free parameters:
;    WL-offset, WL-dispersion (both in Angstrom), the FWHM of a
;    spectral PSF (Gaussian), and spectral straylight contribution
;    (white light).
;
;    The fitting is only done for regions without telluric
;    blends by setting the weight for the computation of the fitness
;    (=1/chi^2) to zero, whenever the telluric FTS spectrum is below 0.9.
;    
;
; INPUTS
;
;    iprof - Stokes I profile
;    lambda - approx. central wavelength of observed
;             region in Angstrom, e.g. lambda=15650.
;    order - grating order. Improtant for computing initial guess for
;             WL-fitting
;    npoly - degree of fitted polynomial 
;
; OUTPUTS
;
;     return value - continuum correction curve (same size as input
;                    Stokes-I profile) 
; Optional output:
;     fit - structure containing the fit paramaters (WL-offset,
;           WL-dispersion, FWHM of spectral Gaussian PSF (in pixel and
;           Angstrom), spectral straylight,), the fitted WL-vector, the
;           input Stokes I-profile normalized to the continuum,
;           the degraded FTS profile, the continuum correction polynomial
;           vector, the polynomial coefficients and the degree of the
;           fitted polynomial
;
; CALLS
;
;    This routine is called from the procedure continuum.pro (GRIS
;    data reduction package):
;    cont=fit2fts(esp,show=show,lambda=lambda,order=order)
;
; HISTORY
;
;    Ver.1, 08-Aug-2017, A. Lagg
;-

function mpfit_fts,par,xval=x,yval=y,errval=err,wgt=wgt,_extra=_extra, $
                   show=show,store=store
  common fts,ftsfull,fts,hdr
  common prof,prof,fitset,fitpar,fitval,fitvalold
  
  if n_elements(y) eq 0 then y=prof
  nwl=n_elements(y)
  if n_elements(x) eq 0 then x=dindgen(nwl)
  
  np=n_elements(par)
;par[0] = WL offset
;par[1] = WL dispersion 
  wlobs=x*par[1]+par[0]
  
;interpolate fts to same grid as observation
  ftsip=interpol(fts[4,*],fts[1,*],wlobs)
  atmip=interpol(fts[3,*],fts[1,*],wlobs)
  wgt=atmip ge 0.95
                                ;remove first & last nb points from
                                ;fitting (at
;  nb=15
;  wgt[0:nb-1]=0
;  wgt[nwl-nb:*]=0
  
  ftsdeg=ftsip
  if np gt 2 then begin
;spectral convolution
;par[2] = Gaussian width of spectral PSF (FWHM)
    fwhm_pix=abs(par[2])/par[1]
    sigma=fwhm_pix/(2*sqrt(2*alog(2.)))
    nw=5
    kernelx=findgen(fix(nw*fwhm_pix/2)*2+1)-fix(nw*fwhm_pix/2)
    kernel=exp(-0.5*(kernelx/sigma)^2)
    kernel=kernel/total(kernel)
    ftscon=convol(ftsip,kernel,/edge_truncate)
    
    if fitset.strayfit then if np gt 3 then begin
;add white light
;par[3] = spectral straylight (constant offset)
      ftsdeg=ftscon*(1.-par[3])+par[3]
    endif
  endif else poly=ftsdeg*0.
  
  if fitset.polyfit eq 1 then begin
    npol=fitset.npoly
    nfin=finite(ftsdeg/prof)
    if total(nfin) lt nwl then message,'WL-range for FTS fitting too small.'
    polpar=poly_fit(wlobs-mean(wlobs),1-ftsdeg/prof,npol, $
                    measure_errors=sqrt(abs(1-ftsdeg/prof)), $
                    yfit=poly,chi=chi,/double,status=status)
    ftsdeg=ftsdeg*(1+poly)
  endif else poly=ftsdeg*0.
  
  if keyword_set(store) or keyword_set(show) then begin
    fitval={wlobs:wlobs,prof:prof, $
            ftsdeg:ftsdeg,poly:poly+1,fitness:0., $
            wloff:par[0],wlbin:par[1],fwhm_a:par[2],fwhm_pix:fwhm_pix, $
            stray:par[3],polypar:polpar,npoly:npol}
  endif

  return,ftsdeg
end

function sclpar,par,reverse=reverse
  common pi,pi
  
  if keyword_set(reverse) then $
    retpar=par*(pi.limits[1]-pi.limits[0])+pi.limits[0] $
  else $
    retpar=(par-pi.limits[0])/(pi.limits[1]-pi.limits[0])
  
  return, retpar
end

function fts_chi,par,show=show,store=store
  common fts,ftsfull,fts,hdrc
  common prof,prof,fitset,fitpar,fitval,fitvalold
  
  ftsdeg=mpfit_fts(par,wgt=wgt,show=show,store=store)
  chisqr=total(wgt*(prof-ftsdeg)^2,/nan)
  if keyword_set(show) then begin
    fitness=1./chisqr
    print,'Fitness: ',fitness
    fitval.fitness=fitness
    !p.multi=[0,1,2]
    yrg=minmaxp([fitval.prof,fitval.ftsdeg])
    plot,fitval.wlobs,fitval.prof,/xst,/yst,thick=1,xtitle='Wavelength [A]', $
         title='FTS-fit, Fitness '+string(fitval.fitness,format='(f10.4)'), $
         yrange=yrg
    oplot,linestyle=2,fitval.wlobs,fitval.ftsdeg,thick=2
    
    yrg=minmaxp(fitval.poly)
    if n_elements(fitvalold) ne 0 then yrg=minmaxp([fitval.poly,fitvalold.poly])
    plot,fitval.wlobs,fitval.poly,/xst,/yst,thick=1,xtitle='Wavelength [A]', $
         title='Continuum correction (polynomial fit order '+ $
         strcompress(/remove_all,(string(fitval.npoly)))+')',yrange=yrg
    if n_elements(fitvalold) ne 0 then begin
      oplot,fitvalold.wlobs,fitvalold.poly,linestyle=1
      xyouts,0,0,/normal,'dotted = previous FTS fit'
    endif
  endif
  return,chisqr
end

function func,n,spar
  return,1/fts_chi(sclpar(spar,/reverse))
end

pro pikcall,par,statpar,f,status
  common prof,prof,fitset,fitpar,fitval,fitvalold
  
;================================================================  
                                ;call pikaia
;================================================================  
  common share2,iv,iy,idum2
  
  pikaia,/init
;
;     Sample driver program for pikaia.f
;
  NTAB=32 
  idum2=123456789
  iv=lonarr(NTAB) 
  iv(*)= 0
  iy=0
;  read, "Random number seed (I*4)? ",seed
  seed=long(randomu(undefinedvar)*1e5)
  seed = -abs(long(seed))
  urand_init,seed
;
;     Set control variables (use defaults)
;
  ctrl = fltarr(12)
  ctrl(*) = -1
  
                                ;set some values optimized for the
                                ;he-line problem
  ctrl(0)=fitset.npop
  ctrl(1)=fitset.niter
  ctrl(2)=5                     ;significant bits (6 is default and should be
                                ;sufficient for real precission
  ctrl(3)=0.85                  ;0.60 ;0.85
  ctrl(4)=2                     ;idl version only allows imut=1 or 2,
                                ;5 seems to converge better, but is
                                ;only available for the fortran version!
  ctrl(5)=0.01                  ;0.005 ;0.01
  ctrl(6)=0.001                 ;0.0005 ;0.001
  ctrl(7)=0.05                  ; 0.01 ;0.05
;    ctrl(11)=2   ;verbose
  ctrl(11)=0                    ;verbose
  
  n=n_elements(par)
  spar=sclpar(par)
  f=func(n,spar)
  pikaia,n,ctrl,spar,statpar,f,status
  par=sclpar(spar,/reverse)

end
;PURPOSE: fits a polynomial as continuum correction curve
;From the order and the lambda keyword 
;INPUT:
;iprof - Stokes I profile
;lambda - approx. central wavelength of observed
;         region in Angstrom, e.g. lambda=15650.
;order - grating order. Improtant for computing initial guess for
;        WL-fitting
;npoly - degree of fitted polynomial 
;OUTPUT:
;return value - continuum correction curve (same size as input Stokes-I profile)
;Optional output:
;fit - structure containing the fit paramaters (WL-offset,
;      WL-dispersion, FWHM of spectral Gaussian PSF (in pixel and
;      Angstrom), spectral straylight,), the fitted WL-vector, the
;      input Stokes I-profile normalized to the continuum,
;      the degraded FTS profile, the continuum correction polynomial
;      vector, the polynomial coefficients and the degree of the
;      fitted polynomial

function fit2fts,iprof,show=show,lambda=lambda,order=order,npoly=npoly, $
                 fit=retfitval
  common fts,ftsfull,fts,hdrc
  common prof,prof,fitset,fitpar,fitval,fitvalold
  common pi,pi
  common oldfit,oldfitpar
  
  show=keyword_set(show)
  print,strjoin(replicate('=',80))
  print,'FTS-fitting, please wait...'
  if keyword_set(show) eq 0 then $
    print,'   (use keyword gris_v5,...,/show to watch fitting procedure)'
  
  
  if n_elements(npoly) eq 0 then npoly=29
  if n_elements(fitset) eq 0 then $
    fitset={wlfit:1b,psffit:1b,strayfit:1b,polyfit:1b, $
            niter:200,npop:512,npoly:npoly}
  
  if n_elements(fts) eq 0 then begin
    ftsdir=getenv('FTSDIR')
    if ftsdir eq '' then ftsdir='.'
    fits=ftsdir+'/fts_combined.fits'
    print,'Reading FTS spectrum: ',fits
    ftsfull=readfits(fits,hdr)    
  endif
  wl=lambda                     ;
  wlrg=50.                      ;    
  wlbin=(18e-3)/redir_gregor_func(lambda,order,show=show)

  infts=where(ftsfull[1,*] ge wl-wlrg and ftsfull[1,*] le wl+wlrg)
  if infts[0] eq -1 then message,'wavelength not coverd by FTS spectrum'
  fts=ftsfull[*,infts]
  
  par=[wl,wlbin ,0.15,.15]
  lim=[20.,0.002,.10,.10]  
  prof=iprof/get_cont(iprof)    ;1.03;mean(histo_scale(iprof,perc=50))
  
  
  
  nwl=n_elements(prof)
  x=dindgen(nwl)
  functargs={xval:x,yval:prof,errval:prof}
  pi = replicate({name:'',value:0.D, fixed:0,limits:[0.D,0]}, $
                 n_elements(par))
  pi.name=['WL-offset','WL-dispersion','spectral resolution (FWHM)', $
           'spectral straylight']
  pi.value=par
  pi.limits=transpose([[(par-lim)>0.01],[par+lim]])
  
                                ;store for faster fitting
                                ;of subsequent flatfields
  if n_elements(oldfitpar) eq 0 then oldfitpar={order:0,wl:0d,pi:pi}
  
                                ;skip the first fitting step if
                                ;already done for same WL and order
  if oldfitpar.order eq order and abs(oldfitpar.wl-lambda) le 0.5 then begin
    if keyword_set(show) then $
      print,'Using fit results from previous fitting:', $
            oldfitpar.wl,oldfitpar.order
    pi=oldfitpar.pi
    pi[0].limits=pi[0].value+[-1,1]/2.   ; narrow bounds for WL-off fitting
    pi[1].limits=pi[1].value+[-1,1]*1e-3 ; narrow bounds for WL-disp fitting
  endif else begin
    
                                ;first fit with low order polynomial
                                ;(increases likelihood to find correct
                                ;WL-offset and dispersion)  
    fitset.polyfit=1b
    fitset.npoly=5
    fitset.niter=100             ;pikaia iterations
    fitset.npop=256             ;pikaia populations
    par=pi.value
    pikcall,par,statpar,f,status
    pi.value=par
    fitness=1./fts_chi(par,show=show)
                                ;then fit the full thing
    pi[0].limits=pi[0].value+[-1,1]      ; narrow bounds for WL-off fitting
    pi[1].limits=pi[1].value+[-1,1]*1e-3 ; narrow bounds for WL-disp fitting
  endelse
  fitset.polyfit=1b
  fitset.niter=60               ;pikaia iterations
  fitset.npop=512               ;pikaia populations
  fitset.npoly=npoly            ;degree of polynomial used for fitting
  pikcall,par,statpar,f,status
  pi.value=par
  fitness=1./fts_chi(par,show=show,/store)
  print,strjoin(replicate('=',80))
  print,'FTS-fitting information:'
  for i=0,n_elements(pi)-1 do $
    print,format='(a30,'':'',f18.8)',pi[i].name,pi[i].value
  print,format='(a30,'':'',f18.8)','Fitness',fitness
  print,strjoin(replicate('=',80))
  
  if fitness lt 5 then begin
    print,'FTS fitting ended with low fitness.'
    print,'Please check the quality of the fit in the displayed figure.'
    print,'Continue y/[n]?'
    fitness=1./fts_chi(par,show=1)    
    repeat tkey=strupcase(get_kbrd()) until tkey ne ''
    if tkey ne 'Y' then stop
    print,strjoin(replicate('=',80))
  endif
  
  oldfitpar.order=order
  oldfitpar.wl=lambda
  oldfitpar.pi=pi
  
  fitvalold=fitval
  
  retfitval=fitval
  return,fitval.poly
end
