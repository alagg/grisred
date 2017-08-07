;calculate continuum level
function get_cont_max,profile
  common cont,cont
  common line,line,cl_idx  
  common clevel,clevel
  
  bdry=(n_elements(profile)*0.02)>5
  if n_elements(clevel) le 1 then begin
    
    allidx=lindgen(n_elements(profile))
    wl=bin2wl(allidx)
    
                                ;forbidden areas for continuum
                                ;calculation (regions of possible
                                ;emission lines), and firs two and
                                ;last two pixels (sometimes noise in
                                ;first pixel)
    goodidx=where(wl lt 10828.809d or $
                  (wl gt 10832 and wl lt 10832.900d $
                   or wl gt 10833.5) and $ ;bump right of telluric line
                  (wl lt 10832.4 or wl ge 10832.7) and  $ ;He-lines 10830
                  (allidx gt bdry and allidx lt n_elements(allidx)-bdry) )
    
    if n_elements(goodidx) gt 5 and goodidx(0) ne -1 then idx=goodidx $
    else idx=allidx
    if total(finite(profile(idx))) ge 5 then $
      c=max(median((profile(idx))(*),3),/nan) $
    else c=max(profile(idx),/nan)
  endif else begin
                                ;add a few pixels around clevel.bin
                                ;values to avoid noise affection of
                                ;clevel determination
    add=10
    pm5bin=transpose(clevel.bin) ## lonarr(2*add+1)
    avgval=fltarr(n_elements(clevel))
    for i=0,n_elements(clevel.bin)-1 do begin
      pm5bin(*,i)=(((clevel(i).bin+lindgen(2*add+1)-5)>bdry)< $
                   (n_elements(profile)-bdry))
      imm=where(profile(pm5bin(*,i)) lt max(profile(pm5bin(*,i)),/nan) and $
                profile(pm5bin(*,i)) gt min(profile(pm5bin(*,i)),/nan))
      if imm(0) eq -1 then tpm5=profile(pm5bin(*,i)) else tpm5=profile(pm5bin(imm,i))
      ntp5=n_elements(tpm5)
      if ntp5 gt 7 then avgval(i)=total(median(tpm5,7))/ntp5 $
      else  avgval(i)=total(tpm5)/ntp5
    endfor
    c=max(avgval)
;    if total(finite(profile(pm5bin))) ge 9 then $
;      c=max(median((profile(pm5bin))(*),7),/nan) $
;    else c=max(profile(pm5bin),/nan)
    
;    plot,profile & plots,psym=2,color=1, pm5bin, profile(pm5bin)  
;    plot,profile(pm5bin),/yst,/xst & oplot,color=1,median((profile(pm5bin))(*),7) & for ic=0,n_elements(clevel)-1 do oplot,color=2+ic,!x.crange,[0,0]+avgval(ic)
  endelse
  
  c=c * cont.perc/100
  return,c
end

;use histogram for cont calculation
;new way to calc. cont level, (June 7th 2006)
;should be more robust and also cover the profiles with emission lines.
function get_cont,profile
  common cont,cont
  
                                ;since the continuum is usually a flat
                                ;line, the histogram should return the
                                ;continuum level
  
  n=n_elements(profile)
  nbin=(n)<200

                                ;take only points between 5% and 95%
                                ;of range (TIP1: 2 and 98%)
  if n lt 300 then nuse=fix([n*0.02,n*0.98]) $
  else nuse=fix([n*0.05,n*0.95])
  nuse=0>nuse<(n_elements(profile)-1)
  
  hist=histogram(median(profile(nuse(0):nuse(1)),3), $
                 nbin=nbin,/nan,omax=hmax,omin=hmin,loc=xhist)
;  xhist=findgen(nbin)/(nbin+1.)*(hmax-hmin)+hmin
  dummy=max(hist,imax)
  
                                ;determine where level falls to 90% of
                                ;max val and use this as imax
  maxlow=0.5*max(hist)
  imaxadd=min(where(hist(imax:*) le maxlow))
  imaxadd=max(where(hist(imax:*) gt maxlow))
  if imaxadd ne -1 then imax=imax+imaxadd
  
                                ;make a gauss fit to the histogram and
                                ;take the value 1/2 width above the
                                ;maximum
;  gf=gaussfit(xhist,hist,gpar,nterms=3)
                                ;add half width of Gaussian to max
                                ;position (only if fit looks
                                ;reasonable, i.e. the half width is
                                ;less than half of the entire range)
;  if gpar(0) le nbin/2. then imax=imax+gpar(0)/10.
  
;  imax=imax<(n_elements(xhist)-1)
  
  c=xhist(imax)
;  !p.multi=[0,1,2] & plot,profile & oplot,!x.crange,[0,0]+c,linestyle=1 & plot,xhist,hist,xrange=[0.8,1.2]*c & oplot,color=1,[0,0]+c,!y.crange & wait,.1
  if n_elements(cont) ne 0 then c=c * cont.perc/100
  return,c
  
end
