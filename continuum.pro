function continuum,esp_in,ngrad,fts=fts,lambda=lambda,order=order,show=show

  norm=mean(esp_in)
  esp=esp_in/norm

  if keyword_set(fts) then begin
    cont=fit2fts(esp,show=show,lambda=lambda,order=order)
  endif else begin
    imin=where(esp eq min(esp))
    imin=imin(0)
    n=n_elements(esp)
    z=(findgen(n)-imin)/n

    x=z
    esp2=esp
    if(ngrad gt 0) then begin
      coef=poly_fit(x,esp2,ngrad,yfit)
      cont=poly(z,coef)
    endif else begin
      cont=fltarr(n)+mean(esp2) 
      yfit=cont  
    endelse

    cota=0.98
    zbad=where(esp lt cont-1.5*std(esp2-yfit) and esp lt cota*cont)

    xbad=fltarr(n)+1

    maxiter=30
    iter=0
    zbadold=0

;print,zbad(0)
    while(zbad(0) ne -1 and max(abs(zbad-zbadold)) ne 0 and iter le maxiter) do begin

      iter=iter+1
      xbad(zbad)=-1
      zgood=where(xbad eq 1)
      xgood=x(zgood)
      esp2=esp(zgood)
      if(ngrad gt 0) then begin
        coef=poly_fit(xgood,esp2,ngrad,yfit)
        cont=poly(z,coef)
      endif else begin
        yfit=fltarr(n_elements(esp2))+mean(esp2)   
        cont=fltarr(n)+mean(esp2)   
      endelse

      zbadold=zbad
      zbad=where(esp lt cont-1.5*std(esp2-yfit) and esp lt cota*cont)

    endwhile
  endelse
                                ;plot,esp & oplot,color=1,esp/(cont*norm) & oplot,color=2,esp/(cont2*norm)
  return,cont*norm
end   
