function lstsqfit,x,y,yfit,tol=tol

if(keyword_set(tol) eq 0) then tol=1.e-6

d=transpose(x)#x

if(n_elements(d) eq 1) then begin
   coef=total(x*y)/total(x*x)
   yfit=coef*x
   ndat=n_elements(y)
   cov=total((yfit-y)*(yfit-y))/ndat
   coef=[coef,sqrt(cov/total(x*x))]
   return,coef
endif else begin
   d=invert_svd(d,tol=tol)
   tam=size(d)
   npar=tam(1)
   ndat=n_elements(y)
   sig=fltarr(npar)
   for j=0,npar-1 do sig(j)=d(j,j)
   coef=d#transpose(x)#y
   yfit=x#coef
   cov=total((yfit-y)*(yfit-y)/(ndat-npar+1))
   coef=[[coef],[sqrt(sig*cov)]]
endelse

return,coef
end

