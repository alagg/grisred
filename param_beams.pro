pro param_beams,im_in,data=data,corr1=dd1,corr2=dd2,lambda=lambda

if(keyword_set(data) eq 0) then data=0
if(keyword_set(dd1) eq 0) then dd1=0
if(keyword_set(dd2) eq 0) then dd2=0

im=im_in/mean(im_in)
if(n_elements(data) ne 8) then data=get_limits(im,lambda)
tam=size(im)

get_beams,im,sub1,sub2,data=data,lambda=lambda

tam1=size(sub1)
tam2=size(sub2)

nbox=50
esp1=total(sub1(*,tam1(2)/2-nbox:tam1(2)/2+nbox),2)/(2.*nbox+1.)
esp2=total(sub2(*,tam2(2)/2-nbox:tam2(2)/2+nbox),2)/(2.*nbox+1.)
esp1=median(esp1,7)
esp2=median(esp2,7)

esp1=esp1/continuum(esp1,5)
esp2=esp2/continuum(esp2,5)

esp=(esp1+esp2)/2.
esp=esp/continuum(esp,5)

np=7
s=findgen(2*np+1)-np
c=fltarr(2*np+1)
dd1=fltarr(tam1(2))
x=findgen(tam1[1])
xx=fltarr(3,3)
xx(*,0)=1.
xx(*,1)=[-1.,0.,1.]
xx(*,2)=[1.,0.,1.]
maxiter=30
;;;; print,'beam 1'
for j=0,tam1(2)-1 do begin
   sub1n= sub1[*,j]/continuum(sub1[*,j],5)
   if(j ne 0) then dd1(j)=dd1(j-1)
   iter=0
   dref=1
   while(abs(dref) ge 0.01 and iter le maxiter) do begin
      iter=iter+1
      for k=0,2*np do begin
         c(k)=correlate(esp1,shift(interpolate(sub1n,x-dd1[j]),s(k)))
;         c(k)=correlate(esp1,shift(desp1d(sub1n,dd1(j)),s(k)))
      endfor
      z=where(c eq max(c))
      z=z(0)
      if(z eq 0 or z eq 2*np) then begin
         dref=0
      endif else begin
         coef=lstsqfit(xx,c(z-1:z+1))
         dref=s[z]-coef(1,0)/2./coef(2,0)
      endelse
      dd1(j)=dd1(j)+dref
   endwhile
endfor

;dd1=-(findgen(tam1(2))-tam1(2)/2)*0.012

dd2=fltarr(tam2(2))
;;;; print,'beam 2'
for j=0,tam2(2)-1 do begin
   sub2n= sub2(*,j)/continuum(sub2[*,j],5)
   if(j ne 0) then dd2(j)=dd2(j-1)
   iter=0
   dref=1
   while(abs(dref) ge 0.01 and iter le maxiter) do begin
      iter=iter+1
      for k=0,2*np do begin
         c(k)=correlate(esp1,shift(interpolate(sub2n,x-dd2[j]),s(k)))
;         c(k)=correlate(esp1,shift(desp1d(sub2n,dd2(j)),s(k)))
      endfor
      z=where(c eq max(c))
      z=z(0)
      if(z eq 0 or z eq 2*np) then begin
         dref=0
      endif else begin
         coef=lstsqfit(xx,c(z-1:z+1))
         dref=s[z]-coef(1,0)/2./coef(2,0)
      endelse
      dd2(j)=dd2(j)+dref
   endwhile
endfor
;dd2b=-(findgen(tam2(2))-tam2(2)/2)*0.002

dd1=median(dd1,11)
dd2=median(dd2,11)

coef=poly_fit(findgen(n_elements(dd1)),dd1,5,yfit1)
dd1=yfit1
coef=poly_fit(findgen(n_elements(dd2)),dd2,5,yfit2)
dd2=yfit2

sub1r=fltarr(tam1(1),tam1(2))
sub2r=fltarr(tam2(1),tam2(2))

for j=0,tam1(2)-1 do sub1r(*,j)=interpolate(sub1(*,j),x-dd1(j))
for j=0,tam2(2)-1 do sub2r(*,j)=interpolate(sub2(*,j),x-dd2(j))
;for j=0,tam1(2)-1 do sub1r(*,j)=desp1d(sub1(*,j),dd1(j))
;for j=0,tam2(2)-1 do sub2r(*,j)=desp1d(sub2(*,j),dd2(j))
esp1r=total(sub1r,2)/tam1(2)
esp2r=total(sub2r,2)/tam2(2)
esp1r=median(esp1r,7)
esp2r=median(esp2r,7)

;;;; print,'both beams'

iter=0
dref=1
while(abs(dref) ge 0.01 and iter le maxiter) do begin
   iter=iter+1
   for k=0,2*np do c(k)=correlate(esp1r,shift(esp2r,s(k)))
   z=where(c eq max(c))
   z=z(0)
   if(z eq 0 or z eq 2*np) then begin
      dref=0
   endif else begin
      coef=lstsqfit(xx,c(z-1:z+1))
      dref=s[z]-coef(1,0)/2./coef(2,0)
   endelse
   dd2=dd2+dref
   for j=0,tam2(2)-1 do sub2r(*,j)=interpolate(sub2(*,j),x-dd2(j))
;   for j=0,tam2(2)-1 do sub2r(*,j)=desp1d(sub2(*,j),dd2(j))
   esp2r=total(sub2r,2)/tam(2)
   esp2r=median(esp2r,7)
endwhile

return
end
