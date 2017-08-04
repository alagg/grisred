pro get_flat,fileff,flat,data=data,time=time,corr1=dd1,corr2=dd2,$
   lambda=lambda,order=order,ffnorm=ffnorm,periodfr=periodfr
  
if(keyword_set(data) eq 0) then data=0
if(keyword_set(dd1) eq 0) then dd1=0
if(keyword_set(dd2) eq 0) then dd2=0

file=fileff
dc=rfits_im(file,1,dd,hdr,nrhdr,/badpix)>0
time=param_fits(hdr,'UT      =',delimiter=':',vartype=1)
time=time(*,0)+(time(*,1)+time(*,2)/60.)/60.


tam=size(dc)
get_lun,unit
openr,unit,file

if(dd.bitpix eq 8) then begin
   datos=assoc(unit,bytarr(dd.naxis1,dd.naxis2),long(2880)*nrhdr)
endif else if(dd.bitpix eq 16) then begin
   datos=assoc(unit,intarr(dd.naxis1,dd.naxis2),long(2880)*nrhdr)
endif else if(dd.bitpix eq 32) then begin
   datos=assoc(unit,lonarr(dd.naxis1,dd.naxis2),long(2880)*nrhdr)
endif

ndc=8
for j=2,ndc do dc=dc+(rfits_im2(datos,dd,j,/badpix)>0)

dc=dc/ndc

npos=(dd.naxis3-8)/4
;npos=npos<15
time=mean(time[0:npos-1])

imorig=fltarr(4,tam(1),tam(2))
format=['(i2,$)','(i3,$)','(i4,$)','(i5,$)','(i6,$)']
print,npos
for i=0,npos-1 do begin
   print,i+1,format=format(fix(alog10(i+1)))
   for j=0,3 do begin
      dum=(rfits_im2(datos,dd,4*i+j+9,/badpix)-dc)>0
      imorig(j,*,*)=imorig(j,*,*)+dum
   endfor
endfor

free_lun,unit

print,' '

im=total(imorig,1)/4.
;for j=0,tam(2)-1 do im(*,j)=median(im(*,j),3)

if(n_elements(data) eq 1 or n_elements(dd1) eq 1 or n_elements(dd2) eq 1) then begin
   param_beams,im,data=data,corr1=dd1,corr2=dd2,lambda=lambda
endif

get_beams,im,sub1,sub2,data=data,corr1=dd1,corr2=dd2,lambda=lambda,/align 

tam1=size(sub1)
tam2=size(sub2)
norm1=median(sub1)
norm2=median(sub2)
norm=(norm1+norm2)/2
sub1=sub1/norm
sub2=sub2/norm
imorig=imorig/norm

naccum=param_fits(hdr,'ACCUMULA=',vartype=1)
ffnorm=norm/naccum/1.e5/npos  ;1.e5 = to get a QS cont. at ~1.e5 counts
                                ; in reduced data at disk center

;esp1=total(sub1,2)/tam1(2)
;esp2=total(sub2,2)/tam2(2)

esp1=fltarr(tam1[1])
esp2=fltarr(tam2[1])
for j=0,tam1[1]-1 do esp1[j]=median(sub1[j,*])
for j=0,tam2[1]-1 do esp2[j]=median(sub2[j,*])

esp=(esp1+esp2)/2.
cont=continuum(esp,5,/fts,lambda=lambda,order=order)
esp=esp*mean(cont)/cont
;;mcv
;fesp=abs(fft(esp-mean(esp)))
;z=where(fesp[70:300] eq max(fesp[70:300]))
;z=z[0]+70
;zmcv=z

;A. Lagg, Sep 2015
;do not take the first nb and the last nb elements for the fourier
;analysis. This should make the frequency determination for the fringe removal more reliable
nb=100
fesp=abs(fft(esp[nb:tam[1]-nb-1]-mean(esp[nb:tam[1]-nb-1])))
fmax=max(fesp[70:300-nb],imax)
z=float(imax+70)/float(tam[1]-2*nb)*tam[1]
;print,'Fringe wavelength [pix] (AL, MCV): ',tam[1]/[z,zmcv]
;plot,esp[400:500],/yst   & oplot,0.01*sin((findgen(100))*2*!pi/tam[1]*z)+1.01 & stop

periodfr=z                           
for k=-2,2 do esp=esp-franjas(esp,float(tam1[1])/(z+k),1)
esp[tam[1]-11:*]=esp[tam[1]-10]
esp[0:9]=esp[10]

ff1=fltarr(tam1(1),tam1(2))
ff2=fltarr(tam2(1),tam2(2))
x=findgen(tam[1])
for k=0,tam1(2)-1 do ff1(*,k)=interpol(esp,x,x+dd1(k))
for k=0,tam2(2)-1 do ff2(*,k)=interpol(esp,x,x+dd2(k))
;for k=0,tam1(2)-1 do ff1(*,k)=desp1d(esp,-dd1(k))
;for k=0,tam2(2)-1 do ff2(*,k)=desp1d(esp,-dd2(k))

flat=fltarr(4,tam(1),tam(2))+1
for j=0,3 do flat[j,*,data[2]:data[3]]=imorig[j,*,data[2]:data[3]]/ff1
for j=0,3 do flat[j,*,data[6]:data[7]]=imorig[j,*,data[6]:data[7]]/ff2

z=where(flat le 0)
if(z(0) ne -1) then flat(z)=1.

flat[*,*,data[2]-2]=flat[*,*,data[2]]
flat[*,*,data[2]-1]=flat[*,*,data[2]]
flat[*,*,data[3]+1]=flat[*,*,data[3]]
flat[*,*,data[3]+2]=flat[*,*,data[3]]
flat[*,*,data[6]-2]=flat[*,*,data[6]]
flat[*,*,data[6]-1]=flat[*,*,data[6]]
flat[*,*,data[7]+1]=flat[*,*,data[7]]
flat[*,*,data[7]+2]=flat[*,*,data[7]]

return

ny=tam1(2)
;ny=round(data(3))-round(data(2))+1
for j=0,3 do begin
   for k=0,tam1(1)-1 do begin

      ystart=data[2]
      yend=data[3]
      y=ystart+findgen(ny)*(yend-ystart)/(ny-1.)
      y1=round(ystart)-1
      y2=round(yend)+1
      ynew=y1+indgen(y2-y1+1)
      flat(j,k,y1:y2)=interpol(ff1(k,*),y,ynew)
      flat(j,k,y1:y2)=imorig(j,k,y1:y2)/flat(j,k,y1:y2)
   endfor
   for k=0,tam2(1)-1 do begin
      ystart=data[6]
      yend=data[7]
      y=ystart+findgen(ny)*(yend-ystart)/(ny-1.)
      y1=round(ystart)-1
      y2=round(yend)+1
      ynew=y1+indgen(y2-y1+1)
      flat(j,k,y1:y2)=interpol(ff2(k,*),y,ynew)
      flat(j,k,y1:y2)=imorig(j,k,y1:y2)/flat(j,k,y1:y2)
   endfor
endfor

z=where(flat le 0)
if(z(0) ne -1) then flat(z)=1.

return
end
