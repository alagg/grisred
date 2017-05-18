pro get_calib,filecal,dmod1,dmod2,data=data,lambda=lambda,zbad=zbad,$
   timecal=timecal,plot=plot,date=date,pzero=pzero,$
   flat1=flat1,flat2=flat2,timeflat1=timeff1,timeflat2=timeff2

if(keyword_set(plot) eq 1) then plot=1 else plot=0
if(keyword_set(pzero) eq 0) then pzero=-1.2

dc=rfits_im(filecal,1,dd,hdr,nrhdr,/badp)
for j=2,8 do dc=dc+rfits_im(filecal,j,/badp)
dc=dc/8.
tam=size(dc)

date=param_fits(hdr,'DATE-OBS=',delimiter='-',vartype=1)
dum=date[0]
date[0]=date[2]
date[2]=dum

tam=size(dc)

;VTT
;lamref=[10830.,11500.,12500.,15648.,16500.,17500.]
;deltaref=[316.9+6,274.1,222.0,90.6,63.9,40.3]
;if(keyword_set(delta) eq 0) then delta=interpol(deltaref,lamref,lambda)
;print,delta

i1=data(0)
i2=data(1)
j1=data(2)
j2=data(3)  
i3=data(4)
i4=data(5)
j3=data(6)
j4=data(7)

npos=(dd.naxis3-8)/4

cal_type=param_fits(hdr,'MEASURE =')
pos_tel=strpos(cal_type,'telescope')
pos_ins=strpos(cal_type,'instrumental')

if(pos_tel eq -1 and pos_ins ne -1) then begin	; instrumental calibration
   thpol=param_fits(hdr,'INSPOLAR=',vartype=3)+pzero
   thret=param_fits(hdr,'INSRETAR=',vartype=3)

   delta=-409.288 + 786.267e4/lambda

endif else if (pos_tel ne -1 and pos_ins eq -1) then begin	; telescope calibration
   thpol=param_fits(hdr,'TELPOLAR=',vartype=3)+pzero
   thret=param_fits(hdr,'TELRETAR=',vartype=3)

   lamref=[10830.,10938.,12818.,15000.,15650.,16000.,16500.,17000.,17500.]
   deltaref=[88.8,91.1,90.7,86.8,85.7-1.5,83.1,81.0,79.5,76.8]
   coef=poly_fit(lamref[1:*],deltaref[1:*],2)
   delta=poly(lambda,coef)
endif else begin
   print,'unknown calibration file. Please check!
   stop
endelse

timecal=param_fits(hdr,'UT      =',delimiter=':',vartype=1) 
timecal=timecal(*,0)+(timecal(*,1)+timecal(*,2)/60.)/60.

if(keyword_set(flat1) eq 0 and keyword_set(flat2) eq 0) then begin
   flat=fltarr(4,tam[1],tam[2])+1.
endif else if(keyword_set(flat1) ne 0 and keyword_set(flat2) eq 0) then begin
   flat=flat1
endif else if(keyword_set(flat1) eq 0 and keyword_set(flat2) ne 0) then begin
   flat=flat2
endif else begin
   factor=(mean(timecal)-timeff1)/(timeff2-timeff1)
   flat=(1-factor)*flat1+factor*flat2
endelse

q=0.

luzin=fltarr(npos,4)
uno=[1,0,0,0]


elaz=get_elaz(date[2],date[1],date[0],timecal[0])
tel_el0=elaz[0]
tel_az0=elaz[1]

for j=0,npos-1 do begin
   elaz=get_elaz(date[2],date[1],date[0],timecal[j])
   tel_el=elaz[0]
   tel_az=elaz[1]
   
   if(pos_tel eq -1 and pos_ins ne -1) then begin	; instrumental calibration
      delta_angle=-tel_az
   endif else if (pos_tel ne- 1 and pos_ins eq -1) then begin	; telescope calibration
      delta_angle=(tel_el-tel_el0)*0+(tel_az-tel_az0)
   endif
   
   luzin(j,*)=retarder(thret(j)-delta_angle,delta)# $
              impolariz2(q,thpol(j)-delta_angle)#uno
endfor   

nslit=j2-j1+1
nx=i2-i1+1

cuad1=fltarr(i2-i1+1,j2-j1+1)
cuad2=fltarr(i4-i3+1,j4-j3+1)
haz1=fltarr(npos,4,nslit)
haz2=fltarr(npos,4,nslit)

get_lun,unit
openr,unit,filecal
      
if(dd.telescope eq 'SVST') then begin
   if(dd.bitpix eq 8) then begin
      datos=assoc(unit,bytarr(dd.naxis2,dd.naxis1),long(2880)*nrhdr)
   endif else if(dd.bitpix eq 16) then begin   
      datos=assoc(unit,intarr(dd.naxis2,dd.naxis1),long(2880)*nrhdr)
   endif else if(dd.bitpix eq 32) then begin   
      datos=assoc(unit,lonarr(dd.naxis2,dd.naxis1),long(2880)*nrhdr)
   endif
endif else begin   
   if(dd.bitpix eq 8) then begin
      datos=assoc(unit,bytarr(dd.naxis1,dd.naxis2),long(2880)*nrhdr)
   endif else if(dd.bitpix eq 16) then begin   
      datos=assoc(unit,intarr(dd.naxis1,dd.naxis2),long(2880)*nrhdr)
   endif else if(dd.bitpix eq 32) then begin   
      datos=assoc(unit,lonarr(dd.naxis1,dd.naxis2),long(2880)*nrhdr)
   endif
endelse 

for i=0,npos-1 do begin
   for j=0,3 do begin
      im=(rfits_im2(datos,dd,4*i+9+j,/badpix)-dc)>0
      im=im/flat[j,*,*]
      cuad1=im(i1:i2,j1:j2)
      cuad2=im(i3:i4,j3:j4)
      haz1(i,j,*)=total(cuad1,1)/nx
      haz2(i,j,*)=total(cuad2,1)/nx
   endfor
endfor

close,unit
free_lun,unit

mm=fltarr(4,4,nslit)
mm1=fltarr(4,4,nslit)
mm2=fltarr(4,4,nslit)
dmod=fltarr(4,4,nslit)
dmod1=fltarr(4,4,nslit)
dmod2=fltarr(4,4,nslit)
aj1=fltarr(npos,4,nslit)
aj2=fltarr(npos,4,nslit)

sigma1=fltarr(nslit)
sigma2=fltarr(nslit)
cor=[[1,1,1,1],[-1,-1,-1,-1],[-1,-1,-1,-1],[-1,-1,-1,-1]]

if(keyword_set(zbad) eq 0) then $
   zbad=getbadpoints(filecal,data=data)

zgood=indgen(npos)
if(zbad[0] ne -999) then begin 
   zgood[zbad]=-999
   zgood=zgood(where(zgood ne -999))
endif

if(plot eq 1) then loadct,2,/silent

for slit=0,nslit-1 do begin
   for j=0,3 do begin
      coef=lstsqfit(luzin[zgood,*],haz1[zgood,j,slit],yfit1)

      if(plot eq 1) then begin
         plot,haz1[*,j,slit],tit=strtrim(slit,2) & oplot,zgood,yfit1,lin=2,color=80
         wait,0.2
      endif

      mm1[j,*,slit]=coef[*,0]
      aj1[zgood,j,slit]=yfit1
      coef=lstsqfit(luzin[zgood,*],haz2[zgood,j,slit],yfit2)

      if(plot eq 1) then begin
         plot,haz2[*,j,slit],tit=strtrim(slit,2) & oplot,zgood,yfit2,lin=2,color=80
         wait,0.2
      endif

      mm2[j,*,slit]=coef[*,0]
      aj2[zgood,j,slit]=yfit2
   endfor

   sigma1[slit]=std(haz1[zgood,*,slit]-aj1[zgood,*,slit])
   sigma2[slit]=std(haz2[zgood,*,slit]-aj2[zgood,*,slit])

   mm1[*,*,slit]=mm1[*,*,slit]/max(mm1[*,0,slit])
   mm2[*,*,slit]=mm2[*,*,slit]/max(mm2[*,0,slit])
   mm2[*,*,slit]=mm2[*,*,slit]*cor
   mm[*,*,slit]=(mm1[*,*,slit]+mm2[*,*,slit])/2.
   
   dmod1[*,*,slit]=invert(mm1[*,*,slit])
   norm=total(dmod1[0,*,slit])
   mm1[*,*,slit]=mm1[*,*,slit]*norm
   dmod1[*,*,slit]=dmod1[*,*,slit]/norm
   for j=1,3 do dmod1[j,*,slit]=dmod1[j,*,slit]- $
      total(dmod1[j,*,slit])/total(dmod1[0,*,slit])*dmod1[0,*,slit]

   dmod2[*,*,slit]=invert(mm2[*,*,slit])
   norm=total(dmod2[0,*,slit])
   mm2[*,*,slit]=mm2[*,*,slit]*norm
   dmod2[*,*,slit]=dmod2[*,*,slit]/norm
   for j=1,3 do dmod2[j,*,slit]=dmod2[j,*,slit]- $
      total(dmod2[j,*,slit])/total(dmod2[0,*,slit])*dmod2[0,*,slit]

   dmod(*,*,slit)=invert(mm[*,*,slit])
   norm=total(dmod[0,*,slit])
   mm[*,*,slit]=mm[*,*,slit]*norm
   dmod[*,*,slit]=dmod[*,*,slit]/norm

endfor

if(plot eq 1) then loadct,0,/silent

;loadct,2,/silent
;plot,sigma2 & oplot,sigma1,color=80
;loadct,0,/silent

sigma=mean([sigma1,sigma2])
print,'sigma1,sigma2,media = ',mean(sigma1),mean(sigma2),mean([sigma1,sigma2])

avdmod=total(dmod,3)/nslit
avdmod1=total(dmod1,3)/nslit
avdmod2=total(dmod2,3)/nslit
eps=fltarr(5,nslit)
eps1=fltarr(5,nslit)
eps2=fltarr(5,nslit)
for j=0,nslit-1 do eps(*,j)=effic(mm[*,*,j])
for j=0,nslit-1 do eps1(*,j)=effic(mm1(*,*,j))
for j=0,nslit-1 do eps2(*,j)=effic(mm2(*,*,j))
aveps=total(eps,2)/nslit
aveps1=total(eps1,2)/nslit
aveps2=total(eps2,2)/nslit 
print,'avdmod'
print,transpose(avdmod),format='(4f8.4)'
print,'avdmod1'
print,transpose(avdmod1),format='(4f8.4)'
print,'avdmod2'
print,transpose(avdmod2),format='(4f8.4)'
print,'aveps : ',aveps
print,'aveps1: ',aveps1
print,'aveps2: ',aveps2

;for j=0,nslit-1 do begin
;   dmod1[1:3,*,j]=2.*eps1[4,j]*dmod1[1:3,*,j]/(eps1[4,j]+eps2[4,j])
;   dmod2[1:3,*,j]=2.*eps2[4,j]*dmod2[1:3,*,j]/(eps1[4,j]+eps2[4,j])
;endfor

return
end
