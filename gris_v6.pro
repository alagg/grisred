pro gris_v6,map_in_base,fileff,filecalib,display=display,$
teles=teles,bitpix=bitpix,desp=desp,factor=factor,mast=mast,theta=theta,$
filt=filt,data=data,option=option,step=step,old=old,getdata=getdata,$
lambda=lambda,order=order,noxtalk=noxtalk,filebad=filebad, $
pzero=pzero,xtau=xtau,dust=dust,$
maxshift=maxshift,checksync=checksync,$	ang_PCU2_az=ang_PCU2_az,$
rotator_offset=rotator_offset, rotator_track=rotator_track,stray=stray,$
fts=fts, show=show,denoise_fft=denoise_fft

ang_PCU2_az=-14.12
ang_derot_PCU2=-47.28
ang_derot_az=ang_derot_PCU2+ang_PCU2_az

if(keyword_set(display) eq 0) then display=0
if(keyword_set(teles) eq 0) then teles=0
if(keyword_set(desp) eq 0) then desp=0
if(keyword_set(factor) eq 0) then factor=1.
if(keyword_set(teles) ne 0 and n_elements(theta) eq 0) then theta = 0
if(keyword_set(filt) eq 0) then filt=0
if(keyword_set(option) eq 0) then option=0
if(keyword_set(step) eq 0) then step=1
if(keyword_set(old) eq 0) then old=0
if(keyword_set(getdata) eq 0) then getdata=0
if(keyword_set(lambda) eq 0) then begin
print,'keyword lambda (in A) is required'
return
endif
if n_elements(fts) eq 0 then begin
print,'continuum will be fitted using FTS spectrum.'
fts=1
endif
if(keyword_set(order) eq 0) and fts then begin
order=grating_angle(lambda)
endif
if(keyword_set(denoise_fft) eq 0) then denoise_fft=0
if(keyword_set(pzero) eq 0) then pzero=53.6 ;41.8	;53.6
if(keyword_set(noxtalk) eq 0) then noxtalk=0
if(keyword_set(dust) eq 0) then dust=0
if(keyword_set(maxshift) eq 0) then maxshift=5
if(keyword_set(checksync) eq 0) then checksync=0
if(keyword_set(stray) eq 0) then stray=0
if(keyword_set(rotator_track) eq 1) then rotator_track=1 else rotator_track=0
if(keyword_set(rotator_offset) eq 0) then rotator_offset=0
ang_derot_az=ang_derot_az+rotator_offset

;get git revision
gitri=routine_info('gris_v6',/source)
gitfi=file_info(file_dirname(gitri.path)+'/grisred.version')
if gitfi.exists then begin
gitrev=strarr(2)
openr,unit,/get_lun,gitfi.name & readf,unit,gitrev & free_lun,unit
endif else gitrev=['n/a','n/a n/a']
print,'GIT-revision: ',gitrev

;use file_search (A. Lagg, May 05)
files=file_search('./level0/'+map_in_base+'*',count=cnt)
if (cnt eq 0) then begin
files=file_search(map_in_base+'*',count=cnt)
pathin="./"
endif else begin
dirout=file_search('level1',count=cntdir)
pathin='./level0/'
if(cntdir eq 0) then spawn,'mkdir ./level1'
endelse

if(cnt eq 0) then begin
message,'No TIP maps found.'
return
endif

iraw=where(strpos(files,'c',/reverse_search) ne strlen(files)-1 and $
           strpos(files,'m',/reverse_search) ne strlen(files)-1)
if iraw(0) eq -1 then begin
message,'No TIP raw data files found.'
return
endif

files=files(iraw)
files=files(sort(files))
nmap_in=n_elements(files)
message,/cont,'Found raw data files: '+strjoin(files,', ')

dum=rfits_im(files[0],1,dd,hdr)
date=param_fits(hdr,'DATE-OBS=',delimiter='-',vartype=1)
year=date[0]

if(year le 2016) then begin
print,'2016'
parameter_mirrors_2016,lambda,x4,tau4,x567,tau567,x8910,tau8910
endif else begin
print,'2017'
parameter_mirrors_2017,lambda,x4,tau4,x567,tau567,x8910,tau8910
endelse

;print,'x4,tau4,x567,tau567,x8910,tau8910 = ',x4,tau4,x567,tau567,x8910,tau8910

if(keyword_set(xtau) eq 0) then begin
xtau=[x4,tau4,x567,tau567,x8910,tau8910]
endif else begin
if(n_elements(xtau) eq 2) then begin
xtau=[xtau,x567,tau567,x8910,tau8910]
endif else if(n_elements(xtau) eq 4) then begin
xtau=[xtau,x8910,tau8910]
endif else begin
xtau=-999
endelse
endelse

if(getdata eq 1) then begin
if(n_elements(data) ne 8) then begin
dcff=rfits_im(pathin+fileff(0),1,/badpix)>0
for j=2,8 do dcff=dcff+rfits_im(pathin+fileff(0),j,/badpix)>0
dcff=dcff/8.
im1=(rfits_im(pathin+fileff(0),9,desp=desp,/badpix)>0)-dcff
im2=(rfits_im(pathin+fileff(0),10,desp=desp,/badpix)>0)-dcff
im3=(rfits_im(pathin+fileff(0),11,desp=desp,/badpix)>0)-dcff
im4=(rfits_im(pathin+fileff(0),12,desp=desp,/badpix)>0)-dcff
data=get_limits(median((im1+im2+im3+im4)/4.,3),lambda)

      if(n_elements(fileff) eq 2) then begin
      dcff=rfits_im(pathin+fileff(1),1,/badpix)>0
      for j=2,8 do dcff=dcff+rfits_im(pathin+fileff(0),j,/badpix)>0
      dcff=dcff/8.
      im1=(rfits_im(pathin+fileff(1),9,desp=desp,/badpix)>0)-dcff
      im2=(rfits_im(pathin+fileff(1),10,desp=desp,/badpix)>0)-dcff
      im3=(rfits_im(pathin+fileff(1),11,desp=desp,/badpix)>0)-dcff
      im4=(rfits_im(pathin+fileff(1),12,desp=desp,/badpix)>0)-dcff
      data2=get_limits(median((im1+im2+im3+im4)/4.,3),lambda)
      
      data[2]=data[2] > data2[2]
      data[3]=data[3] < data2[3]
      data[6]=data[6] > data2[6]
      data[7]=data[7] < data2[6]
      endif
      
      if((data(3)-data(2)) ne (data(7)-data(6))) then begin
      print,'**********************************************************'
      print,'PROBLEMS IN THE AUTOMATIC DETERMINATION OF THE BEAM LIMITS'
      print,'Run the manual routine (acum2iquv4)'
      print,'**********************************************************'
      endif  
      endif
      return
    endif
    
    print,' '
    print,'FLATFIELD'
    print,' '
    pos=0
    pendiente=0
    
    nff=n_elements(fileff)
    
    get_flat,pathin+fileff(0),ff1,data=data,time=timeff1,corr1=dd1,corr2=dd2,$
    lambda=lambda,order=order,ffnorm=ffnorm1,periodfr=periodfr,show=show, $
    fts=fts,fit=ftsfit1
    
    i1=data(0)
    i2=data(1)
    j1=data(2)
    j2=data(3)
    i3=data(4)
    i4=data(5)
    j3=data(6)
    j4=data(7)
    
    nslit=j4-j3+1
    
    ff1=ff1>0.1
    timeff1=mean(timeff1)
    if(nff eq 1) then begin
    ff2=ff1
    timeff2=timeff1
    ffnorm2=ffnorm1
    endif else if(nff eq 2 and fileff(0) eq fileff(1)) then begin
    ff2=ff1
    timeff2=timeff1
    ffnorm2=ffnorm1
    endif else begin   
    get_flat,pathin+fileff(1),ff2,data=data,time=timeff2,corr1=dd1,corr2=dd2,$
    lambda=lambda,order=order,ffnorm=ffnorm2,show=show,fts=fts, $
    fit=ftsfit2
    ff2=ff2>0.1
    timeff2=mean(timeff2)
    endelse   
    
    if(dust eq 1) then begin
    ff1tot=total(ff1,1)/4.
    ff2tot=total(ff2,1)/4.
    dust1_down=total(ff1tot[data[0]:data[1],data[2]:data[3]],1)/(data[1]-data[0]+1)
    dust1_up=total(ff1tot[data[4]:data[5],data[6]:data[7]],1)/(data[5]-data[4]+1)
    dust2_down=total(ff2tot[data[0]:data[1],data[2]:data[3]],1)/(data[1]-data[0]+1)
    dust2_up=total(ff2tot[data[4]:data[5],data[6]:data[7]],1)/(data[5]-data[4]+1)
    
    dust1_down=dust1_down/continuum(dust1_down,5)
    dust1_up=dust1_up/continuum(dust1_up,5)
    dust2_down=dust2_down/continuum(dust2_down,5)
    dust2_up=dust2_up/continuum(dust2_up,5)
    for k=0,3 do begin
    for j=data[0],data[1] do begin
    ff1[k,j,data[2]:data[3]]=ff1[k,j,data[2]:data[3]]/dust1_down
    ff1[k,j,data[6]:data[7]]=ff1[k,j,data[6]:data[7]]/dust1_up
    ff2[k,j,data[2]:data[3]]=ff2[k,j,data[2]:data[3]]/dust2_down
    ff2[k,j,data[6]:data[7]]=ff2[k,j,data[6]:data[7]]/dust2_up
    endfor
    endfor
    endif
    
    ffnorm=(ffnorm1+ffnorm2)/2.
    
    ;ff1=ff1*0.1     ; to have the final files with a factor 10 more counnts
    ;ff2=ff2*0.1
    
    print,' '
    print,'CALIBRATION'
    print,' '
    
    ; ask if derotator is inserted
    
    filecalib=filecalib[0]
    dumcal=rfits_im(pathin+filecalib[0],1,ddcal,hdrcal,nrhdr)
    
    rotcode=param_fits(hdrcal,'ROTCODE =',vartype=1)
    if(min(abs(rotcode)) eq 0) then begin
    rotangle=param_fits(hdrcal,'ROTANGLE=',vartype=3)
    derotator_angle=rotangle[0]+ang_derot_az
    print,'CALIBRATION 1: derotator IN ',rotangle
    endif else begin
    derotator_angle=-999
    print,'CALIBRATION 1: derotator OUT'
    endelse
    ;pzero=0
    ;rzero=0
    if(lambda eq 0) then begin
    lambda=param_fits(hdr,'WAVELENG=',vartype=1)*10. 
    print,'ATTENTION: wavelength taken from header = ',lambda, ' Angstroem'
    endif
    
    ;delta=773.*1.e4/lambda-401.8
    ;delta=786.5*1.e4/lambda-408.3	; NEW APRIL 2002
    ;delta=789.5*1.e4/lambda-412.4	; NEW MAY 2002
    
    ;calib_slit4,filecalib,ff1,dmod,data=data,lambda=lambda
    nfcal=n_elements(filecalib)
    if(keyword_set(filebad) eq 0) then zbad=-999 else restore,filebad[0]
    if(nff eq 1) then begin
    get_calib_v5,pathin+filecalib[0],dmod1a,dmod2a,data=data,lambda=lambda,$
    zbad=zbad,timecal=timecal1,date=date1,pzero=pzero, $
    flat1=ff1,timeflat1=timeff1
    endif else if(nff eq 2 and fileff(0) eq fileff(1)) then begin 
    get_calib_v5,pathin+filecalib[0],dmod1a,dmod2a,data=data,lambda=lambda,$
    zbad=zbad,timecal=timecal1,date=date1,pzero=pzero, $
    flat1=ff1,timeflat1=timeff1
    endif else begin   
    get_calib_v5,pathin+filecalib[0],dmod1a,dmod2a,data=data,lambda=lambda,$
    zbad=zbad,timecal=timecal1,date=date1,pzero=pzero, $
    flat1=ff1,timeflat1=timeff1,flat2=ff2,timeflat2=timeff2
    endelse
    
    timecal1=timecal1[0]
    elaz=get_elaz(date1[2],date1[1],date1[0],timecal1,decsun,p0)
    tel_ref1=gregormod_v6(elaz[0],elaz[1],derotator_angle,decsun,p0,/POLCALIB,xtau=xtau,$
                          lambda=lambda/10.)
    if(nfcal eq 1) then begin
    dmod1b=dmod1a
    dmod2b=dmod2a
    timecal2=timecal1
    tel_ref2=tel_ref1
    endif else if(nff eq 2 and fileff(0) eq fileff(1)) then begin
    dmod1b=dmod1a
    dmod2b=dmod2a
    timecal2=timecal1
    tel_ref2=tel_ref1
    endif else begin   
    if(keyword_set(filebad) eq 0) then zbad=-999 else restore,filebad[1]
    
    ; ask if derotator is inserted
    dumcal=rfits_im(pathin+filecalib[1],1,ddcal,hdrcal,nrhdr)
    
    rotcode=param_fits(hdrcal,'ROTCODE =',vartype=1)
    if(min(abs(rotcode)) eq 0) then begin
    rotangle=param_fits(hdr,'ROTANGLE=',vartype=3)
    derotator_angle=rotangle[0]+ang_derot_az
    print,'CALIBRATION 2: derotator IN ',rotangle
    endif else begin
    derotator_angle=-999.
    print,'CALIBRATION 2: derotator OUT '
    endelse
    
    if(nff eq 1) then begin
    get_calib_v5,pathin+filecalib[1],dmod1b,dmod2b,data=data,lambda=lambda,$
    zbad=zbad,timecal=timecal2,date=date2,pzero=pzero, $
    flat1=ff1,timeflat1=timeff1
    endif else if(nff eq 2 and fileff(0) eq fileff(1)) then begin 
    get_calib_v5,pathin+filecalib[1],dmod1b,dmod2b,data=data,lambda=lambda,$
    zbad=zbad,timecal=timecal2,date=date2,pzero=pzero, $
    flat1=ff1,timeflat1=timeff1
    endif else begin   
    get_calib_v5,pathin+filecalib[1],dmod1b,dmod2b,data=data,lambda=lambda,$
    zbad=zbad,timecal=timecal2,date=date2,pzero=pzero, $
    flat1=ff1,timeflat1=timeff1,flat2=ff2,timeflat2=timeff2
    endelse
    
    timecal2=timecal2[0]
    elaz=get_elaz(date1[2],date1[1],date1[0],timecal2,decsun,p0)
    
    tel_ref2=gregormod_v6(elaz[0],elaz[1],derotator_angle,decsun,p0,/POLCALIB,xtau=xtau, $
                          lambda=lambda/10.)
    endelse  
    
    if(checksync eq 1) then begin
    dmod1am=total(dmod1a,3)/(data[3]-data[2]+1)
    dmod1bm=total(dmod1b,3)/(data[3]-data[2]+1)
    dmod2am=total(dmod2a,3)/(data[7]-data[6]+1)
    dmod2bm=total(dmod2b,3)/(data[7]-data[6]+1)
    endif
    
    print,' '
    print,'FLATFIELDING + DEMODULATING + MERGING BEAMS'
    print,' '
    
    nrhdr=n_elements(hdr)
    bitpix=32
    
    ; ask if derotator is inserted
    dum_start=rfits_im(files[0],1,dd_start,hdr_start,nrhdr_start)
    rotcode=param_fits(hdr_start,'ROTCODE =',vartype=1)
    ;rotcode=0		; 
    if(min(abs(rotcode)) eq 0) then begin
    
    rotangle_start=param_fits(hdr_start,'ROTANGLE=',vartype=3)
    rotangle_start=rotangle_start[0]
    derotator_angle_start=rotangle_start+ang_derot_az
    time_start=param_fits(hdr_start,'UT      =',delimiter=':',vartype=1) 
    time_start=time_start(*,0)+(time_start(*,1)+time_start(*,2)/60.)/60.
    time_start=time_start[0]
    
    dum_final=rfits_im(files[nmap_in-1],1,dd_final,hdr_final,nrhdr_final)
    rotangle_end=param_fits(hdr_final,'ROTANGLE=',vartype=3)
    rotangle_end=(reverse(rotangle_end))[0]
    
    if((rotangle_end eq 0 or rotangle_end lt rotangle_start) and rotator_track eq 0) then begin
    rotangle_end=rotangle_start
    print,' '
    print,'******************'
    print,'NOTE: DE-ROTATOR IN FIXED POSITION. PLEASE CHECK!!!'
    print,'******************'
    print,' '
    endif
    
    derotator_angle_end=rotangle_end+ang_derot_az
    time_end=param_fits(hdr_final,'UT      =',delimiter=':',vartype=1) 
    time_end=time_end(*,0)+(time_end(*,1)+time_end(*,2)/60.)/60.
    time_end=(reverse(time_end))[0]
    
    print,'MAPS: derotator IN ',rotangle_start,':',rotangle_end
    
    endif else begin
    derotator_angle_start=-999.
    derotator_angle_end=-999.
    
    print,'MAPS: derotator OUT '
    
    endelse 
    
    
    
    for jj=0,nmap_in-1 do begin
    map_in=files(jj)
    print,map_in
    map_hdr = map_in+'cc'
    dum=rfits_im(map_in,1,dd,hdr,nrhdr,/badpix)
    tam=size(dum)
    close,1
    
    if(jj eq 0) then begin 
    offset=8 
    dc=dum
    for j=2,8 do begin
    dc=dc + rfits_im(map_in,j,desp=desp,/badpix)>0
    endfor
    dc=dc/8.        
    z=where(abs(dc-mean(dc)) gt 5*std(dc))
    if(z(0) ne -1) then dc(z)=mean(dc)
    
    naccum=param_fits(hdr,'ACCUMULA=',vartype=1)
    factor=ffnorm*naccum
    endif else begin
    offset=0
    endelse
    npos=(dd.naxis3-offset)/4 
    
    get_lun,unit
    openr,unit,map_in
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
    
    im_a=fltarr(4,tam(1),tam(2))
    im_b=fltarr(4,tam(1),tam(2))
    im2_a=fltarr(4,i2-i1+1,j2-j1+1)
    im2_b=fltarr(4,i2-i1+1,j2-j1+1)
    im2=fltarr(4,i2-i1+1,j2-j1+1)
    fr=fltarr(4,i2-i1+1,j2-j1+1)
    toti=fltarr(npos,j2-j1+1)
    totq=fltarr(npos,j2-j1+1)
    totu=fltarr(npos,j2-j1+1)
    totv=fltarr(npos,j2-j1+1)
    
    if(display eq 0) then begin
    
    map_out=map_in
    
    pos=strpos(map_out,'level0') 
    if(pos ne -1) then strput,map_out,string(format='(a6)',"level1"),pos 
    
    get_lun,unit_out
    openw,unit_out,map_out+'cc'
    
    ; The new keywords are added at the end of the code in order to 
    ; solve the invalid character bug.
    ; Sebastian Castellanos Duran - May 2018
    
    for j=0L,nrhdr-1 do begin
    header=hdr(j)
    
    pos=strpos(hdr(j),'BITPIX  =')
    if(pos ne -1) then strput,header,string(format='(i20)',fix(bitpix)),pos+10
    
    pos=strpos(header,'NAXIS1  =')
    if(pos ne -1) then strput,header,string(format='(i20)',i2-i1+1),pos+10
    
    pos=strpos(header,'NAXIS2  =')
    if(pos ne -1) then strput,header,string(format='(i20)',j2-j1+1),pos+10
    
    pos=strpos(header,'NAXIS3  =')
    if(pos ne -1) then strput,header,string(format='(i20)',4*fix(npos/step)),pos+10
    
    pos=strpos(hdr(j),'TELESCOP=')
    if(dd.telescope eq 'SVST' and pos ne -1) then $
    strput,header,string(format='(a8)','SVST_COR'),pos+11   
    
    pos=strpos(header,'WAVELENG=')
    if(pos ne -1) then strput,header,string(format='(i20)',round(lambda/10.)),pos+10
    
    pos=strpos(hdr(j),'DATAVERS=')
    if(pos ne -1) then begin
    strput,header,string(format='(i20)',fix(1)),pos+10
    tit='data version (I->Q,U,V crosstalk corrected)'
    format_tit='(A'+strtrim(strlen(tit),2)+')'
    strput,header,string(format=format_tit,tit),pos+33
    endif
    
    hdr(j)=header
    endfor
    writeu,unit_out,byte(hdr)
    
    if(bitpix eq 8) then begin
    dat_out=assoc(unit_out,bytarr(i2-i1+1,j2-j1+1),long(2880)*nrhdr)
    endif else if(bitpix eq 16) then begin   
    dat_out=assoc(unit_out,intarr(i2-i1+1,j2-j1+1),long(2880)*nrhdr)
    endif else if(bitpix eq 32) then begin   
    dat_out=assoc(unit_out,lonarr(i2-i1+1,j2-j1+1),long(2880)*nrhdr)
    endif
    endif
    
    time=param_fits(hdr,'UT      =',delimiter=':',vartype=1) 
    time=time(*,0)+(time(*,1)+time(*,2)/60.)/60.
    istep=param_fits(hdr,'ISTEP   =',vartype=1)
    if(old eq 0) then begin
    date=param_fits(hdr,'DATE-OBS=',delimiter='-',vartype=1)
    dum=date[0]
    date[0]=date[2]
    date[2]=dum
    endif else begin
    date=param_fits(hdr,'DATE    =',delimiter='/',vartype=1)
    date[2]=date[2]+1900
    endelse
    
    format=['(i2,$)','(i3,$)','(i4,$)','(i5,$)','(i6,$)']
    print,npos
    
    for i=0L,npos-1,step do begin
    
    print,i+1,format=format(fix(alog10(i+1)))
    for j=0,3 do begin
    im_a(j,*,*)=((rfits_im2(datos,dd,4*i+j+offset+1,desp=desp,/badpix)>0)-dc)	
    endfor   
    
    if(teles ne 0 and dd.telescope eq 'VTT') then begin
    mat_tel=vtt(n,k,date[2],date[1],date[0],time(i),theta,mast,lambda=lambda)
    mat_tel=invert(mat_tel)
    endif else if(dd.telescope eq 'GREGOR') then begin
    elaz=get_elaz(date[2],date[1],date[0],time[i],decsun,p0)
    
    if(jj eq 0 and i eq 0) then begin
    el_start=elaz[0]
    az_start=elaz[1]
    parall_start=elaz[2]
    p0_start=p0
    endif
    
    if(min(abs(rotcode)) eq 0) then begin
    if(derotator_angle_end ne derotator_angle_start or rotator_track eq 1) then begin
    ;               cte=(time(i)-time_start)/(time_end-time_start)
    ;               derotator_angle=(1.-cte)*derotator_angle_start+cte*derotator_angle_end
    derotator_angle=derotator_angle_start+(elaz[0]-el_start)/2.+(elaz[1]-az_start)/2.+ $
    (elaz[2]-parall_start)/2.
    ;               print,' DEROTATOR_ANGLE= ',derotator_angle-ang_derot_az
    endif else begin
    derotator_angle=derotator_angle_start
    endelse
    endif else begin
    derotator_angle=-999.
    endelse
    
    ;         if(min(abs(rotcode)) eq 0) then begin
    ;            delta_elev=elaz[0]-el_start
    ;            delta_az=elaz[1]-az_start
    ;            delta_parall=elaz[2]-parall_start
    ;            delta_p0=p0-p0_start
    ;            derotator_angle=derotatorangle_start+(delta_parall+delta_p0-delta_elev-delta_az)/2.
    ;         endif else begin
    ;            derotator=-999.
    ;         endelse
    
    mat_tel=gregormod_v6(elaz[0],elaz[1],derotator_angle,decsun,p0,xtau=xtau,lambda=lambda/10.)
    mat_tel=invert(mat_tel)
    
    endif else if (teles ne 0 and dd.telescope eq 'SVST') then begin
    x=[0.969392,0.969392,1.01156]
    tau=[157.036,157.036,145.447]
    windows=[1.54097,7.94586,306.629,0.816829]
    mat_tel=svst_xtau(x,tau,65.2,date[2],date[1],date[0],time(i),windows=windows)  
    mat_tel=invert(mat_tel)
    endif else begin
    mat_tel=fltarr(4,4)
    for j=0,3 do mat_tel(j,j)=1.
    endelse      
    
    if(timecal2 ne timecal1) then begin
    fact_time=(timecal2-time(i))/(timecal2-timecal1)
    endif else begin
    fact_time=1
    endelse
    
    if(checksync eq 1) then begin
    ffsync=fact_time*ff1+(1-fact_time)*ff2
    dum1=mat_tel#tel_ref1#dmod1am#total(im_a[*,*,j1:j2]/ffsync[*,*,j1:j2],3)/(j2-j1+1)
    dum2=mat_tel#tel_ref1#dmod1am#total(shift(im_a[*,*,j1:j2],1,0,0)/ffsync[*,*,j1:j2],3)/(j2-j1+1)
    dum3=mat_tel#tel_ref1#dmod1am#total(shift(im_a[*,*,j1:j2],2,0,0)/ffsync[*,*,j1:j2],3)/(j2-j1+1)
    dum4=mat_tel#tel_ref1#dmod1am#total(shift(im_a[*,*,j1:j2],3,0,0)/ffsync[*,*,j1:j2],3)/(j2-j1+1)
    dum1m=abs(median(dum1[1:3,*]))
    dum2m=abs(median(dum2[1:3,*]))
    dum3m=abs(median(dum3[1:3,*]))
    dum4m=abs(median(dum4[1:3,*]))
    pos=sort([dum1m,dum2m,dum3m,dum4m])
    im_a=shift(im_a,pos[0],0,0)
    endif
    
    im_b=im_a/ff2
    im_a=im_a/ff1
    
    if(dust eq 1) then begin
    espa_down=total(im_a[*,data[0]:data[1],data[2]:data[3]],1)/4.           
    espa_down=total(espa_down,1)/(data[1]-data[0])
    espa_up=total(im_a[*,data[4]:data[5],data[6]:data[7]],1)/4.
    espa_up=total(espa_up,1)/(data[5]-data[4])                    
    
    espb_down=total(im_b[*,data[0]:data[1],data[2]:data[3]],1)/4.           
    espb_down=total(espb_down,1)/(data[1]-data[0])
    espb_up=total(im_b[*,data[4]:data[5],data[6]:data[7]],1)/4.
    espb_up=total(espb_up,1)/(data[5]-data[4])                    
    
    ncs=2*maxshift+1
    ca_down=fltarr(ncs)
    for j=0,ncs-1 do ca_down[j]=correlate(espa_down,shift(dust1_down,j-maxshift))
    cmax=where(ca_down[1:ncs-2] eq max(ca_down[1:ncs-2]))+1
    cmax=cmax[0]
    coef=poly_fit([cmax-1,cmax,cmax+1],ca_down[cmax-1:cmax+1],2,cfit)
    shifta_down=-coef[1]/2./coef[2]-maxshift
    
    ca_up=fltarr(ncs)
    for j=0,ncs-1 do ca_up[j]=correlate(espa_up,shift(dust1_up,j-maxshift))
    cmax=where(ca_up[1:ncs-2] eq max(ca_up[1:ncs-2]))+1
    cmax=cmax[0]
    coef=poly_fit([cmax-1,cmax,cmax+1],ca_up[cmax-1:cmax+1],2,cfit)
    shifta_up=-coef[1]/2./coef[2]-maxshift
    
    cb_down=fltarr(ncs)
    for j=0,ncs-1 do cb_down[j]=correlate(espb_down,shift(dust2_down,j-maxshift))
    cmax=where(cb_down[1:ncs-2] eq max(cb_down[1:ncs-2]))+1
    cmax=cmax[0]
    coef=poly_fit([cmax-1,cmax,cmax+1],cb_down[cmax-1:cmax+1],2,cfit)
    shiftb_down=-coef[1]/2./coef[2]-maxshift
    
    cb_up=fltarr(ncs)
    for j=0,ncs-1 do cb_up[j]=correlate(espb_up,shift(dust2_up,j-maxshift))
    cmax=where(cb_up[1:ncs-2] eq max(cb_up[1:ncs-2]))+1
    cmax=cmax[0]
    coef=poly_fit([cmax-1,cmax,cmax+1],cb_up[cmax-1:cmax+1],2,cfit)
    shiftb_up=-coef[1]/2./coef[2]-maxshift
    
    for j=0,3 do begin
    for idust=data[0],data[1] do begin
    im_a[j,idust,data[2]:data[3]]=im_a[j,idust,data[2]:data[3]]/ $
    desp1d(dust1_down,shifta_down)
    im_b[j,idust,data[2]:data[3]]=im_b[j,idust,data[2]:data[3]]/ $
    desp1d(dust2_down,shiftb_down)
    endfor
    for idust=data[4],data[5] do begin
    im_a[j,idust,data[6]:data[7]]=im_a[j,idust,data[6]:data[7]]/ $
    desp1d(dust1_up,shifta_up)
    im_b[j,idust,data[6]:data[7]]=im_b[j,idust,data[6]:data[7]]/ $
    desp1d(dust2_up,shiftb_up)
    endfor
    endfor   
    endif
    
    for k=0,nslit-1 do begin
    
    matriz1a=mat_tel#tel_ref1#reform(dmod1a(*,*,k))
    matriz1b=mat_tel#tel_ref2#reform(dmod1b(*,*,k))
    matriz1=fact_time*matriz1a+(1-fact_time)*matriz1b
    
    matriz2a=mat_tel#tel_ref1#reform(dmod2a(*,*,k))
    matriz2b=mat_tel#tel_ref2#reform(dmod2b(*,*,k))
    matriz2=fact_time*matriz2a+(1-fact_time)*matriz2b
    
    im_a(*,*,k+j1)=matriz1#reform(im_a(*,*,k+j1))
    im_a(*,*,k+j3)=matriz2#reform(im_a(*,*,k+j3))
    endfor  
    
    get_beams,im_a,beam1,beam2,data=data,corr1=dd1,corr2=dd2,lambda=lambda,/align
    imdemod1=beam1[*,i1:i2,*]
    imdemod2=beam2[*,i3:i4,*]
    
    im2_a(0,*,*)=(imdemod1(0,*,*)+imdemod2(0,*,*))/2.
    im2_a(1:3,*,*)=(imdemod1(1:3,*,*)-imdemod2(1:3,*,*))/2.
    im2_a[*,*,0]=im2_a[*,*,1]         ; I do not know why row 0 is very noisy
    
    nx=i2-i1+1
    for j=0,j2-j1 do begin
    esp=reform(im2_a[0,*,j])
    for k=-2,2 do esp=esp-franjas(esp,float(nx)/(periodfr+k),1)
    im2_a[0,*,j]=esp
    endfor
    
    if(filt eq 1) then begin
    
    for j=0,j2-j1 do begin
    ;            esp=reform(im2_a(0,*,j))
    ;            im2_a(0,*,j)=esp-franjas(esp,24,4)
    esp=reform(im2_a(1,*,j))
    im2_a(1,*,j)=esp-franjas(esp,71.5,1)
    esp=reform(im2_a(2,*,j))
    im2_a(2,*,j)=esp-franjas(esp,71.5,1)
    ;            esp=reform(im2_a(3,*,j))
    ;            im2_a(3,*,j)=esp-franjas(esp,71.5,1)
    endfor
    endif  
    
    if(timeff1 ne timeff2) then begin
    
    for k=0,nslit-1 do begin
    if(timecal2 ne timecal1) then begin
    fact_time=(timecal2-time(i))/(timecal2-timecal1)
    endif else begin
    fact_time=1
    endelse
    
    matriz1a=mat_tel#tel_ref1#reform(dmod1a(*,*,k))
    matriz1b=mat_tel#tel_ref2#reform(dmod1b(*,*,k))
    matriz1=fact_time*matriz1a+(1-fact_time)*matriz1b
    
    matriz2a=mat_tel#tel_ref1#reform(dmod2a(*,*,k))
    matriz2b=mat_tel#tel_ref2#reform(dmod2b(*,*,k))
    matriz2=fact_time*matriz2a+(1-fact_time)*matriz2b
    
    im_b(*,*,k+j1)=matriz1#reform(im_b(*,*,k+j1))
    im_b(*,*,k+j3)=matriz2#reform(im_b(*,*,k+j3))
    endfor  
    
    get_beams,im_b,beam1,beam2,data=data,corr1=dd1,corr2=dd2,lambda=lambda,/align
    imdemod1=beam1[*,i1:i2,*]
    imdemod2=beam2[*,i3:i4,*]
    
    im2_b(0,*,*)=(imdemod1(0,*,*)+imdemod2(0,*,*))/2.
    im2_b(1:3,*,*)=(imdemod1(1:3,*,*)-imdemod2(1:3,*,*))/2.
    im2_b[*,*,0]=im2_b[*,*,1]         ; I do not know why row 0 is very noisy
    
    nx=i4-i3+1
    for j=0,j2-j1 do begin
    esp=reform(im2_b[0,*,j])
    for k=-2,2 do esp=esp-franjas(esp,float(nx)/(periodfr+k),1)
    im2_b[0,*,j]=esp
    endfor
    
    if(filt eq 1) then begin
    for j=0,j2-j1 do begin
    ;               esp=reform(im2_b(0,*,j))
    ;               im2_b(0,*,j)=esp-franjas(esp,24,4)
    esp=reform(im2_b(1,*,j))
    im2_b(1,*,j)=esp-franjas(esp,71.5,1)
    esp=reform(im2_b(2,*,j))
    im2_b(2,*,j)=esp-franjas(esp,71.5,1)
    ;               esp=reform(im2_b(3,*,j))
    ;               im2_b(3,*,j)=esp-franjas(esp,71.5,1)
    endfor
    endif  
    fact_time=(timeff2-time(i))/(timeff2-timeff1)
    im2=fact_time*im2_a + (1-fact_time)*im2_b
    endif else begin
    im2=im2_a
    endelse    
    
    if(option eq 1 and istep(i) eq 1) then begin
    period=45
    for k=1,3 do begin
    for j=0,j2-j1 do begin
    esp=reform(im2(k,*,j))
    fr(k,*,j)=ajusta_seno(x,esp,period)+ajusta_seno(x,esp,2*period)
    fr(k,*,j)=fr(k,*,j)-mean(fr(k,*,j))
    endfor
    endfor
    endif
    im2=im2-fr
    
    ;added FFT denoising for TIP2 data,
    ;A. Lagg, Oct 2005
    if keyword_set(denoise_fft) then im2=fft_denoise(im2)
    
    for j=0,3 do  im2[j,*,*]=im2[j,*,*]-stray*mean(im2[j,*,*])
    
    if(noxtalk eq 0) then xtalk_v5,im2
    
    ii2=reform(im2(0,*,*))>1
    qq2=reform(im2(1,*,*))/ii2
    uu2=reform(im2(2,*,*))/ii2
    vv2=reform(im2(3,*,*))/ii2
    
    toti(i,*)=total(abs(ii2(1:i2-i1-1,*)),1)/(i2-i1-1)
    totq(i,*)=total(abs(qq2(1:i2-i1-1,*)),1)/(i2-i1-1)
    totu(i,*)=total(abs(uu2(1:i2-i1-1,*)),1)/(i2-i1-1)
    totv(i,*)=total(abs(vv2(1:i2-i1-1,*)),1)/(i2-i1-1)
    
    if(display ne 0) then begin
    tamv2=size(vv2)
    nv2=tamv2(1)
    mv2=tamv2(2)
    ii2=ii2(1:nv2-2,1:mv2-2)
    qq2=qq2(1:nv2-2,1:mv2-2)
    uu2=uu2(1:nv2-2,1:mv2-2)
    vv2=vv2(1:nv2-2,1:mv2-2)
    qq2=qq2-mean(qq2)
    uu2=uu2-mean(uu2)
    vv2=vv2-mean(vv2)
    qq2=qq2/max(abs(qq2))
    uu2=uu2/max(abs(uu2))
    vv2=vv2/max(abs(vv2))
    mm=transpose([transpose(vv2),transpose(uu2),transpose(qq2)])
    mmm=max(mm)
    ii2=mmm*(ii2-min(ii2))/(max(ii2)-min(ii2))
    
    tvwinp,compress([reform(5*im2[1,*,*]),5*reform(im2[2,*,*]),reform(im2[3,*,*])],2)
    endif  else begin
    
    if(!version.arch eq "alpha" or strmid(!version.arch,0,3) eq "x86") then begin   
    if(bitpix eq 8) then begin
    for j=0,3 do dat_out(4*i+j) = byte(im2(j,*,*)/factor)
    endif else if(bitpix eq 16) then begin   
    for j=0,3 do begin
    dum = fix(im2(j,*,*)/factor)
    byteorder,dum
    dat_out(4*i+j) = dum
    endfor   
    endif else if(bitpix eq 32) then begin   
    for j=0,3 do begin
    dum = long(im2(j,*,*)/factor)
    byteorder,dum,/lswap
    dat_out(4*i+j) = dum
    endfor   
    endif
    endif else begin
    if(bitpix eq 8) then begin
    for j=0,3 do dat_out(4*i+j) = byte(im2(j,*,*)/factor)
    endif else if(bitpix eq 16) then begin   
    for j=0,3 do dat_out(4*i+j) = fix(im2(j,*,*)/factor)
    endif else if(bitpix eq 32) then begin   
    for j=0,3 do dat_out(4*i+j) = long(im2(j,*,*)/factor)
    endif
    endelse
    endelse      
    endfor
    
    free_lun,unit
    close,unit
    totl=sqrt(totu*totu+totq*totq)
    totp=sqrt(totl*totl+totv*totv)
    free_lun,unit
    if(display eq 0) then free_lun,unit_out
    save,filename=map_out+'cm',toti,totq,totu,totv,totl,totp,hdr
    print,' '
    print,' '
    
    endfor
    

  
    ; add git revision to the header
    ; patch added to temporary solve the header character bug. 
    ; It just added the new keywords and the "END" on the header.
    ; Sebastian Castellanos Duran - May 2018
    hdrstr='' & data =0.
    
    if strmatch(map_in_base,'*level0/*') eq 1 then begin 
    pos = strpos(map_in_base,'level0/')
    map_in_base=strmid(map_in_base,pos+7)
    endif
    ccfiles=file_search('./level1/'+map_in_base+'*cc',count=nfiles)
    
    for i=0,nfiles-1 do begin
    print,'Writing keywords on the header --> '+ccfiles[i]
    data = readfits(ccfiles[i],hdrstr)
    
    if keyword_set(fts) then nfts=1 else nfts=0
    sxaddpar,hdrstr,'FTSFLAT',nfts,' flatfield calibration with FTS spectrum',before='LC1-1'
    if n_elements(ftsfit1) ne 0 then ftsfit=ftsfit1
    if n_elements(ftsfit2) ne 0 then ftsfit=[ftsfit1,ftsfit2]
    for ii=0,n_elements(ftsfit)-1 do begin
    fs='FF'+strcompress(string(ii+1),/remove_all)
    sxaddpar,hdrstr,fs+'WLOFF',ftsfit[ii].wloff,' WL-offset ',before='LC1-1'
    sxaddpar,hdrstr,fs+'WLDSP',ftsfit[ii].wlbin,' WL-dispersion '+fs,before='LC1-1'
    sxaddpar,hdrstr,fs+'FWHMA',ftsfit[ii].fwhm_a,' spectral FWHM [A] '+fs,before='LC1-1'
    sxaddpar,hdrstr,fs+'FWHMP',ftsfit[ii].fwhm_pix,' spectral FWHM [pix] '+fs,before='LC1-1'
    sxaddpar,hdrstr,fs+'STRAY',ftsfit[ii].stray,' spectral straylight '+fs,before='LC1-1'
    sxaddpar,hdrstr,fs+'NPOLY',ftsfit[ii].npoly,' order of fitted polynomial '+fs,before='LC1-1'
    sxaddpar,hdrstr,fs+'FITNS',ftsfit[ii].fitness,' fitness of PIKAIA fit to FTS '+fs,before='LC1-1'
    endfor
    
    sxaddpar,hdrstr,'DATEREDU',systime(0),' date of data reduction (yyyy-mm-dd)',before='FILENAME'

    gitrev_tmp=string(gitrev[0])
    sxaddpar,hdrstr,'GITREV',gitrev_tmp,' grisred git revision',before='FILENAME'

    gitrev_tmp=(strsplit(gitrev[1],/extract))[1]  
    sxaddpar,hdrstr,'GITREP',gitrev_tmp,' grisred git repository',before='FILENAME'
    
    ;;remove blank spaces
    blank_str=string(bytarr(80)+32B)
    zblank=where(hdrstr eq blank_str,nblank)
    if nblank gt 0 then for j=0,nblank-1 do hdrstr=remove_blankline(hdrstr)
    
    WRITEFITS, ccfiles[i], data, hdrstr
    endfor
    
    return
  end
