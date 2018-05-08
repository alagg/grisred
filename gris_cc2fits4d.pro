pro gris_cc2fits4d,ccmask,outdir=outdir,noflip=noflip

  if n_params() eq 0 then help=1b else begin
  ccfiles=file_search(ccmask,count=cnt)
  
  help=0b
  if cnt eq 0 then help=1b else begin
  ccpos=strlen(ccfiles)-strpos(ccfiles,'cc',/reverse_search)
  if total(ccpos)/2 ne cnt then help=1b
  endelse
  endelse
  
  if help then begin
  print,'Creating 4D-Fits files from GRIS cc-files'
  print,'Usage: The argument must be a file mask to find all relevant cc-files, e.g.:'
  print,'  gris_cc2fits4d,''level1/14jun17.001*cc''[,outdir=''~/tmp/'']'
  print,'This will create 4D-file 14jun17.001.4d.fits'
  return
endif

  print,'cc-files found: ',cnt
  print,ccfiles
  
  fout=strmid(ccfiles[0],0,strpos(ccfiles[0],'.',/reverse_search)+4)+'.4d.fits'
  froot=(reverse(str_sep(fout,'/')))[0]
  if keyword_set(outdir) then fout=outdir+froot
  print,'Output file: ',fout
  ;check write permission
  openw,unit,/get_lun,fout,error=error,/delete
  if error eq 0 then free_lun,unit else begin
  print,'Cannot write '+fout
  print,'define output directory with write permission with keyword outdir=''./'''
  retall
  endelse
  
  
  print,'Reading cc-files:'
  for i=0,cnt-1 do begin
  print,'reading ',ccfiles[i]
  ;    f=rfits_im(ccfiles[i],1,dd,hdr,nrhdr)>0
  stri=string(i,format='(i2.2)')
  dvar='data'+stri
  hvar='hdr'+stri
  exstr=dvar+'=readfits('''+ccfiles[i]+''','+hvar+') & sz=size('+dvar+',/dim)'
  dummy=execute(exstr)
  if i eq 0 then szcc=sz else szcc[2]+=sz[2]
  endfor
  nx=szcc[2]/4
  ny=szcc[1]
  nwl=szcc[0]
  nstk=4
  mkhdr,hdr4d,4,[nwl,nstk,nx,ny],/extend
  sxaddpar,hdr4d,'CTYPE1','NWL','number of wavelength points'
  sxaddpar,hdr4d,'CTYPE2','NSTK','number of Stokes parameters'
  sxaddpar,hdr4d,'CTYPE3','NX','number of scan positions (x-axis)'
  sxaddpar,hdr4d,'CTYPE4','NY','number of pixels along the slit'
  sxaddpar,hdr4d,'STOKES','IQUV','order of Stokes vector'
  hdr4d=hdr4d[0:(where(strmid(hdr4d,0,3) eq 'END'))[0]-1]
  
  dat4d=fltarr(szcc[0],szcc[1],szcc[2])
  off=0
  
  for i=0,cnt-1 do begin
  stri=string(i,format='(i2.2)')
  dvar='data'+stri
  hvar='hdr'+stri
  exstr='data='+dvar+' & hdr='+hvar
  dummy=execute(exstr)
  sz=size(data,/dim)
  dat4d[*,*,off:off+sz[2]-1]=data
  off+=sz[2]
  ;construct header
  utpos=where(strmid(hdr,0,8) eq 'UT      ')
  dut=utpos[1]-utpos[0]
  if i eq 0 then begin
  naxpos=max(where(strmid(hdr,0,5) eq 'NAXIS'))
  hdr4d=[hdr4d,hdr[naxpos+1:max(utpos)+dut-1]]
  endif else if i lt cnt-1 then begin
  hdr4d=[hdr4d,hdr[min(utpos):max(utpos)+dut]]
  endif else begin
  cmtpos=where(strmid(hdr,0,8) eq 'COMMENT ')
  if cmtpos[0] eq -1 then cmtpos=n_elements(hdr)-1
  hdr4d=[hdr4d,hdr[min(utpos):min(cmtpos)]]
  endelse
  endfor
  
  ;re-order array: NWL, NSTOKES, NX, NY
  print,'rearranging cube: nwl,nstk,nx,ny = ',nwl,nstk,nx,ny
  cube=fltarr(nwl,nstk,nx,ny)
  for i=0,nstk-1 do $
  cube[*,i,*,*]=transpose(dat4d[*,*,indgen(nx)*4+i],[0,2,1])
  
  
  ;write out WL cal
  ff1wloff=double(sxpar(hdr4d,'FF1WLOFF'))
  ff1wldsp=double(sxpar(hdr4d,'FF1WLDSP'))
  ff2wloff=double(sxpar(hdr4d,'FF2WLOFF'))
  ff2wldsp=double(sxpar(hdr4d,'FF2WLDSP'))
  ffwl=0b
  if ff1wloff ge 1. then begin
  wloff=ff1wloff
  wldisp=ff1wldsp
  ffwl=1b
  endif
  if ff2wloff ge 1. then begin
  wloff=(ff1wloff+ff2wloff)/2.
  wldisp=(ff1wldsp+ff2wldsp)/2.
  ffwl=1b
  endif
  if ffwl eq 1 then begin
  wl_vec=dindgen(nwl)*wldisp+wloff
  endif
  
  
  ;write out continuum image
  icont=fltarr(nx,ny)
  for ix=0,nx-1 do for iy=0,ny-1 do begin
  icont(ix,iy)=get_cont(cube[*,0,ix,iy])
  endfor
  ;get_histocont
  maxic=max(smooth(icont,5))
  iin=where(icont ge 0.5*maxic)
  hist=histogram(icont[iin],nbin=100,loc=x)
  dummy=max(hist,imax)
  icimg=x[imax]
  ;determine IC_HRSA
  tp=reform(sqrt(total(cube[*,1,*,*]^2+cube[*,2,*,*]^2+cube[*,3,*,*]^2,1))/nwl)
  wd=50<(min([nx,ny])/10)
  medtp=median(tp,wd)           ;boxcar to find most quiet region
  medall=medtp*0+!values.f_nan
  medall[wd:nx-wd-1,wd:ny-wd-1]=medtp[wd:nx-wd-1,wd:ny-wd-1]
  icmed=median(icont,wd)
  iqs=where(icmed ge 0.9*max(icmed))
  dummy=min(medall[iqs],/nan,imin)
  ixy=array_indices(medall,iqs[imin])
  icont_hsra=(median(icont,wd))[ixy[0],ixy[1]]
  
    ; added flip keyword in order to have the maps with the same orientation as 
    ; SDO - Sebastian Castellanos Duran - May 2018
    flipped = 0
    if ~KEYWORD_SET(noflip) then begin
    if double(sxpar(hdr4d,'STEPANGL')) lt 90 then begin
        print,'Flipping data'
        ss = size(cube)
        for ii=0,ss[1]-1 do begin
          for jj=0,ss[2]-1 do begin
            cube[ii,jj,*,*] = reverse(reform(cube[ii,jj,*,*]))
          endfor
        endfor
        icont = reverse(icont)
        flipped = 1
      endif
    endif
    sxaddpar,hdr4d,'DATAFLIP',flipped,' data is flipped? yes:1, no:0',after='STEPANGL'
    sxaddpar,headeric,'DATAFLIP',flipped,' data is flipped? yes:1, no:0',after='STEPANGL'
    
  
  ;write 4d cube
  writefits,fout,cube,hdr4d
  
  ;write wl-extension
  if n_elements(wl_vec) ne 0 then begin
  mkhdr,headerwl,5,nwl,/image
  sxaddpar,headerwl,'WLVEC','','wavelength vector'
  sxaddpar,headerwl,'WLDISP',wldisp,'wavelength dispersion'
  sxaddpar,headerwl,'WLOFF',wloff,'wavelength offset'
  print,'WLDISP',wldisp,' / wavelength dispersion'
  print,'WLOFF',wloff ,' / wavelength offset'
  writefits,fout,wl_vec,headerwl,append=1
  endif else begin
  message,/cont,'No WL-info found. Please reduce the data using a newer' + $
  ' version of grisred before creating the 4D fits file.'
  endelse
  
  ;write continuum extension
  if n_elements(wl_vec) ne 0 then begin
  mkhdr,headeric,4,[nx,ny],/image
  sxaddpar,headeric,'ICONTIMG',icimg,'continuum image'
  sxaddpar,headeric,'IC_HSRA',icont_hsra,'average QS continuum level'
  print,'ICONTIMG',icimg,' / average continuum level'
  print,'IC_HSRA',icont_hsra,' / average QS continuum level'
  writefits,fout,icont,headeric,append=1
  endif else begin
  message,/cont,'No continuum image found. Please reduce the data using a newer' + $
  ' version of grisred before creating the 4D fits file.'
  endelse
  
  print,'Wrote 4D cube to '+fout
  
  end