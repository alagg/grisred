pro get_beams,im_in,beam1,beam2,data=data,$
   corr1=corr1,corr2=corr2,lambda=lambda,align=align

if(keyword_set(align) eq 0) then align=0

lamref=[10830.,15650.]
dxref=[507.5,513.5]
coef=poly_fit(lamref,dxref,1) 
ddx=poly(lambda,coef)
;dpx=ddx-fix(ddx)        ; for subpixel aligment between the two beams
dpx=0.                   ; alignment with precion of 1 pixel

im=im_in
tam=size(im)

if(tam(0) eq 2) then begin
   ny=round(data(3))-round(data(2))+1
   beam1=fltarr(tam(1),ny)
   beam2=fltarr(tam(1),ny)
   for j=0,tam(1)-1 do begin
      ystart=data[2]
      yend=data[3]
      y=ystart+findgen(ny)*(yend-ystart)/(ny-1.)
      beam1(j,*)=interpolate(im[j,*],y)
   endfor
   for j=0,tam(1)-1 do begin
      ystart=data[6]
      yend=data[7]
      y=ystart+findgen(ny)*(yend-ystart)/(ny-1.)
      dum=(1-dpx)*im[j,*]+dpx*im[(j+1)<(tam[1]-1),*]
      beam2(j,*)=interpolate(dum,y)
   endfor
   if(align ne 0) then begin
      x=findgen(tam[1])
      for j=0,ny-1 do begin
         beam1[*,j]=interpolate(beam1[*,j],x-corr1(j))
;         beam1[*,j]=desp1d(beam1[*,j],corr1(j))
      endfor
      for j=0,ny-1 do begin
         beam2[*,j]=interpolate(beam2[*,j],x-corr2(j))
;         beam2[*,j]=desp1d(beam2[*,j],corr2(j))
       endfor
   endif
;   beam1=beam1(round(data(0)):round(data(1)),*)
;   beam2=beam2(round(data(0)):round(data(1)),*)
endif else if(tam(0) eq 3) then begin
   ny=round(data(3))-round(data(2))+1
   beam1=fltarr(tam(1),tam(2),ny)
   beam2=fltarr(tam(1),tam(2),ny)
   for k=0,tam(1)-1 do begin
      for j=0,tam(2)-1 do begin
         ystart=data[2]
         yend=data[3]
         y=ystart+findgen(ny)*(yend-ystart)/(ny-1.)
         beam1(k,j,*)=interpolate(im[k,j,*],y)
      endfor
      for j=0,tam(2)-1 do begin
         ystart=data[6]
         yend=data[7]
         y=ystart+findgen(ny)*(yend-ystart)/(ny-1.)
         dum=(1-dpx)*im[k,j,*]+dpx*im[k,(j+1)<(tam[2]-1),*]
         beam2(k,j,*)=interpolate(dum,y)
      endfor
      if(align ne 0) then begin
         x=findgen(tam[2])
         for j=0,ny-1 do begin
            beam1[k,*,j]=interpolate(beam1[k,*,j],x-corr1(j))
;            beam1[k,*,j]=desp1d(reform(beam1[k,*,j]),corr1(j))
         endfor
         for j=0,ny-1 do begin
            beam2[k,*,j]=interpolate(beam2[k,*,j],x-corr2(j))
;            beam2[k,*,j]=desp1d(reform(beam2[k,*,j]),corr2(j))
         endfor
      endif
   endfor
;   beam1=beam1(*,round(data(0)):round(data(1)),*)
;   beam2=beam2(*,round(data(0)):round(data(1)),*)
endif else begin
   print,'THIS ROUTINE IS ONLY VALID FOR IMAGES WITH 2 OR 3 DIMENSIONS'
   stop
endelse
return
end

