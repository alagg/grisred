function getbadpoints,file,data=data

dc=rfits_im(file,1,dd)
for j=2,8 do dc=dc+rfits_im(file,j,/badp)
dc=dc/8.

npos=(dd.naxis3-8)/4

dum=rfits_im(file,9,/badp)>0-dc
dum=dum+(rfits_im(file,10,/badp)>0-dc)
dum=dum+(rfits_im(file,11,/badp)>0-dc)
dum=dum+(rfits_im(file,12,/badp)>0-dc)

if(keyword_set(data) eq 0) then data=limits_tip(dum)
slit=data[3]-data[2]+1
tam=size(dc)

haz1=fltarr(npos,4)
haz2=fltarr(npos,4)
for j=0,npos-1 do begin
   for k=0,3 do begin
      dum=rfits_im(file,4*j+k+9,/badp)-dc
      haz1[j,k]=total(dum[*,data[2]:data[3]])/tam[1]/slit
      haz2[j,k]=total(dum[*,data[6]:data[7]])/tam[1]/slit
   endfor
endfor

nbad=0
loadct,2
xx=indgen(npos)
for j=0,3 do begin
   if(nbad eq 0) then begin
      plot,xx,haz1[*,j],psym=-1
   endif else if(nbad eq 1) then begin
      plot,xx,haz1[*,j],psym=-1 
      oplot,[zbad,zbad],[haz1[zbad,j],haz1[zbad,j]],psym=1,color=80
   endif else begin
      plot,xx,haz1(*,j),psym=-1 & oplot,zbad,haz1[zbad,j],psym=1,color=80
   endelse

   !err=0
   while(!err ne 4) do begin
      print,'click for new point'
      cursor,x,y,3,/down
      if (!err eq 1) then begin
         nbad=nbad+1
         dd=((xx-x)/max(xx))^2+((haz1[*,j]-y)/max(haz1[*,j]))^2
         z=where(dd eq min(dd))
         z=z[0]
         if(nbad eq 1) then begin
            zbad=xx[z]
            plot,xx,haz1(*,j),psym=-1 
            oplot,[zbad,zbad],[haz1[zbad,j],haz1[zbad,j]],psym=1,color=80
         endif else begin 
            zbad=[zbad,xx[z]]
            plot,xx,haz1(*,j),psym=-1 & oplot,zbad,haz1[zbad,j],psym=1,color=80
         endelse
      endif else if(!err eq 2) then begin
         nbad=nbad-1
         zbad=zbad[0:nbad-2]
         if(nbad eq 0) then begin
            plot,xx,haz1(*,j),psym=-1 
         endif else if(nbad eq 1) then begin
            plot,xx,haz1(*,j),psym=-1 
            oplot,[zbad,zbad],[haz1[zbad,j],haz1[zbad,j]],psym=1,color=80
         endif else begin 
            plot,xx,haz1(*,j),psym=-1 & oplot,zbad,haz1[zbad,j],psym=1,color=80
         endelse
      endif
   endwhile
endfor

for j=0,3 do begin
   if(nbad eq 0) then begin
      plot,xx,haz2[*,j],psym=-1
   endif else if(nbad eq 1) then begin
      plot,xx,haz2[*,j],psym=-1 
      oplot,[zbad,zbad],[haz2[zbad,j],haz2[zbad,j]],psym=1,color=80
   endif else begin
      plot,xx,haz2(*,j),psym=-1 & oplot,zbad,haz2[zbad,j],psym=1,color=80
   endelse

   !err=0
   while(!err ne 4) do begin
      print,'click for new point'
      cursor,x,y,3,/down
      if !err eq 1 then begin
         nbad=nbad+1
         dd=((xx-x)/max(xx))^2+((haz2[*,j]-y)/max(haz2[*,j]))^2
         z=where(dd eq min(dd))
         z=z[0]
         if(nbad eq 1) then begin
            zbad=xx[z]
            plot,xx,haz2(*,j),psym=-1 
            oplot,[zbad,zbad],[haz2[zbad,j],haz2[zbad,j]],psym=1,color=80
         endif else begin 
            zbad=[zbad,xx[z]]
            plot,xx,haz2(*,j),psym=-1 & oplot,zbad,haz2[zbad,j],psym=1,color=80
         endelse
      endif else if(!err eq 2) then begin
         nbad=nbad-1
         zbad=zbad[0:nbad-2]
         if(nbad eq 0) then begin
            plot,xx,haz1(*,j),psym=-1 
         endif else if(nbad eq 1) then begin
            plot,xx,haz2(*,j),psym=-1 
            oplot,[zbad,zbad],[haz2[zbad,j],haz2[zbad,j]],psym=1,color=80
         endif else begin 
            plot,xx,haz2(*,j),psym=-1 & oplot,zbad,haz2[zbad,j],psym=1,color=80
         endelse
      endif
   endwhile
endfor

if(nbad eq 0) then begin
   zbad=-999
endif else begin
   zbad=zbad(sort(zbad))
   zbad2=zbad
   for j=0,nbad-1 do begin
      z=where(abs(zbad-zbad[j]) eq (npos-1)/2)
      if(z[0] ne -1) then zbad2[z]=-zbad[z]
   endfor     
   zbad=zbad2
endelse

print,'zbad = ',zbad 
return,zbad
end

