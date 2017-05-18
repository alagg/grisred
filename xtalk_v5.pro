pro xtalk_v5,im

tam=size(im)
nslit=tam[3]

imi=reform(im[0,*,*])>1
imq=reform(im[1,*,*])
imu=reform(im[2,*,*])
imv=reform(im[3,*,*])

step=2.e-4	;1.e-5	;2.e-4 	;5.e-6
xh=findgen(2501)*step
xh=xh-max(xh)/2.

for k=0,nslit-1 do begin
   hq=smooth(histogram(imq[*,k]/imi[*,k],binsize=step,min=min(xh),max=max(xh)),5)
   hu=smooth(histogram(imu[*,k]/imi[*,k],binsize=step,min=min(xh),max=max(xh)),5)
   hv=smooth(histogram(imv[*,k]/imi[*,k],binsize=step,min=min(xh),max=max(xh)),5)

   xi2q=xh(min(where(hq eq max(hq))))+step/2.
   xi2u=xh(min(where(hu eq max(hu))))+step/2.
   xi2v=xh(min(where(hv eq max(hv))))+step/2.

   im[1,*,k]=imq[*,k]-xi2q*imi[*,k]
   im[2,*,k]=imu[*,k]-xi2u*imi[*,k]
   im[3,*,k]=imv[*,k]-xi2v*imi[*,k]

endfor

return
end


