pro fov_img_v2, image,file,hdr=hdr

if ~keyword_set(hdr) then scale=0.135 else scale = sxpar(hdr,'STEPSIZE')

sd = size(image)

factor = sd[1]/float(sd[2])

!p.charsize=1.6
!p.thick = 5
!x.thick = 5
!y.thick = 5
!p.charthick = 5
!p.ticklen = -.02
!p.font=7
!x.omargin= [1,0]
xtit = 'X (arcsec)'
ytit = 'Y (arcsec)'
xpaper = 20

axesk={xtickformat:'(A1)',ytickformat:'(A1)',xs:13,ys:13,charsize:!p.charsize}
fname = file+'.eps'

;OPEN PS
set_plot,'ps'
device,xsize=xpaper*1.1,ysize=xpaper,xoffset=0., yoffset=0.,/color,$ 
/encapsulate,bits_per_pixel=8,filename=fname,/helvetica

loadct,0,/sil
tit = strmid(file,strpos(file,'level1/')+7)

image_cont_al,/cut,/aspect,/exact,contour=0,image,zrange=minmaxp(image,perc=99.5),title=tit

;cgimage,image,/axis,/keep,tit=tit, AXKEYWORDS=axesk
;axis,xaxis=0,xrange=(!X.CRANGE)*scale,/xs,xtit=xtit,xtickint=round(sd[1]*scale/4)
;axis,yaxis=0,yrange=(!Y.CRANGE)*scale,/ys,ytit=ytit,ytickint=round(sd[2]*scale/4)
;axis,xaxis=1,xrange=(!X.CRANGE)*scale,/xs,xtickformat='(A1)',xtickint=round(sd[1]*scale/4)
;axis,yaxis=1,yrange=(!Y.CRANGE)*scale,/ys,ytickformat='(A1)',ytickint=round(sd[2]*scale/4)

;CLOSE PS
device,/close
set_plot,'x'

end
