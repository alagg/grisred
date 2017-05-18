pro calib2file_v5,filecal,fileout,lambda

data=0
gris_v5,filecal,filecal,filecal,lambda=lambda,data=data,/get

get_calib_v5,filecal,dmod1,dmod2,data=data,lambda=lambda,zbad=-999

tam=size(dmod1)
dmod1_av=total(dmod1,3)/tam(3)

openw,1,fileout
printf,1,dmod1_av,format='(4f9.4)'
close,1

return
end
