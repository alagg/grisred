function effic,mat1,nocomment=nocomment

if(keyword_set(nocomment) eq 1) then nocomment=1 else nocomment=0
tam=size(mat1)
nmed=tam(1)

d=transpose(double(mat1))#mat1

if(determ(d,/check) ne 0) then begin
   d=invert(d)
   eps=fltarr(tam(2))
   for j=0,tam(2)-1 do eps(j)=1./sqrt(nmed*d(j,j))
   eps=[eps,sqrt(total(eps(1:tam(2)-1)*eps(1:tam(2)-1)))]
endif else begin
   if(nocomment ne 1) then print,'Modulacion ineficiente'
   eps=fltarr(5)
endelse

return,eps
end
