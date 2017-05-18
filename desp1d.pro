function desp1d,datin,deltax

fdatin=fft(reform(datin),-1)

tam=size(reform(datin))
if(tam(1) MOD 2 eq 0) then begin
   nux=findgen(tam(1)/2+1)/tam(1)
   nux=[nux,-reverse(nux(1:tam(1)/2-1))]
   fdatin(tam(1)/2)=0.
endif else begin
   nux=findgen(tam(1)/2+1)/tam(1)
   nux=[nux,-reverse(nux(1:tam(1)/2))]
   fdatin(tam(1)/2:tam(1)/2+1)=0.
endelse

phi=-2*!pi*(deltax*nux)
i=complex(0.,1.)
phi=cos(phi)+i*sin(phi)
fdatout=fdatin*phi

datout=float(fft(fdatout,1))

return,datout
end
