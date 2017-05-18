pro get_teles, filecal, pzero=pzero, lambda=lambda, xtau=xtau, $
    display=display,verbose=verbose

if(keyword_set(display) eq 1) then display=1 else display=0
if(keyword_set(verbose) eq 1) then verbose=1 else verbose=0
if(keyword_set(pzero) eq 0) then pzero=fltarr(n_elements(filecal))-999.

nfiles=n_elements(filecal)

if(nfiles lt 2) then begin
   print,'at least 2 telescope calibration files are needed'
   return
end

npzero=n_elements(pzero)
if(nfiles ne npzero) then begin
   print,'the same number of calibration files and offset angles is required'
   return
end

; leemos los ficheros de calibracion, separando los dos haces
; y calculamos los parametros pertinentes

time=fltarr(nfiles)
date=intarr(3,nfiles)

th4=90.
th5=fltarr(nfiles)
th7_8=fltarr(nfiles)

r4=rotacion(th4)
r5=fltarr(4,4,nfiles)
r78=fltarr(4,4,nfiles)

inci3=0.
inci4=45.                      ; incidence angle at M4 
inci5=30.                      ; incidence angle at M5
inci6=7.                       ; incidence at M6
inci7=8.                       ; incidence at M7

format=['(i2,$)','(i3,$)','(i4,$)','(i5,$)','(i6,$)']

for j=0,nfiles-1 do begin

   read_filecal,filecal[j],cal_up,cal_down,time=timecal,date=datecal,$
      pzero=pzero[j],lambda=lambda,luzin=input

   if(j eq 0) then begin
      tam=size(cal_up)
      beam_up=fltarr(tam[1],tam[2],nfiles)
      beam_down=fltarr(tam[1],tam[2],nfiles)
      tam=size(input)
      luzin=fltarr(tam[1],tam[2],nfiles)
   endif

   beam_up[*,*,j]=cal_up
   beam_down[*,*,j]=cal_down
   luzin[*,*,j]=input

   time[j]=mean(timecal)
   date[*,j]=datecal
      
;   calculamos los ángulos del telescopio 
 
   elaz=get_elaz(date[0,j],date[1,j],date[2,j],time[j])

   th5[j]= 90. + elaz[0]       ; rotation at M5
   r5[*,*,j]=rotacion(th5[j])

   th7_8[j]= -elaz[1]            ; rotation AFTER M7
   r78[*,*,j]=rotacion(th7_8[j])
endfor

; calculamos la matriz de demodulacion de referencia

mod_up  =beam_up[*,*,0]  #invert_svd(luzin[*,*,0])
mod_down=beam_down[*,*,0]#invert_svd(luzin[*,*,0])

dmod_up  =invert(mod_up)
dmod_down=invert(mod_down)

dmod_up  =dmod_up  /total(dmod_up[0,*])
dmod_down=dmod_down/total(dmod_down[0,*])

; calculamos la matrices del telescopio y sus derivadas

tam=size(beam_up)
teor_up          = fltarr(tam[1],tam[2],nfiles)
teor_down        = fltarr(tam[1],tam[2],nfiles)
dteor_up_x4      = fltarr(tam[1],tam[2],nfiles)
dteor_down_x4    = fltarr(tam[1],tam[2],nfiles)
dteor_up_tau4    = fltarr(tam[1],tam[2],nfiles)
dteor_down_tau4  = fltarr(tam[1],tam[2],nfiles)
dteor_up_x567    = fltarr(tam[1],tam[2],nfiles)
dteor_down_x567  = fltarr(tam[1],tam[2],nfiles)
dteor_up_tau567  = fltarr(tam[1],tam[2],nfiles)
dteor_down_tau567= fltarr(tam[1],tam[2],nfiles)
dteor_up_cte     = fltarr(tam[1],tam[2],nfiles,nfiles)
dteor_down_cte   = fltarr(tam[1],tam[2],nfiles,nfiles)

lambda=15650.
nk=n_k_mirror(lambda)

xtau3=get_x_tau(nk[0],nk[1],inci3)
xtau4=get_x_tau(nk[0],nk[1],inci4)
xtau5=get_x_tau(nk[0],nk[1],inci5)
xtau6=get_x_tau(nk[0],nk[1],inci6)
xtau7=get_x_tau(nk[0],nk[1],inci7)

m3=espejo(xtau3[0],xtau3[1])

x4=xtau4[0]
tau4=xtau4[1]
x567=xtau5[0]*xtau6[0]*xtau7[0]
tau567=xtau5[1]+xtau6[1]+xtau7[1]
cte=fltarr(nfiles)+1.

delta_cte=fltarr(nfiles)+1
delta_x4=1.
delta_tau4=1.
delta_x567=1.
delta_tau567=1.

error=1.

sigma= 0

if(verbose eq 1) then begin
   print,'x4,tau4,x567,tau567 = ',x4,tau4,x567,tau567         ;,sigma,cte[1:*]
endif

while (error ge 1.e-5) do begin

   m4=espejo(x4,tau4)
   dm4_x=despejo_x(x4,tau4)
   dm4_tau=despejo_tau(x4,tau4)

   m567=espejo(x567,tau567)
   dm567_x=despejo_x(x567,tau567)
   dm567_tau=despejo_tau(x567,tau567)

   telref           =  r78[*,*,0]# m567    #r5[*,*,0]# m4    #r4#m3
   invtelref        =  invert(telref)
   dtelref_x4       =  r78[*,*,0]# m567    #r5[*,*,0]#dm4_x  #r4#m3
   dtelref_tau4     =  r78[*,*,0]# m567    #r5[*,*,0]#dm4_tau#r4#m3
   dtelref_x567     =  r78[*,*,0]#dm567_x  #r5[*,*,0]# m4    #r4#m3
   dtelref_tau567   =  r78[*,*,0]#dm567_tau#r5[*,*,0]# m4    #r4#m3
   dinvtelref_x4    = -invtelref#dtelref_x4    #invtelref
   dinvtelref_tau4  = -invtelref#dtelref_tau4  #invtelref
   dinvtelref_x567  = -invtelref#dtelref_x567  #invtelref
   dinvtelref_tau567= -invtelref#dtelref_tau567#invtelref

   for j=1,nfiles-1 do begin
      tel        =r78[*,*,j]# m567    #r5[*,*,j]# m4    #r4#m3
      dtel_x4    =r78[*,*,j]# m567    #r5[*,*,j]#dm4_x  #r4#m3
      dtel_tau4  =r78[*,*,j]# m567    #r5[*,*,j]#dm4_tau#r4#m3
      dtel_x567  =r78[*,*,j]#dm567_x  #r5[*,*,j]# m4    #r4#m3
      dtel_tau567=r78[*,*,j]#dm567_tau#r5[*,*,j]# m4    #r4#m3

      teor_up[*,*,j]  =mod_up  #invtelref#tel#luzin[*,*,j]
      teor_down[*,*,j]=mod_down#invtelref#tel#luzin[*,*,j]


      dteor_up_x4[*,*,j]       = mod_up    #   $
               (dinvtelref_x4#tel + invtelref#dtel_x4) #luzin[*,*,j]
      dteor_down_x4[*,*,j]     = mod_down  #   $
               (dinvtelref_x4#tel + invtelref#dtel_x4) #luzin[*,*,j]

      dteor_up_tau4[*,*,j]     = mod_up    #   $
               (dinvtelref_tau4#tel + invtelref#dtel_tau4) #luzin[*,*,j]
      dteor_down_tau4[*,*,j]   = mod_down  #   $
               (dinvtelref_tau4#tel + invtelref#dtel_tau4) #luzin[*,*,j]

      dteor_up_x567[*,*,j]     = mod_up    #   $
               (dinvtelref_x567#tel + invtelref#dtel_x567) #luzin[*,*,j]
      dteor_down_x567[*,*,j]   = mod_down  #   $
               (dinvtelref_x567#tel + invtelref#dtel_x567) #luzin[*,*,j]

      dteor_up_tau567[*,*,j]   = mod_up    #   $
               (dinvtelref_tau567#tel + invtelref#dtel_tau567) #luzin[*,*,j]
      dteor_down_tau567[*,*,j] = mod_down  #   $
               (dinvtelref_tau567#tel + invtelref#dtel_tau567) #luzin[*,*,j]


      dteor_up_cte[*,*,j,j]   = teor_up[*,*,j]
      dteor_down_cte[*,*,j,j] = teor_down[*,*,j]

      teor_up[*,*,j]      =cte[j]*teor_up[*,*,j]
      teor_down[*,*,j]    =cte[j]*teor_down[*,*,j]

      dteor_up_x4[*,*,j]      =cte[j]*dteor_up_x4[*,*,j]
      dteor_down_x4[*,*,j]    =cte[j]*dteor_down_x4[*,*,j]
      dteor_up_tau4[*,*,j]    =cte[j]*dteor_up_tau4[*,*,j]
      dteor_down_tau4[*,*,j]  =cte[j]*dteor_down_tau4[*,*,j]
      dteor_up_x567[*,*,j]    =cte[j]*dteor_up_x567[*,*,j]
      dteor_down_x567[*,*,j]  =cte[j]*dteor_down_x567[*,*,j]
      dteor_up_tau567[*,*,j]  =cte[j]*dteor_up_tau567[*,*,j]
      dteor_down_tau567[*,*,j]=cte[j]*dteor_down_tau567[*,*,j]

   endfor

   desv_up   = beam_up[*,*,1]  -teor_up[*,*,1]
   for j=2,nfiles-1 do desv_up  = [ [desv_up],[beam_up[*,*,j] - teor_up[*,*,j]] ]
   desv_down = beam_down[*,*,1]-teor_down[*,*,1]
   for j=2,nfiles-1 do desv_down= [ [desv_down],[beam_down[*,*,j] - teor_down[*,*,j]] ]
   desv=[[desv_up],[desv_down]]
   desv=reform(desv,n_elements(desv))

   xx_up   = [ [[dteor_up_x4[*,*,1]]]   , [[dteor_up_tau4[*,*,1]]],$
               [[dteor_up_x567[*,*,1]]] , [[dteor_up_tau567[*,*,1]]], $
               [[dteor_up_cte[*,*,1:*,1]]] ]

   for j=2,nfiles-1 do begin
      dum  = [ [[dteor_up_x4[*,*,j]]]   , [[dteor_up_tau4[*,*,j]]],$
               [[dteor_up_x567[*,*,j]]] , [[dteor_up_tau567[*,*,j]]], $
               [[dteor_up_cte[*,*,1:*,j]]] ]
      xx_up= [ [xx_up],[dum] ]
   endfor
   
   tam=size(xx_up)
   xx_up=reform(xx_up,tam[1]*tam[2],tam[3])

   xx_down = [ [[dteor_down_x4[*,*,1]]]   , [[dteor_down_tau4[*,*,1]]],$
               [[dteor_down_x567[*,*,1]]] , [[dteor_down_tau567[*,*,1]]], $
               [[dteor_down_cte[*,*,1:*,1]]] ]

   for j=2,nfiles-1 do begin
      dum  = [ [[dteor_down_x4[*,*,j]]]   , [[dteor_down_tau4[*,*,j]]],$
               [[dteor_down_x567[*,*,j]]] , [[dteor_down_tau567[*,*,j]]], $
               [[dteor_down_cte[*,*,1:*,j]]] ]
      xx_down= [ [xx_down],[dum] ]
   endfor
   
   xx_down=reform(xx_down,tam[1]*tam[2],tam[3])

   xxx=[xx_up,xx_down]

   corr=lstsvd(xxx,desv,desvfit)
   error=max(abs(corr/[x4,tau4,x567,tau567,cte[1:*]]))

   x4=x4+corr[0]
   tau4=tau4+corr[1]
   x567=x567+corr[2]
   tau567=tau567+corr[3]
   cte[1:*]=cte[1:*]+corr[4:*]

   sigma=std(desv-desvfit)

   if(verbose eq 1) then begin
      print,'x4,tau4,x567,tau567 = ',x4,tau4,x567,tau567         ;,sigma,cte[1:*]
   endif
endwhile

if(display eq 1) then begin
   window,0,xsize=1200,ysize=600
   !p.multi=[0,2,2]
   loadct,2,/silent

   plot,beam_up[0,*,1:*],/xsty & oplot,teor_up[0,*,1:*],color=80
   plot,beam_up[1,*,1:*],/xsty & oplot,teor_up[1,*,1:*],color=80
   plot,beam_up[2,*,1:*],/xsty & oplot,teor_up[2,*,1:*],color=80
   plot,beam_up[3,*,1:*],/xsty & oplot,teor_up[3,*,1:*],color=80

   !p.multi=0
   loadct,0,/silent
endif

xtau=[x4,tau4,x567,tau567]

return
end
