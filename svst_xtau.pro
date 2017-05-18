function svst_xtau,x,tau,th,yy,mm,dd,hh,windows=windows

if(keyword_set(windows) eq 0) then windows=[0,0,0,0]

; Matriz de Mueller de la torre solar sueca t sus derivadas con x(*) y tau(*)
; Modelamos el telescopio como un sistema de tres espejos y las rotaciones
; correspondientes. 
; La primera rotacion la quito porque utilizo esta rutina con un polarizador lineal


; INPUT
; x : cociente entre amplitudes (paralela/perpend) de la onda reflejada 
;     Vector de 3 elementos
; tau : desfase (vector de tres elementos)
; yy : anyo (ej.: 1991)
; mm : mes (ej.: 4)
; dd : dia (ej.: 28)
; hh : hora del dia (ej.: 11.5  --once y media--), TU
; OUTPUT
; mat   : matriz 4x4 del teslecopio SVST 

; th: angulo que forma el banco optico de la sala de observacion
;  th=78 para la mesa del espectrografo
;  th=56 para la mesa del filtro



r_frame_asp,0,0,yy,mm,dd,hh,0,0,d1,d2,d3,d4,d5,d6,d7,d8,raSun,decSun,d9, $
   b0,p,d10,d11,par,haSun,/lapalma


phi=(28.+45./60.)*!pi/180.
eq2altaz,haSun,decSun,phi,altSun,azSun
th1=asin(sin(haSun)*cos(phi)/cos(altsun))
if(haSun lt 0) then th1=-th1+3.*!pi/2. else th1=th1+3.*!pi/2.

azN2=azSun-!pi/2.
; FORMA ALTERNATIVA DE CALCULAR TH1 
;altN1=asin(sin(!pi/4.)*sin(altsun))
;azN1=azN2-acos(cos(!pi/4)/cos(altN1))

;altaz2eq,altN1,azN1,phi,haN1,decN1

;sinth1=-sin(haN1-haSun)*cos(decN1)/sin(!pi/4.)
;th1=asin(sinth1)

th2=!pi/2-altSun

th3=th*!pi/180.-azN2

fac=180./!pi
th1=th1*fac
th2=th2*fac
th3=th3*fac

slit=rotacion(-th1+th2-th3)
mat=espejo(x(0),tau(0))#retarder(windows(0),windows(1))#rotacion(th1)#slit
mat=rotacion(th3)#espejo(x(1),tau(1))#rotacion(th2)#mat
mat=retarder(windows(2),windows(3))#espejo(x(2),tau(2))#mat

mat=mat/mat(0,0,0)
return,mat
end
