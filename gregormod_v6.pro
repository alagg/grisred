; ************************   GREGORMOD   ***************************************
;
; IDL-function to model the polarimetric properties of the GREGOR-telescope.
;
; H. Balthasar, October 2009, following ideas of VTT.PRO of MCV.
;

;   Input 
;                    lambda: wavelength (nm), will be changed to Angstrom.
;         (from telescope:)
;                    telescope elevation telescope azimuth and derotator 
;                    orientation.
;                    The derotator is ignored if angle is larger than 360. 
;         (from epehmeris calculation:)
;                    solar declination, zenith distance, an solar axis angle P_0 
;
;   Output:  telescope polarization matrix
;
;   Keywords:     FOLLOW_AO          : to be set if AO-train must be modelled too.
;                 FOLLOW_GRIS        ; to be set if GRIS is used.
;                 FOLLOW_GFPI        : to be set if GFPI is used
;                 FOLLOW_BLUECHANNEL ; to be set if GFPI blue channel is used.   
;                 POLCALIB           : to be set if Polarization Calibration Unit
;                                      is inserted, ignore effects of M1 and M2
;
; ******************************************************************************
; ******************************************************************************
;
;  Auxilliary function N_K_MIRROR from MCV.
;
      function n_k_mirror,lam_ref

; determines the refractive index  n + ik of a mirror
; at a certain wavelength lambda [A] (between 1000 and 50000 A)
; and the reflection coefficient R for normal incidence.
; Output is a 3-component vector [n,k,R]
;
; data reference:
;   CRC Handbook of Chemistry & Physics, D.R. Lide,
;      74th edition, pg 12-109 (1993)
;
; (Aluminum coating)


     E_eV=[0.250,0.300,0.350,0.400,0.500,0.600,0.700,0.800,0.900,1.000,$
           1.100,1.200,1.300,1.400,1.500,1.600,1.700,1.800,1.900,2.000,$
           2.200,2.400,2.600,2.800,3.000,3.200,3.400,3.600,3.800,4.000,$
           4.200,4.400,4.600,4.800,5.000,6.000,6.500,7.000,7.500,8.000,$
           8.500,9.000,9.500,10.000,10.500,11.000,11.500,12.000,12.500]  
; en eV

     n=[8.586,6.759,5.438,4.454,3.072,2.273,1.770,1.444,1.264,1.212,1.201,$
        1.260,1.468,2.237,2.745,2.625,2.143,1.741,1.488,1.304,1.018,0.826,$
        0.695,0.598,0.523,0.460,0.407,0.363,0.326,0.294,0.267,0.244,0.223,$
        0.205,0.190,0.130,0.110,0.095,0.082,0.072,0.063,0.056,0.049,0.044,$
        0.040,0.036,0.033,0.033,0.034]

     k=[48.235,40.960,35.599,31.485,25.581,21.403,18.328,15.955,14.021,12.464,$
        11.181,10.010, 8.949, 8.212, 8.309, 8.597, 8.573, 8.205, 7.821, 7.479,$
        6.846, 6.283, 5.800, 5.385, 5.024, 4.708, 4.426, 4.174, 3.946, 3.740,$
        3.552, 3.380, 3.222, 3.076, 2.942, 2.391, 2.173, 1.983, 1.814, 1.663,$
        1.527, 1.402, 1.286, 1.178, 1.076, 0.979, 0.883, 0.791, 0.700]

     R=[0.9858,0.9844,0.9834,0.9826,0.9817,0.9806,0.9794,0.9778,0.9749,0.9697,$
        0.9630,0.9521,0.9318,0.8852,0.8678,0.8794,0.8972,0.9069,0.9116,0.9148,$
        0.9200,0.9228,0.9238,0.9242,0.9241,0.9243,0.9245,0.9246,0.9247,0.9248,$
        0.9248,0.9249,0.9249,0.9249,0.9244,0.9257,0.9260,0.9262,0.9265,0.9269,$
        0.9272,0.9277,0.9282,0.9286,0.9293,0.9298,0.9283,0.9224,0.9118]

     h=6.62626176e-34
     c=2.99792458e8
     q=1.602189e-19

     factor=h*c*1.e10/q

     lambda=factor/E_eV

     n_ref=interpol(n,lambda,lam_ref)
     k_ref=interpol(k,lambda,lam_ref)
     R_ref=interpol(R,lambda,lam_ref)

     RETURN,[n_ref,k_ref,R_ref]
     END
;
;*******************************************************************************
;*******************************************************************************
;
;  Auxilliary function ROTACION to rotate polarization plane (following MCV)
;
     function ROTACION,theta

     ang=theta*!dtor

     rotmatrix=dblarr(4,4)

     c2=cos(2*ang)
     s2=sin(2*ang)

     rotmatrix[0,0]=1.d0
     rotmatrix[3,3]=1.d0

     rotmatrix[1,1]=c2
     rotmatrix[1,2]=s2

     rotmatrix[2,1]=-s2
     rotmatrix[2,2]=c2

     RETURN,rotmatrix
     END
;
;*******************************************************************************
;*******************************************************************************
;
;    Auxilliary function MIRROR_POL to determine the M"ullermatrix of a mirror
;    following Capitani et al, 1989, Sol. Phys., 120, 173
;
     function MIRROR_POL,n,k,ang,rho,sigma

; ENTRIES:
; N (refractive index = n + i*k
; ang = incidence angle  (in degrees)
;
; OUTPUT:
; mueller =  Mueller matrix
; rho = factor determining linear polarization
; sigma = phase shift (in degrees)


    rad=double(ang)*!dtor
    p=double(n)^2.-double(k)^2.-sin(rad)^2.
    q=4.*n^2.*k^2.

    f2=1.d0/2.*(p+sqrt(p^2.+q))
    g2=1.d0/2.*(-p+sqrt(p^2.+q))

    r=2.d0*sqrt(f2)*sin(rad)*tan(rad)
    s=sin(rad)^2.*tan(rad)^2.

    x2=(f2+g2-r+s)/(f2+g2+r+s)
    rho=sqrt(x2)

    tantau=2.*sqrt(g2)*sin(rad)*tan(rad)/(s-f2-g2)

    sigma=atan(tantau)/!dtor
    if (sigma le 0.) then sigma=sigma+180.

    mueller=dblarr(4,4)

    mueller[0,0]=rho^2.+1.d0
    mueller[1,0]=rho^2.-1.d0
    mueller[0,1]=mueller[1,0]
    mueller[1,1]=mueller[0,0]
    mueller[2,2]=2.d0*rho*cos(sigma*!dtor)
    mueller[3,2]=-2.d0*rho*sin(sigma*!dtor)
    mueller[2,3]=-mueller[3,2]
    mueller[3,3]=mueller[2,2]
    mueller=mueller/2.

    RETURN,mueller
    END
;
; ******************************************************************************
;*******************************************************************************
;
;    Auxilliary function WINDOW_POL to determine the M"ullermatrix of a glass 
;    window with unidirectional stress distribution
;    following Beck et al, 2005, A&A 443, 1053
;    Case of small retardance .
;    calls ROTACION
;
     function WINDOW_POL,angle,retard

; ENTRIES:

; angle = direction of stress  (in degrees)
; retard = retardance (in degrees)
;
; OUTPUT:
; mueller =  Mueller matrix

;    rang = angle * !dtor
    rret = double(retard) * !dtor

    mueller=dblarr(4,4)

    mueller[0,0]=1.d0
    mueller[1,0]=0.
    mueller[0,1]=0.
    mueller[1,1]=mueller[0,0]
    mueller[2,2]=cos(rret)
    mueller[3,2]=sin(rret)
    mueller[2,3]=-mueller[3,2]
    mueller[3,3]=mueller[2,2]
;
    mueller =  rotacion(-angle) # mueller # rotacion(angle)

    RETURN,mueller
    END
;
; ******************************************************************************
; ******************************************************************************
; ******************************************************************************
;
;   Now  function   GREGORMOD
;
     function GREGORMOD_V6, tel_elevation, tel_azimuth, derotatorangle, $
              decSun, p0,  FOLLOW_AO=FOLLOW_AO, FOLLOW_GRIS=FOLLOW_GRIS, $
              FOLLOW_GFPI=FOLLOW_GFPI, FOLLOW_BLUECHANNEL=FOLLOW_BLUECHANNEL, $
              POLCALIB = POLCALIB, LAMBDA=LAMBDA, XTAU=XTAU
;
 
    if(keyword_set(lambda) eq 0 and keyword_set(xtau) eq 0) then begin
       print,'either lambda or xtau is needed in GREGORMOD'
       stop
    endif

    if(keyword_set(lambda) eq 0 and keyword_set(xtau) eq 1 and $
          xtau[0] eq -999) then begin
       print,'6 values (x4,tau4,x567,tau567,x8910,tau8910) are required for xtau'
       print,stop
    endif 

    aocheck =0 
     if (keyword_set(FOLLOW_AO) or keyword_set(FOLLOW_GRIS) or keyword_set(FOLLOW_GFPI) $
         or keyword_set(FOLLOW_BLUECHANNEL) ) then aocheck = 1

     lamb10=10.*lambda
     zenitdist = 90. - tel_elevation
     onemat=fltarr(4,4)
     for i=0,3 do onemat[i,i]=1.
;
     refr_ind=n_k_mirror(lamb10)
     n=refr_ind[0]
     k=refr_ind[1]
; n : real part of refractive index 
; k : imaginary part of refractive index 
;
; positive sense of rotation: counterclockwise observing the Sun 
;
;DEFINITIONS
     latitude=28.+18./60.        ; Izana latitude (degrees)
;
     phir=(90.-latitude)*!dtor
     cosphi=cos(phir)
     sinphi=sin(phir)
     deltar=(90.-decSun)*!dtor
     cosd=cos(deltar)
     sind=sin(deltar)
     zdr=zenitdist*!dtor     
     sinzd=sin(zdr)
     coszd=cos(zdr)
;
; start calculations
;
;  TH1 : angle from solar N-S to telescope elevation plane
;        (rotation angle at M1)
;
     costh1=(cosphi -coszd*cosd)/(sinzd*sind)
     costh1=-1>costh1<1
     th1r=acos(costh1)
     if (tel_azimuth gt 0.) then th1r = -th1r
     th1 =-th1r/!dtor - p0     

;     print,'th1=', th1,th1r/!dtor,p0

;     stop
;
     mm12=onemat
     if NOT KEYWORD_SET(POLCALIB) then begin
 
        inci1=0      ; incidence angle at main mirror M1
        r1=rotacion(th1)
        m1=mirror_pol(n,k,inci1)
;
        inci2=0                     ; incidence angle at M2
        m2=mirror_pol(n,k,inci2)
;
        mm12 = m2 # m1 # r1
;        print,mm12
     endif

     inci3=0.                       ; incidence angle at M3 
     m3=mirror_pol(n,k,inci3)
;
     inci4=45.                      ; incidence angle at M4 
     th4=90.
     r4=rotacion(th4)
;
;www   entrance window of coud'e-train
;www   
;      alpha1 = 15.  ;tbm/tbc! 
;      retard1 = 30.     ; 3. * 6000./lam10    tbm/tbc.
;      w1 = window_pol(alpha1, retard1)
;www
; 
     inci5=30.                      ; incidence angle at M5
     th5= tel_elevation+90      ; rotation at M5
     r5=rotacion(th5)
;
     inci6=7.                       ; incidence at M6
     th6=0                          ; rotation at M6
;     r6=rotacion(th6)
;
     inci7=8.                       ; incidence at M7
     th7=0                          ; rotation at M7
   
     if(keyword_set(xtau) eq 0 or xtau[0] eq -999) then begin
        m4=mirror_pol(n,k,inci4)
        m5=mirror_pol(n,k,inci5)
        m6=mirror_pol(n,k,inci6)
        m7=mirror_pol(n,k,inci7)
        m567=m7#m6#m5
     endif else begin
        m4=espejo(xtau[0],xtau[1])
        m567=espejo(xtau[2],xtau[3])
     endelse

;     r7=rotacion(th7)
;
;www   exit window of coud'e-train
;www   
      alpha2 = -15.      ;tbm/tbc! 
      retard2 = 30      ;3.*6000./lam10    tbm/tbc.
      w2 = window_pol(alpha2, retard2)
;www
; 
     th7_8= -90-tel_azimuth            ; rotation AFTER M7
     r78=rotacion(th7_8)
;
     gregormatrix =  r78 # m567 # r5 # m4 # r4 # m3 # mm12
;activate next line to include vacuum windows!
;     gregormatrix =  r78 # w2 # m567 # w1 # r5  # m4 # r4 # m3 # mm12

;  Now entrance to image derotator if derotatorangle le 360.
;  Assuming North-South direction inside building as zero angle of the derotator.
;   (incidence-excidence plane of the beam)

     if derotatorangle ne -999 then begin

        th8 = derotatorangle             ; rotation of derotator
        r8 = rotacion(th8)
        th10 = -derotatorangle           ; rotation after derotator   
        r10=rotacion(th10)

        if(keyword_set(xtau) eq 0 or xtau[0] eq -999) then begin
           inci8 = 60.                      ; incidence to M8

           m8 = mirror_pol(n,k,inci8)
;
           inci9 = 30.                      ; incidence at M9
           th9 = 0.
           m9 = mirror_pol(n,k,inci9)
;        r9=rotacion(th9)
;
           inci10 = 60.                     ; incidence at M10 
           m10 = mirror_pol(n,k,inci10)
           m8910=m10#m9#m8
        endif else begin
           m8910=espejo(xtau[4],xtau[5])
        endelse
;   exit derotator
;
        gregormatrix = r10 # m8910 # r8 # gregormatrix
;
     endif
; XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
;
;    Interruption here because of possible second calibration unit.
;
; now following AO-train if keyword is set..
;
    if aocheck eq 1  then begin

;      Now deflection mirror M11
        inci11 = 45.                 ; incidence at M11
        th11 = 150.                   ; rotation at M11
        m11=mirror_pol(n,k,inci11)
        r11=rotacion(th11)
;
;      AO-collimator    
        inci12 = 1.5                 ; incidence at M12
        th12 = 0.                  ; rotation at M12
        m12 = mirror_pol(n,k,inci12)
        r12 = rotacion(th12)
;
;      Deformable mirror
        inci13 = 15.                 ; incidence at M13 
        th13 = 0
        m13 = mirror_pol(n,k,inci13)
        r13 = rotacion(th13)
;
;      Tip-tilt
        inci14 = 15.                ; incidence  at M14
        th14 = 0
        m14 = mirror_pol(n,k,inci14)
        r14 = rotacion(th14)
;
;      AO-camera mirror
        inci15=5.5
        th15=0
        m15=mirror_pol(n,k,inci15)
        r15=rotacion(th15)
;
        gregormatrix = m15#r15#m14#r14#m13#r13#m12#r12#m11#r11#gregormatrix
;
     endif

     if keyword_set(FOLLOW_GRIS) then begin

;      PF-feeding mirror M16
        sidea = 5.5 * !dtor
        sideb = (27. + (180.-th11)) * !dtor
        sidec = acos(cos(sidea)*cos(sideb))
        inci16 = (0.5/!dtor) * sidec
        th16 = (1./!dtor) * asin(sin(sideb)/sin(sidec))
        th17 = (1./!dtor) * asin(sin(sidea)/sin(sidec))
        m16 = mirror_pol(n,k,inci16)
        r16 = rotacion(th16)
        r17 = rotacion(th17)
;           r17: plane parallel horizontal slit.
;
        gregormatrix = r17#m16#r16#gregormatrix
     endif

     if keyword_set(FOLLOW_GFPI) then begin
;        effects of lenses, straight passed beamsplitters 
;        and FPIs ignored so far 

;      PF-feeding mirror M16
        sidea = 5.5 * !dtor
        sideb = (27. + (180.-th11)) * !dtor
        sidec = acos(cos(sidea)*cos(sideb))
        inci16 = (0.5/!dtor) * sidec
        th16 = (1./!dtor) * asin(sin(sideb)/sin(sidec))
        th17 = (1./!dtor) * asin(sin(sidea)/sin(sidec))
        m16 = mirror_pol(n,k,inci16)
        r16 = rotacion(th16)
        r17 = rotacion(th17)
;        rotation to have the reflecting plane of following elemnnts
;        horizontal.

;      first pentaprism treated as two Al-mirrors (zero-approach!!!!!) 
        inci17 = 22.5
        m17 = mirror_pol(n,k,inci17)
;        m17  twice!

;      45 degree mirror
        th18 = 0.
        inci18 = 45.        
;        r18 = rotacion(th18)
        m18 = mirror_pol(n,k,inci18)

; mirrors 19 - 22 reserved for FPIs

;       last folding mirror
;        th23 = 90.
        inci23 = 16.845         ; check this angle!
        m23 = mirror_pol(n,k,inci23)
;        r23 = rotacion(th23)

        gregormatrix = m23#m18#m17#m17#r17#m16#r16#gregormatrix

     endif

     if keyword_set(FOLLOW_BLUECHANNEL) then begin
;        effects of lenses, straight passed beamsplitters 
;        ignored so far 

;      PF-feeding mirror M16
        sidea = 5.5 * !dtor
        sideb = (27. + (180.-th11)) * !dtor
        sidec = acos(cos(sidea)*cos(sideb))
        inci16 = (0.5/!dtor) * sidec
        th16 = (1./!dtor) * asin(sin(sideb)/sin(sidec))
        th17 = (1./!dtor) * asin(sin(sidea)/sin(sidec))
        m16 = mirror_pol(n,k,inci16)
        r16 = rotacion(th16)
        r17 = rotacion(th17)
;        rotation to have the reflecting plane of following elemnnts
;        horizontal.

;      first pentaprism treated as two Al-mirrors (zero-approach!) 
        inci17 = 22.5
        m17 = mirror_pol(n,k,inci17)
;        m17 twice!
;
;      second pentaprism treated as two Al-mirrors (zero-approach!) 
        inci18 = 22.5
        m18 = mirror_pol(n,k,inci18)
;        m18 twice!

        gregormatrix = m18#m18#m17#m17#r17#m16#r16#gregormatrix

     endif
;
     RETURN,gregormatrix
     END

