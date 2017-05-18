;
; ******************************************************************************
; ******************************************************************************
; ******************************************************************************
;
;   From  function   GREGORMOD
;

  function get_elaz,yy,mm,dd,hh,decsun,p0

   r_frame_asp,0,0,yy,mm,dd,hh,0,0,d1,d2,d3,d4,d5,d6,d7,d8,raSun,decSun,$
      d9,b0,p0,d10,d11,par,haSun,/tenerife

   phi=28.299083
   delta0=decSun
   h0=haSun

   phir=phi*!pi/180.
   cosphi=cos(phir)
   sinphi=sin(phir)
   deltar=delta0; 	*!pi/180.
   cosd=cos(deltar)
   sind=sin(deltar)
   h0rad=h0;	/15.
   cosh0=cos(h0rad)
   sinh0=sin(h0rad)
   sinalt0=cosd*cosh0*cosphi+sind*sinphi
   alt0=asin(sinalt0)
   cosalt0=cos(alt0)
   sinaz0=sinh0*cosd/cosalt0
   cosaz0=(cosd*cosh0*sinphi-sind*cosphi)/cosalt0
   az0=atan(sinaz0,cosaz0)

   tel_elev=alt0*180./!pi
   tel_az=az0*180./!pi
   decsun=decsun*180./!pi
   p0=p0*180./!pi

   par=par*180./!pi
   RETURN,[tel_elev,tel_az,par]
END

