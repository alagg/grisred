function redir_gregor_func,lambda,n_orden,show=show

  deg2rad=!pi/180.
  rad2deg=180/!pi

  phi=63.4*deg2rad
  offset=2.655
  gamma=offset*deg2rad

  d=1./316e-7                   ;groove distance in A
  focal=3000.	;effective camera mirror focal length (mm). Includes demagnification of reimaging optics.

  lambdab=2*d*sin(phi)*cos(gamma)/n_orden
  disp=focal*n_orden/d/cos(phi-gamma)

  if keyword_set(show) then begin
    print,'   orden   lambda blaze  dispersion(mm/A)'
    print,n_orden,lambdab,disp
  endif
  
  sinbeta=n_orden*float(lambda)/d/2./cos(gamma)
  if(abs(sinbeta) gt 1) then begin
    print,lambda, ' A not visible in this order'
  endif else begin
    beta=(asin(sinbeta)+gamma)
    disp=focal*n_orden/d/cos(beta)
    lambda1=2*d*sin(beta+gamma)*cos(gamma)/(n_orden-1)
    lambda2=2*d*sin(beta+gamma)*cos(gamma)/(n_orden+1)
    if keyword_set(show) then begin
      print,'   order   lambda        angle   dispersion(mm/A)'
      print,n_orden,lambda,beta*rad2deg-offset,disp
      print,' '
      print,'wavelenghts of previous and next orders'
      print,lambda1,lambda2
    endif
  endelse

  return,disp
end
