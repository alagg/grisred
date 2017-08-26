function grating_angle,lambda

; lambda [A]

deg2rad=!pi/180.
rad2deg=180/!pi

phi=63.4*deg2rad
gamma=2.655*deg2rad

d=1./316e-7     ;groove distance in A

order=fix(2.*d*sin(phi)*cos(gamma)/lambda)

return,order
end

