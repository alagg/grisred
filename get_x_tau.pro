function get_x_tau,n,k,ang

rad=ang*!pi/180.

p=n^2.-k^2.-sin(rad)^2.
q=4.*n^2.*k^2.

f2=(p+sqrt(p^2.+q))/2.
g2=(-p+sqrt(p^2.+q))/2.

if(k eq 0) then r=0 else r=2.*sqrt(f2)*sin(rad)*tan(rad) 
s=sin(rad)^2.*tan(rad)^2.

x2=(f2+g2-r+s)/(f2+g2+r+s)
x=sqrt(x2)

tau=atan(2.*sqrt(g2)*sin(rad)*tan(rad),(s-f2-g2))*180./!pi

return,[x,tau]
end
