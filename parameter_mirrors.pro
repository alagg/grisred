pro parameter_mirrors,lambda,x4,tau4,x567,tau567,x8910,tau8910

lamref=[10830.,15650.]
x4_ref=[0.973,0.9922]
tau4_ref=[171.07,173.5]

x4=interpol(x4_ref,lamref,lambda)
tau4=interpol(tau4_ref,1./lamref,1./lambda)
x567_ref=[1.005,0.9964]
tau567_ref=[536.79,529.41]+[0.,-3.]

x567=interpol(x567_ref,lamref,lambda)
tau567=interpol(tau567_ref,1./lamref,1./lambda)

lamref=[10830.,11000., 12000.,13000.,15650.,16000.,17000.]
x8910_ref=  [1.000, 1.000, 1.000, 1.000, 1.000, 0.999, 0.997]
tau8910_ref=[496.6, 497.4, 501.4, 504.6, 510.8, 511.4, 513.2]

x8910=interpol(x8910_ref,lamref,lambda)
;tau8910=interpol(tau8910_ref,1./lamref,1./lambda)
tau8910=542.276 - 4.92761e5/lambda

return
end




