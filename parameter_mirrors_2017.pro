pro parameter_mirrors_2017,lambda,x4,tau4,x567,tau567,x8910,tau8910

lamref=[10830.,15650.]

x4_ref=[0.997,0.992]
tau4_ref=[171.133,170.652]
x567_ref=[0.978,1.039]
tau567_ref=[538.888,523.456]
x8910_ref=[1.076,1.054]
tau8910_ref=[497.095,505.807]

x4=interpol(x4_ref,lamref,lambda)
tau4=interpol(tau4_ref,1./lamref,1./lambda)
x567=interpol(x567_ref,lamref,lambda)
tau567=interpol(tau567_ref,1./lamref,1./lambda)
x8910=interpol(x8910_ref,lamref,lambda)
tau8910=interpol(tau8910_ref,1./lamref,1./lambda)

return
end




