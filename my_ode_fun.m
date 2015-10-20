function dy=my_ode_fun(t,y)
global IF_CHIRP;
global eta;
global tao;
global dur;
global xi;
global omega_L;
global omega_0;
global rabbi_0;
global CHIRP_PHASE;
CHIRP_PHASE=-IF_CHIRP*eta*tanh(t/tao);
rabbi_tmp=rabbi_0*exp(-4.0 * log(2.0)*t*t/dur/dur)*cos(omega_L*t + CHIRP_PHASE);
dy=zeros(3,1);
dy(1)= -omega_0*y(2) + 2.0*xi*rabbi_tmp * y(2);
dy(2)=  omega_0*y(1) - 2.0*xi*rabbi_tmp * y(1) + 2.0*rabbi_tmp * y(3);
dy(3)= -2.0*rabbi_tmp*y(2);
end


