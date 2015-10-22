function dy=my_ode_fun(t,y)
global nvar;
global omega_L;
global rabbi_0;
global laserchoice;
global dur;
global IF_CHIRP;
global eta;
global tao;
global delayFlag;
global xi;
global omega_0;
global T1;
global T2;
global w_init;
global CHIRP_PHASE;
CHIRP_PHASE=-IF_CHIRP*eta*tanh(t/tao);
if(laserchoice==0)
    if(t>=-dur/2.0 && t<=dur/2.0)
        rabbi_tmp=rabbi_0*sin(omega_L*t);
    else rabbi_tmp=0.0;
    end
else
if(laserchoice==1)
    rabbi_tmp=rabbi_0*exp(-4.0 * log(2.0)*t*t/dur/dur)*cos(omega_L*t + CHIRP_PHASE);
else 
if(laserchoice==2)
    if(t>=-dur && t<=dur)
        rabbi_tmp=rabbi_0*power(sin(pi*(t + dur) / (2.0*dur))...
            , 2)*cos(omega_L*t + CHIRP_PHASE);
    end
end
end
end

dy=zeros(nvar,1);
dy(1)= -omega_0*y(2) + 2.0*xi*rabbi_tmp * y(2)-delayFlag*y(1)/T2;
dy(2)=  omega_0*y(1) - 2.0*xi*rabbi_tmp * y(1) + 2.0*rabbi_tmp * y(3)-delayFlag*y(2)/T2;
dy(3)= -2.0*rabbi_tmp*y(2)-delayFlag*(y(3)-w_init)/T1;
end


