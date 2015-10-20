%ode45 ���Two-level���⣬��Two-Level Bloch equations
%Ŀ������C++��rk4 & stepperdopr853 �Ա�

clear;
global IF_CHIRP;
global eta;
global tao;
global dur;
global xi;
global omega_L;
global omega_0;
global rabbi_0;
global CHIRP_PHASE;


nvar=3;
mu=1.0e-29 / 1.60217653e-19 / 5.291772108e-11;
eta=6.25;
tao=120.0;
omega_L=0.056;
T=2.0*pi/omega_L;
atol=1.0e-10;rtol=atol;

y=zeros(nvar,1);
y=input('�ֱ����벼��շ���u��v��w�ĳ�ʼֵ: ');
omega_0=input('�����ܼ�ԾǨƵ��omega_0: ');
mu_11_22=input('�ֱ��������ż����mu_11/mu_22(��λ��ż��ԾǨ����Ԫmu�ı���): ');
mu_11=mu_11_22(1);
mu_22=mu_11_22(2);
mu_11=mu*mu_11;
mu_22=mu*mu_22;
xi=(mu_22-mu_11)/(2.0*mu);
rabbi_0=input('����������������Ƶ��rabbi_0: ');
cycles_g=input('���������FWHM������cycles_g��');
dur=cycles_g*T;
cycles=cycles_g*3;  %for Gaussian pulse
IF_CHIRP=input('�Ƿ�����ౣ�����0û�У�����1�У�');

n=13;
N=power(2,n);
NN=N*cycles;
NT=NN+1;
tstart=-cycles/2.0*T;
tend=-tstart;

t=linspace(tstart,tend,NT);
if(IF_CHIRP==0)
   rabbi= rabbi_0*exp(-4.0 * log(2.0)*t.*t / dur / dur).*cos(omega_L*t);
else
   rabbi= rabbi_0*exp(-4.0 * log(2.0)*t.*t / dur / dur).*cos(omega_L*t-eta*tanh(t/tao));
end
laser_field=rabbi/mu;

options=odeset('RelTol',rtol,'AbsTol',atol);
tspan=t;
[T,Y]=ode45(@my_ode_fun,tspan,y,options);

figure;
dipole=mu*Y(:,1);
plot(T,dipole);

    