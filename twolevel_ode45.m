%ode45 解决Two-level问题，即Two-Level Bloch equations
%目的是与C++的rk4 & stepperdopr853 对比

clear;
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

%some const varibles
nvar=3;
omega_L=0.056;
T=2.0*pi/omega_L;
eta=6.25;
tao=120.0;
T1=1.0E-12 / 2.418884326505E-17;
T2=0.5E-12 / 2.418884326505E-17;
mu=1.0e-29 / 1.60217653e-19 / 5.291772108e-11;
Natoms=7.5E24*power(5.29E-11, 3);
atol=1.0e-10;rtol=atol;

%varibles needed to assign values
rabbi_0=input('输入驱动脉冲拉比频率rabbi_0: ');
laserchoice=input('选择脉冲，矩形输入0，高斯输入1，sin-square输入2: ');
IF_CHIRP=input('是否有啁啾，数字0没有，数字1有：');
cycles_g=input('输入脉冲的FWHM周期数cycles_g：');
delayFlag=input('是否考虑衰减，数字0不考虑，数字1考虑: ');
omega_0=input('输入能级跃迁频率omega_0: ');
mu_11_22=input('分别输入固有偶极矩mu_11/mu_22(单位以偶极跃迁矩阵元mu的倍数): ');
mu_11=mu_11_22(1);
mu_22=mu_11_22(2);
mu_11=mu*mu_11;
mu_22=mu*mu_22;
xi=(mu_22-mu_11)/(2.0*mu);
y=input('分别输入布洛赫方程u、v和w的初始值: ');
w_init=y(3);
CollectFlag=input('偶极矩的计算是否考虑集体效应，数字0不考虑，数字1考虑: ');
if(CollectFlag==0)
    Natoms=1.0;
else if(CollectFlag==1)
        
    else disp('Error Input for "CollectFlag"!!!');
    end
end

%laser field
dur=cycles_g*T;
if laserchoice==0
    cycles=cycles_g+2.0;
else if laserchoice==1
        cycles=cycles_g*4.0;
    else if laserchoice==2
            cycles=cycles_g*2.0+2.0;
        else disp('Error Input for "laserchoice"!!!');
        end
    end
end
n=13;
N=power(2,n);
dt=T/N;
NN=N*cycles;
NT=NN+1;
tstart=-cycles/2.0*T;
tend=-tstart;

t=linspace(tstart,tend,NT)';
if(IF_CHIRP==0)
    if(laserchoice==0)
        index1=find(t>=-dur/2.0 & t<=dur/2.0);
        rabbi(index1)= rabbi_0*sin(omega_L*t(index1));
        index2=find(t<-dur/2.0 | t>dur/2.0);
        rabbi(index2)=0.0;
    else if(laserchoice==1)
            rabbi= rabbi_0*exp(-4.0 * log(2.0)*t.*t / dur / dur)...
                .*cos(omega_L*t);
        else if(laserchoice==2)
                index1=find(t >= -dur & t <= dur);
                rabbi(index1)=rabbi_0*power(sin(pi*(t(index) + dur)...
                    /(dur*2.0)), 2)*cos(omega_0*t(index1));
                index2=find(t < -dur | t > dur);
                rabbi(index2)=0.0;
            end
        end
    end
else if(IF_CHIRP==1)
        if(laserchoice==0)
        index1=find(t>=-dur/2.0 & t<=dur/2.0);
        rabbi(index1)= rabbi_0*sin(omega_L*t(index1)-eta*tanh(t(index1)/tao));
        index2=find(t<-dur/2.0 | t>dur/2.0);
        rabbi(index2)=0.0;
    else if(laserchoice==1)
            rabbi= rabbi_0*exp(-4.0 * log(2.0)*t.*t / dur / dur)...
                .*cos(omega_L*t-eta*tanh(t/tao));
        else if(laserchoice==2)
                index1=find(t >= -dur & t <= dur);
                rabbi(index1)=rabbi_0*power(sin(pi*(t(index) + dur)...
                    /(2.0*dur)), 2)*cos(omega_0*t(index1)-eta*tanh(t(index1)/tao));
                index2=find(t < -dur | t > dur);
                rabbi(index2)=0.0;
            end
        end
        end
    end
end
laser_field=rabbi/mu;

options=odeset('RelTol',rtol,'AbsTol',atol);
tspan=t;
[tt,yy]=ode45(@my_ode_fun,tspan,y,options);
dipole=Natoms*mu*(yy(:,1)+xi*yy(:,3));

%write some results to txt files
ftime=fopen('res\time.txt','wt');
fprintf(ftime,'%12.10e\n',tt/T);
fclose(ftime);
ffield=fopen('res\field.txt','wt');
fprintf(ffield,'%12.10e\n',laser_field);
fclose(ffield);
fdipole=fopen('res\dipole.txt','wt');
fprintf(fdipole,'%12.10e\n',dipole);
fclose(fdipole);

%data process and results display
figure;
plot(tt/T,laser_field,'linewidth',2)
title('Laser Field');
xlabel('Time (T)','fontsize',14);
ylabel('E(t) (a.u.)','fontsize',14);

figure;
plot(tt/T,dipole);
title('dipole');
xlabel('Time (T)','fontsize',14);
ylabel('Intensity (arb.unit)','fontsize',14);

len=length(dipole);
wmg=2.0*pi/dt;
fre=(0:round(len/2)-1)/len*wmg/omega_L;
FFA=fft(dipole);
hff=FFA;
FFA=abs(FFA(1:round(len/2))/length(dipole));%
FFA=2.0*log10(abs(FFA)+eps);
figure
plot(fre,FFA,'r-','linewidth',2);
title('HHG');
xlabel('Harmonic Order(\omega/\omega_L)','fontsize',14);
ylabel('Harmonic Intensity(arb.units)','fontsize',14);

%wavelet transform using matlab
freqL=input('Input the Lower limit of freqrange (unit in order): ');
freqU=input('Input the Upper limit of freqrange (unit in order): ');
freqL=freqL*omega_L/(2.0*pi*2.418884326505E-17);  %a.u. to Hz
freqU=freqU*omega_L/(2.0*pi*2.418884326505E-17);  %a.u. to Hz
SamplingRate=power(2,3);
SamplingIndex=1:SamplingRate:length(t);
tSample=t(SamplingIndex);
dSample=dipole(SamplingIndex);
DT=dt*SamplingRate;
Flag=input('Freqs or Scales is linear spaced, 0 is Freqs, 1 is Scales: ');
cwtdipole=mywavelet(DT,dSample,freqL,freqU,512,Flag);

figure;
subplot(2,1,1);
harmonic_energy=sqrt((2.0*xi*rabbi-omega_0).^2+4.0*rabbi.^2)/omega_L;
plot(harmonic_energy(SamplingIndex),t(SamplingIndex));
subplot(2,1,2);
mapsize=256;
colormap(pink(mapsize));
[Tcenter,Freqs]=meshgrid(cwtdipole.frequency,tSample);
% % surf(Tcenter,Freqs,wcodemat(cwtdipole.cfs,mapsize)');shading('interp');view(0,90);
% % image(t,freqs,wcodemat(abs(Coeffs),mapsize));
imagesc(cwtdipole.frequency.*(2*pi*2.41888E-17)/omega_L,tSample,...
    wcodemat(cwtdipole.cfs',mapsize));
% colorbar;

    