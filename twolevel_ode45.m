%ode45 ���Two-level���⣬��Two-Level Bloch equations
%Ŀ������C++��rk4 & stepperdopr853 �Ա�

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
rabbi_0=input('����������������Ƶ��rabbi_0: ');
laserchoice=input('ѡ�����壬��������0����˹����1��sin-square����2: ');
IF_CHIRP=input('�Ƿ�����ౣ�����0û�У�����1�У�');
cycles_g=input('���������FWHM������cycles_g��');
delayFlag=input('�Ƿ���˥��������0�����ǣ�����1����: ');
omega_0=input('�����ܼ�ԾǨƵ��omega_0: ');
mu_11_22=input('�ֱ��������ż����mu_11/mu_22(��λ��ż��ԾǨ����Ԫmu�ı���): ');
mu_11=mu_11_22(1);
mu_22=mu_11_22(2);
mu_11=mu*mu_11;
mu_22=mu*mu_22;
xi=(mu_22-mu_11)/(2.0*mu);
y=input('�ֱ����벼��շ���u��v��w�ĳ�ʼֵ: ');
w_init=y(3);
CollectFlag=input('ż���صļ����Ƿ��Ǽ���ЧӦ������0�����ǣ�����1����: ');
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
                rabbi(index1)=rabbi_0*power(sin(pi*(t(index1) + dur)...
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
                rabbi(index1)=rabbi_0*power(sin(pi*(t(index1) + dur)...
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
dipole=Natoms*mu*(yy(:,1)+xi*yy(:,3));  %here DC component is omitted

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
%%laser filed 
if(laserchoice==0)
    %TODO
else if(laserchoice==1)
        envelope_upper=rabbi_0/mu*exp(-4.0*log(2.0)*(t/dur).^2);
        envelope_lower=-envelope_upper;
    else if(laserchoice==2)
            index1=find(t >= -dur & t <= dur);
            envelope_upper(index1)=rabbi_0/mu*power(sin(pi*(t(index1) + dur)...
                    /(2.0*dur)), 2);
            index2=find(t < -dur | t > dur);
            envelope_upper(index2)=0.0;
            envelope_lower=-envelop_upper;
        end
    end
end
figure;
plot(tt/T,laser_field,tt/T,envelope_upper,tt/T,envelope_lower,'linewidth',2)
% title('Laser Field');
set(gca,'linewidth',1);
set(gca,'fontsize',16);
xlabel('Time (\it{T})','fontname','Times New Roman','fontsize',24);
ylabel('\it{E}(\it{t}) (a.u.)','fontname','Times New Roman','fontsize',24);

figure;
plot(tt/T,dipole);
title('dipole');
set(gca,'linewidth',1);
set(gca,'fontsize',16);
xlabel('Time ({\itT})','fontname','Times New Roman','fontsize',24);
ylabel('Intensity (arb.units)','fontname','Times New Roman','fontsize',24);

len=length(dipole);
wmg=2.0*pi/dt;
fre=(0:round(len/2)-1)/len*wmg/omega_L;
FFA=fft(dipole);
hff=FFA;
FFA=abs(FFA(1:round(len/2))/length(dipole));%
FFA=2.0*log10(abs(FFA)+eps);
figure
plot(fre,FFA,'r-','linewidth',1.5);
% title('HHG');
set(gca,'linewidth',1);
set(gca,'fontsize',16);
xlabel('Harmonic Order(\omega_{\itL})','FontName','Times New Roman','fontsize',24);
ylabel('Harmonic Intensity (arb.units)','FontName','Times New Roman','fontsize',24);

%wavelet transform using matlab
freqL=input('Input the Lower limit of freqrange (unit in order): ');
freqU=input('Input the Upper limit of freqrange (unit in order): ');
freqL=freqL*omega_L/(2.0*pi*2.418884326505E-17);  %a.u. to Hz
freqU=freqU*omega_L/(2.0*pi*2.418884326505E-17);  %a.u. to Hz
SamplingRate=power(2,3);
SamplingIndex=1:SamplingRate:length(t);
tSample=t(SamplingIndex);
dSample=dipole(SamplingIndex);
DT=dt*SamplingRate*2.418884326505E-17;
Flag=input('Freqs or Scales is linear spaced, 0 is Freqs, 1 is Scales: ');
wname='cmor1-1';
numScales=512;
cwtdipole=mywavelet(DT,dSample,freqL,freqU,numScales,wname,Flag);

%%%%%
dressed_upper_energy=(-(mu_11+mu_22)/mu*rabbi+sqrt((2.0*xi*rabbi-...
    omega_0).^2+(2.0*rabbi).^2))/2.0;
dressed_lower_energy=(-(mu_11+mu_22)/mu.*rabbi-sqrt((2.0*xi*rabbi-...
    omega_0).^2+(2.0*rabbi).^2))/2.0;
figure;
[hAx,hLine1,hLin2]=plotyy(t/T,[dressed_upper_energy,dressed_lower_energy],...
    t/T,[laser_field,envelope_lower,envelope_upper]);
xlabel('Time (\it{T})','FontName','Times New Roman','FontSize',40);
box off
ax1 = axes('Position',get(hAx(1),'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k');
set(ax1,'YTick', []);
set(ax1,'XTick', []);
box on
ylabel(hAx(1),'\epsilon_{\pm} (a.u.)','FontName','Times New Roman','fontsize',40);
ylabel(hAx(2),'{\it{E}}({\it{t}}) (a.u.)','FontName','Time News Roman','fontsize',40);
set(hAx(1),'linewidth',1);
set(hAx(1),'fontsize',16);
set(hAx(2),'linewidth',1);
set(hAx(2),'fontsize',16);

%%%%%
figure;
subplot(1,2,1);
harmonic_energy=sqrt((2.0*xi*rabbi-omega_0).^2+4.0*rabbi.^2)/omega_L;
plot(harmonic_energy(SamplingIndex),t(SamplingIndex)/T);
set(gca,'Ydir','reverse');
set(gca,'linewidth',1);
set(gca,'fontsize',16);
xlabel('Harmonic Order(\omega_{\itL})','FontName','Times New Roman','fontsize',24);
ylabel('Time ({\itT})','FontName','Times New Roman','fontsize',24);
subplot(1,2,2);
mapsize=256;
colormap(pink(mapsize));
[Tcenter,Freqs]=meshgrid(cwtdipole.frequency.*(2*pi*2.418884326505E-17)...
    /omega_L,tSample/T);
surf(Tcenter,Freqs,wcodemat(cwtdipole.cfs,mapsize)');shading('interp');view(0,90);
set(gca,'Ydir','reverse');
set(gca,'linewidth',1);
set(gca,'fontsize',16);
% % image(t,freqs,wcodemat(abs(Coeffs),mapsize));
% imagesc(cwtdipole.frequency.*(2*pi*2.418884326505E-17)/omega_L,tSample/T,...
%     wcodemat(cwtdipole.cfs',mapsize));

%%synthesize attosecond pulse train
min_order=input('ѡ��г����С���Σ�');   %�ϳ�����ѡ��г�����η�Χ������ʵ�����ȷ��
max_order=input('ѡ��г����󼶴Σ�');
fre=(0:len-1)/len*wmg/omega_L;
order_selected=find(fre<min_order|fre>max_order);
hff2=hff;
hff2(order_selected)=0;
hpulse=abs(ifft(hff2));
hpulse=hpulse.*hpulse;
pulse_intensity=hpulse;
figure
plot(tt(1:length(pulse_intensity))/T,pulse_intensity,'-.b','linewidth',2);%
xlabel('Time(\it{T})','fontsize',14);
ylabel('Pulse Intensity (arb.units)','fontsize',14);
%����ϳ������FWHM
t2=tt(1:length(pulse_intensity));
Max_pulseintensity=max(pulse_intensity);
HalfMax_order=find(abs(pulse_intensity-Max_pulseintensity/2.0)<1e-3&abs(pulse_intensity-Max_pulseintensity/2.0)/Max_pulseintensity<1e-3);
FWHM=abs(2*t2(HalfMax_order(1))*2.41888e-17); %��λ����

    