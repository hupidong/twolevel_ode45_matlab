function Coeffs_Struct=mywavelet(dt,sample,freqL,freqU,numscales,Flag)
my_wavelet_fun='cmor1-1';
fc=centfrq(my_wavelet_fun);   %Hz
freqrange=[freqL, freqU];
DT=dt*2.418884326505E-17;
scalerange=fc./(freqrange.*DT)';
numScales=numscales; 
scaleLower=scalerange(1);
scaleUpper=scalerange(2);
if(Flag==0) %freqs linear spaced
    power_of_scaleLower=log10(scaleLower);
    power_of_scaleUpper=log10(scaleUpper);
    scales=logspace(power_of_scaleLower,power_of_scaleUpper,numScales);
    freqs=scal2frq(scales,my_wavelet_fun,DT);
else
if(Flag==1) %scales linear spaced
   scales=linspace(scaleLower,scaleUpper,numScales);
   freqs=scal2frq(scales,my_wavelet_fun,DT);
else
    disp('Error for Choosing scales or freqs linear spaced!!!');
end
end
Coeffs=cwt(sample,scales,my_wavelet_fun);
Coeffs_Struct=struct('cfs',abs(Coeffs),'scales',scales,'frequency',freqs);
end
