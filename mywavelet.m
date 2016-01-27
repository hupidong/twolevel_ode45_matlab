function Coeffs_Struct=mywavelet(dt,sample,freqL,freqU,numscales,wname,Flag)
my_wavelet_fun=wname;
fc=centfrq(my_wavelet_fun);   %Hz
DT=dt;
numScales=numscales; 
if(Flag==0) %freqs linear spaced
%     power_of_scaleLower=log10(scaleLower);
%     power_of_scaleUpper=log10(scaleUpper);
%     scales=logspace(power_of_scaleLower,power_of_scaleUpper,numScales);
%     freqs=scal2frq(scales,my_wavelet_fun,DT);
    freqs=linspace(freqL,freqU,numScales);
    scales=fc./(freqs.*DT);
else
if(Flag==1) %scales linear spaced
    scaleLower=fc/(freqU*DT);
    scaleUpper=fc/(freqL*DT);
   scales=linspace(scaleLower,scaleUpper,numScales);
   freqs=scal2frq(scales,my_wavelet_fun,DT);
else
    disp('Error for Choosing scales or freqs linear spaced!!!');
end
end
Coeffs=cwt(sample,scales,my_wavelet_fun);
Coeffs_Struct=struct('cfs',abs(Coeffs),'scales',scales,'frequency',freqs);
end
