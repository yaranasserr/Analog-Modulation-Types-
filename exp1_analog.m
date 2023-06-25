clc; clear; 
% 1-Play your  file through Matlab 
[source,FS]=audioread('eric.wav'); 
NS= length(source);    %num of samples 
fprintf('the  file is playing\n'); 
% sound(source,FS);
% pause(8);

%spectrum of this signal 
source_freq= fftshift(fft(source));   
t=linspace(0,NS/FS,NS);     %audio time vector
F=linspace(-FS/2,FS/2,NS);   %audio freq Vector
f = FS/2*linspace(-1,1,NS); %for filter
figure('Name','Original Signal', 'NumberTitle','OFF');  
subplot(2, 1, 1);         plot(t, source); 
title('Original signal - Time domain'); 
subplot(2, 1, 2);         plot(F,abs(source_freq)/NS); 
title('Original signal - Frequency domain'); 

 
%2- ideal filter 
Filtered_signal = LP_Filter(source);

x= real(ifft(ifftshift(Filtered_signal))); 
Filtered_signal=abs(Filtered_signal);
figure('Name','Filtered Signal', 'NumberTitle','OFF');  
subplot(2, 1, 1);         plot(t,x); 
title('filtered signal - Time domain'); 
subplot(2, 1, 2);         plot(F,Filtered_signal/length(x)); 
title('filtered signal - Frequency domain'); 
fprintf('the Filtered Signal is playing\n');
sound(x,FS);
pause(8);

%generate modulation
Fc=100000;
FS_new=5*Fc;


msg_t=resample(x,FS_new,FS);
a=max(msg_t);
msg_f=abs(fftshift(fft(msg_t)));
NS2= length(msg_t);
t2=linspace(0,NS2/FS_new,NS2);
F2=linspace(-FS_new/2,FS_new/2,NS2);

%carrier for DSB-SC and DSB-TC
Carrier = (cos(2*pi*Fc*t2'));

%DSB_SC
SC_t=msg_t .* Carrier;
SC_f=abs(fftshift(fft(SC_t)));
figure('Name','Modulated SC Signal', 'NumberTitle','OFF');  
subplot(2, 1, 1);         plot(t2,SC_t); 
title('modulated signal DSB-SC - Time domain'); 
subplot(2, 1, 2);         plot(F2,SC_f/length(SC_t)); 
title('modulated signal DSB-SC - Frequency domain');

%DSB-TC
DC_bias=2*a;
TC_t=(DC_bias+msg_t) .* Carrier;
TC_f=abs(fftshift(fft(TC_t)));
figure('Name','Original Signal', 'NumberTitle','OFF');  
subplot(2, 1, 1);         plot(t2,TC_t); 
title('modulated TC signal DSB-TC - Time domain'); 
subplot(2, 1, 2);         plot(F2,TC_f/length(TC_t)); 
title('modulated signal DSB-TC - Frequency domain');

%envelop detector
env_SC=abs(hilbert(SC_t));
env_TC=abs(hilbert(TC_t))-2.*a;

%receiving both mod types
%DSB-SC
%envelop detector
rec_SC=resample(env_SC,FS,FS_new);
NS3= length(rec_SC);
t3=linspace(0,NS3/FS,NS3);    
F3=linspace(-FS/2,FS/2,NS3);
rec_SC_f=abs(fftshift(fft(rec_SC)));
figure('Name','Recieved SC Signal', 'NumberTitle','OFF');
subplot(2, 1, 1);         plot(t3,rec_SC); 
title('received signal DSB-SC - Time domain'); 
subplot(2, 1, 2);         plot(F3,rec_SC_f/length(rec_SC)); 
title('received signal DSB-SC - Frequency domain');
fprintf('the Recieved SC Signal Envelop is playing\n');
sound(rec_SC,FS);
pause(8);

%DSB-TC
rec_TC=resample(env_TC,FS,FS_new);
NS3= length(rec_TC);
t3=linspace(0,NS3/FS,NS3);    
F3=linspace(-FS/2,FS/2,NS3);
rec_TC_f=abs(fftshift(fft(rec_TC)));
figure('Name','Recieved TC Signal', 'NumberTitle','OFF');
subplot(2, 1, 1);         plot(t3,rec_TC); 
title('received signal DSB-TC - Time domain'); 
subplot(2, 1, 2);         plot(F3,rec_TC_f/length(rec_TC)); 
title('received signal DSB-TC - Frequency domain');
sound(rec_TC,FS);
pause(8);

%coherent detector of DSB-SC with diff SNR
rec_sc_0=awgn(SC_t,0,'measured');   %0 dB
rec_sc_10=awgn(SC_t,10,'measured'); %10 dB
rec_sc_30=awgn(SC_t,30,'measured'); %30 dB

%Mixer With Carrier
coherent=SC_t.*Carrier.*2; %without noise
SC_0=rec_sc_0.*Carrier.*2;
SC_10=rec_sc_10.*Carrier.*2;
SC_30=rec_sc_30.*Carrier.*2;

%resambling to original FS
coherent=resample(coherent,FS,FS_new);
SC_0=resample(SC_0,FS,FS_new);
SC_10=resample(SC_10,FS,FS_new);
SC_30=resample(SC_30,FS,FS_new);

%LP_Filter for original signal
Filtered_Coh=LP_Filter(coherent);
Filtered_0 = LP_Filter(SC_0);
Filtered_10 = LP_Filter(SC_10);
Filtered_30 = LP_Filter(SC_30);

x_coh=real(ifft(ifftshift(Filtered_Coh))); 
x0 = real(ifft(ifftshift(Filtered_0))); 
x10 = real(ifft(ifftshift(Filtered_10))); 
x30 = real(ifft(ifftshift(Filtered_30))); 

figure('Name','Recieved SC Coherent Signal', 'NumberTitle','OFF');
subplot(2,1,1); plot(t,x_coh(1:length(t))); title(' coherent received signal DSB-SC  - Time domain');
subplot(2,1,2); plot(F3,abs(Filtered_Coh)/length(x_coh)); xlim([-25000 25000]); title(' coherent received signal DSB-SC Frequency domain');
fprintf('the Recieved SC Signal Coherent is playing\n');
sound(x_coh,FS);
pause(8);

figure('Name','Recieved SC Coherent Signal With SNR =:', 'NumberTitle','OFF');
subplot(2,3,1); plot(t,x0(1:length(t))); title('0 dB - Time domain');
subplot(2,3,2); plot(t,x10(1:length(t))); title('10 dB - Time domain');
subplot(2,3,3); plot(t,x30(1:length(t))); title('30 dB - Time domain');
subplot(2,3,4); plot(F3,abs(Filtered_0)/length(x0)); xlim([-25000 25000]); title('0 dB - Frequency domain');
subplot(2,3,5); plot(F3,abs(Filtered_10)/length(x10)); xlim([-25000 25000]); title('10 dB - Frequency domain');
subplot(2,3,6); plot(F3,abs(Filtered_30)/length(x30)); xlim([-25000 25000]); title('30 dB - Frequency domain');
fprintf('the Recieved TC Signal With 0 dB SNR is playing\n');
sound(x0,FS);
pause(8);
fprintf('the Recieved TC Signal With 10 dB SNR is playing\n');
sound(x10,FS);
pause(8);
fprintf('the Recieved TC Signal With 30 dB SNR is playing\n');
sound(x30,FS);
pause(8);

%-------- Frequency Error calling functions
%Frequency:distortion +attenuation 
Frequency_error_SNR0 = FREQ(F3,NS,FS,FS_new,t2,rec_sc_0,'Time Domain Coherent Detection with Frequency error and SNR=0','Spectrum of  Coherent Detection with Frequency error and SNR=0');
Frequency_error_SNR10 = FREQ(F3,NS,FS,FS_new,t2,rec_sc_10,'Time Domain Coherent Detection with Frequency error and SNR=10','Spectrum of  Coherent Detection with Frequency error and SNR=10');
Frequency_error_SNR30 = FREQ(F3,NS,FS,FS_new,t2,rec_sc_30,'Time Domain Coherent Detection with Frequency error and SNR=30','Spectrum of  Coherent Detection with Frequency error and SNR=30');

%---------Phase Error
%Phase Error:attenuation
Phase_error_SNR0 = PH(F3,NS,FS,FS_new,t2,rec_sc_0,'Time Domain Coherent Detection with Phase error and SNR=0','Spectrum of Coherent Detection with Phase error and SNR=0');
Phase_error_SNR10 = PH(F3,NS,FS,FS_new,t2,rec_sc_10,'Time Domain Coherent Detection with Phase error and SNR=10','Spectrum of Coherent Detection with Phase error and SNR=10');
Phase_SNR30 = PH(F3,NS,FS,FS_new,t2,rec_sc_30,'Time Domain Coherent Detection with Phase error and SNR=30','Spectrum of Coherent Detection with Phase error and SNR=30');

%--------functions-----%
function f = LP_Filter(signal)
    NS = length(signal);
    FS = 48000;
    Filter = ones(NS,1); %rect signal ,amp=1 ,length =length of input signal
    source_freq= fftshift(fft(signal)); %transform signal to frequency domain to apply filter
    %Y = fix( X ) rounds each element of X to the nearest integer toward zero
    %to calculate number of wanted samples within cut off frequency
    point1 = fix((NS/2)-(4000*(NS/FS)));      %negative side
    point2 = fix((NS/2)+(4000*(NS/FS)));      %positive  side
    
    Filter([1:point1 point2+1:end]) = 0;      %unwanted samples , out of range areequal zero
    f = source_freq.*Filter;
end

function plotting(x3,y_frequ,timex,timey,timetitle_label,Freqtitle_label)
    figure()
    subplot(2,1,1); plot(timex,timey);title(timetitle_label);
    subplot(2,1,2) ;plot(x3,y_frequ);
    title(Freqtitle_label);	
end

function Freq= FREQ(F3,NS,FS,FS_new,t2,recieved,timetitle_label,Freqtitle_label) 

    Fc_error = 100100;
    Carrier_Freq_error = (cos(2*pi* Fc_error*t2'));
    Frequency_errored = recieved .* Carrier_Freq_error .*2;
    Frequency_errored=resample(Frequency_errored,FS,FS_new);
    Frequency_errored_Filtered = LP_Filter(Frequency_errored);
    Freq = real(ifft(ifftshift(Frequency_errored_Filtered))); 
    t_freq=linspace(0,NS/FS_new,NS);
    y_freq=abs(Frequency_errored_Filtered)/length(Freq);
    time_y=Freq(1:length(t_freq));
    plotting(F3,y_freq,t_freq,time_y,timetitle_label,Freqtitle_label)
    fprintf('distorted + attenuated\n');
    sound(Freq,FS);
    pause(8)
    
end
function Phase= PH(F3,NS,FS,FS_new,t2,recieved,timetitle_label,Freqtitle_label) 

    Fc = 100000;
    Carrier_phased = (cos(2*pi*Fc*t2'+pi/9));
    Phased_errored = recieved .* Carrier_phased .*2;
    Phased_errored=resample(Phased_errored,FS,FS_new);
    Phased_errored_Filtered = LP_Filter(Phased_errored);
    Phase = real(ifft(ifftshift(Phased_errored_Filtered))); 
    
    t_freq=linspace(0,NS/FS_new,NS);
    y_freq=abs(Phased_errored_Filtered)/length( Phase);
    
    time_y= Phase(1:length(t_freq));
    
    plotting(F3,y_freq,t_freq,time_y,timetitle_label,Freqtitle_label)
    fprintf('attenuated\n');
    sound(Phase,FS);
    pause(8)   
end

