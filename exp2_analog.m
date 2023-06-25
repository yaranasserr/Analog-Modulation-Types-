clc; clear; 
% 1-Play your sound file through Matlab  
[source,FS]=audioread('eric.wav'); 
NS= length(source);    %num of samples 
fprintf('the sound file is playing\n'); 
% sound(source,FS);
% pause(8);

%spectrum of this signal 
source_freq= fftshift(fft(source));   
t=linspace(0,NS/FS,NS);     %audio time vector
F=linspace(-FS/2,FS/2,NS);   %audio freq Vector
f = FS/2*linspace(-1,1,NS); %for filter
figure(1);  
subplot(2, 1, 1);         plot(t, source); 
title('Original signal - Time domain'); 
subplot(2, 1, 2);         plot(F,abs(source_freq)/NS); 
title('Original signal - Frequency domain'); 

%2- ideal filter 
Filtered_signal = LP_Filter(source);        %LPF BW=4000
x= real(ifft(ifftshift(Filtered_signal))); 
Filtered_signal=abs(Filtered_signal);
figure();  
subplot(2, 1, 1);         plot(t,x); 
title('filtered signal - Time domain'); 
subplot(2, 1, 2);         plot(F,Filtered_signal/length(x)); 
title('filtered signal - Frequency domain'); 
% sound(x,FS);
% pause(8);

%Generate a DSB-SC modulated signal and plot its spectrum
Fc=100000;           %carrier freq
FS_new=5*Fc;         %FS_new : sampling freq of thr filtered audio signal after resample
msg_t=resample(x,FS_new,FS);     
a=max(msg_t);
msg_f=abs(fftshift(fft(msg_t)));
NS2= length(msg_t);
t2=linspace(0,NS2/FS_new,NS2);          %msg time vector after resample
F2=linspace(-FS_new/2,FS_new/2,NS2);    %msg freq vector after resample

%carrier for DSB-SC and DSB-TC
Carrier = cos(2*pi*Fc*t2');

%DSB_SC
SC_t=msg_t .* Carrier;           %modulated signal of DSB-SC in time
SC_f=fftshift(fft(SC_t));        %modulated signal of DSB-SC in freq
figure(3);
subplot(2, 1, 1);         plot(t2,SC_t); 
title('modulated signal DSB-SC - Time domain'); 
subplot(2, 1, 2);         plot(F2,abs(SC_f)/length(SC_t)); 
title('modulated signal DSB-SC - Frequency domain');

%5- Obtain the SSB-LSB
LSB_f = BP_Filter(SC_t);              %SSB modulated signal in freq
LSB_t=real(ifft(ifftshift(LSB_f)));     %SSB modulated signal in time
figure(4)
subplot(2, 1, 1);         plot(t2,LSB_t); 
title('modulated signal of LSB-SC in Time domain'); 
subplot(2, 1, 2);         plot(F2,abs(LSB_f)/length(LSB_t));
title('modulated signal of LSB-SC in Frequency domain');

%Use coherent detection with no noise interference to get the received signal
Demodulated_t =  LSB_t .* Carrier .* 2; 
Demodulated_t=resample(Demodulated_t,FS,FS_new);
Demodulated_Filtered_f = LP_Filter(Demodulated_t);                      %received LSB signal after deacreasing freq
Demodulated_Filtered_t = real(ifft(ifftshift(Demodulated_Filtered_f))); %received LSB signal after deacreasing freq


figure()  
subplot(2, 1, 1);         plot(t,Demodulated_Filtered_t(1:length(t)));
title('Demodulated signal of Demodulated SSB-SC in Time domain');
subplot(2, 1, 2);         plot(F,abs(Demodulated_Filtered_f(1:length(F)))/length(Demodulated_Filtered_t)); 
title('Demodulated signal of LSB-SC in Frequency domain');
% sound(Demodulated_Filtered_t,FS);
% pause(8);


%7.1-Obtaining the SSB-LSB using a practical 4th order Butterworth filter
[b a] = butter(4,[2*(Fc-4000)/FS_new  2*Fc/FS_new]);   %BPF practical 4th order Butterworth filter
butter_filtered_t = filter(b,a,SC_t);                  %the filteration of DSB-SC to LSB-SC in time
butter_filtered_f = fftshift(fft(butter_filtered_t));%LBS in freq
figure('Name','') 
subplot(2, 1, 1);         plot(t2,butter_filtered_t(1:length(t2)));
title('Modulated signal of Demodulated SSB-SC by Butter in Time domain');
subplot(2, 1, 2);         plot(F2,abs(butter_filtered_f(1:length(F2)))/length(butter_filtered_t)); 
title('Modulated signal of LSB-SC in Frequency domain');

%7.2- coherent detection with no noise by practical 4th Butterworth filter
Demodulated_butter_t =  butter_filtered_t .* Carrier .* 2; 
%Demodulated_butter_t = resample(Demodulated_butter_t,FS,FS_new);
[b,a]=butter(4,4000/(FS_new/2));
Demodulated_butter_t = filter(b,a,Demodulated_butter_t);
Demodulated_butter_f=abs(fftshift(fft(Demodulated_butter_t)));
rec_butter_t = resample(Demodulated_butter_t,FS,FS_new);
rec_butter_f=fftshift(fft(rec_butter_t));
%Demodulated_butter_t = real(ifft(ifftshift(Demodulated_butter_f))); 

figure()  
subplot(2, 1, 1);         plot(t,rec_butter_t(1:length(t)));
title('Demodulated signal of Demodulated SSB-SC Buttered in Time domain');
subplot(2, 1, 2);         plot(F,abs(rec_butter_f(1:length(F)))/length(Demodulated_butter_t)); 
title('Demodulated signal of LSB-SC Buttered in Frequency domain');
% sound(Demodulated_butter_t,FS);
% pause(8);

%9- Noise SNR rakam 8
LSB_t_0 = awgn(LSB_t,0);   
LSB_t_10 = awgn(LSB_t,10);   
LSB_t_30 = awgn(LSB_t,30);

Dem_LSB_t_0 = LSB_t_0 .* Carrier .* 2;
Dem_LSB_t_10 = LSB_t_10 .* Carrier .* 2;
Dem_LSB_t_30 = LSB_t_30 .* Carrier .* 2;

Dem_LSB_t_0=resample(Dem_LSB_t_0,FS,FS_new);
Dem_LSB_t_10=resample(Dem_LSB_t_10,FS,FS_new);
Dem_LSB_t_30=resample(Dem_LSB_t_30,FS,FS_new);

Dem_LSB_f_0 = LP_Filter(Dem_LSB_t_0);
Dem_LSB_f_10 = LP_Filter(Dem_LSB_t_10);
Dem_LSB_f_30 = LP_Filter(Dem_LSB_t_30);

x0=real(ifft(ifftshift(Dem_LSB_f_0)));
x10=real(ifft(ifftshift(Dem_LSB_f_10)));
x30=real(ifft(ifftshift(Dem_LSB_f_30)));

figure()
subplot(2,3,1); plot(t,x0(1:length(t))); title('SSB-SC SNR = 0 dB - Time domain');
subplot(2,3,2); plot(t,x10(1:length(t))); title('SSB-SC SNR = 10 dB - Time domain');
subplot(2,3,3); plot(t,x30(1:length(t))); title('SSB-SC SNR = 30 dB - Time domain');
subplot(2,3,4); plot(F,abs(Dem_LSB_f_0(1:length(F))/length(x0))); xlim([-25000 25000]); title('SSB-SC SNR = 0 dB - Frequency domain');
subplot(2,3,5); plot(F,abs(Dem_LSB_f_10(1:length(F))/length(x10))); xlim([-25000 25000]); title('SSB-SC SNR = 10 dB - Frequency domain');
subplot(2,3,6); plot(F,abs(Dem_LSB_f_30(1:length(F))/length(x30))); xlim([-25000 25000]); title('SSB-SC SNR = 30 dB - Frequency domain');

% sound(x0,FS);
% pause(8);
% sound(x10,FS);
% pause(8);
% sound(x30,FS);
% pause(8);

%9.1- generate SSB-TC
DC_bias = 2 * max(msg_t);   %DC bias added to messege signal
SC_t=(msg_t+DC_bias).* Carrier;       %equation of s(t) DSB-TC
SC_f=fftshift(fft(SC_t));
LSB = BP_Filter(SC_t);
LSB_t=real(ifft(ifftshift(LSB)));
env_TC=abs(hilbert(LSB_t)) - 2 .* DC_bias;   %demodulated signal
rec_TC=resample(env_TC,FS,FS_new);
rec_TC_f=fftshift(fft(rec_TC));
figure();
subplot(2, 1, 1);         plot(t,rec_TC(1:length(t))); 
title('received signal SSB-TC - Time domain'); 
subplot(2, 1, 2);         plot(F,abs(rec_TC_f(1:length(F))/length(rec_TC))); 
title('received signal SSB-TC - Frequency domain');
sound(rec_TC,FS);
pause(8);

function f = LP_Filter(signal)
    NS = length(signal);
    FS = 48000;
    Filter = ones(NS,1);
    source_freq= fftshift(fft(signal));
    point1 = fix((NS/2)-(4000*(NS/FS)));      %cutoff freq front 
    point2 = fix((NS/2)+(4000*(NS/FS)));      %cutoff freq back 
    Filter([1:point1 point2+1:end]) = 0;      %Equation of filter 
    f = source_freq.*Filter;
end

function X = BP_Filter(signal)
    signal=fftshift(fft(signal));
    NS=length(signal);
    BPF=ones(NS,1);
    Fc = 100000;
    FS_new = Fc*5;
    BW=4000;
    p1=fix((NS/2)-(Fc*(NS/FS_new)));       %cutoff freq 1
    p2=fix((NS/2)+(Fc*(NS/FS_new)));       %cutoff freq 4
    p3=fix((NS/2)-((Fc-BW)*(NS/FS_new)));  %cutoff freq 2
    p4=fix((NS/2)+((Fc-BW)*(NS/FS_new)));  %cutoff freq 3
    BPF([1:p1 p3:p4 p2+3:end])=0;          %ideal BPF
    LSB=signal .* BPF .* 2;
    X = LSB;
end
