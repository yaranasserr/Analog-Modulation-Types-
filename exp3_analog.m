clc; clear; 
% 1-Play your sound file through Matlab 
[source,FS]=audioread('eric.wav'); 
NS= length(source);    %num of samples 
fprintf('the sound file is playing\n'); 
sound(source,FS);
pause(8);

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
Filtered_signal = LP_Filter(source);
x= real(ifft(ifftshift(Filtered_signal))); 
%signal in time domain after filter
figure('Name','Filtered Signal');  
subplot(2, 1, 1);         plot(t,x); 
title('filtered signal - Time domain'); 
subplot(2, 1, 2);         plot(F,abs(Filtered_signal)/length(x)); 
title('filtered signal - Frequency domain'); 
fprintf('the filtered sound file is playing\n'); 
% sound(x,FS);
% pause(8);

%generate modulation
Fc=100000;
FS_new=5*Fc;
%resampling the msg to 500000 to prevent overlapping between signal and carrier  after applying
%modulation 
msg_t=resample(x,FS_new,FS);
msg_f=fftshift(fft(msg_t));
NS2= length(msg_t);
%new values for plotting after resampling
t2=linspace(0,NS2/FS_new,NS2);    
F2=linspace(-FS_new/2,FS_new/2,NS2);

Carrier = cos(2*pi*Fc*t2');
kf = 1; %kf change the loudness of the msg , higher kf higher sound , kf=1 no change
A = 1;
% the original equation for nbfm without approximation
modulated_t = cos(2*pi*Fc*t2' + kf*cumsum(msg_t));
modulated_f=fftshift(fft(modulated_t));
figure('Name','Modulated Signal');  
subplot(2,1,1); plot(t2,modulated_t);
title('Modulated Signal - Time Domain')
subplot(2,1,2); plot(F2,abs(modulated_f)/length(modulated_t));
title('Modulated Signal - Frequency Domain')

%demodulate
%first differentiate to remove cumsum(integration effect)
drev = diff(modulated_t/kf);
%use envelop detector to demodulate the signal after differntiating , and
%subtract dc value(mean of signal)
envSignal=abs(hilbert(drev)) - mean(abs(hilbert(drev)));
%to sound the signal resample to old FS (original)
msg=real(resample(envSignal,FS,FS_new));
msg_f=abs(fftshift(fft(msg)));

figure('Name','DeModulated Signal');
subplot(2,1,1); plot(t, msg);
title('demodulated Signal - Time Domain')
subplot(2,1,2); plot(F,abs(msg_f)/length(msg));
title('demodulated Signal - Frequency Domain')
fprintf('the demodulated sound file is playing\n'); 
sound(msg,FS)

function f = LP_Filter(signal)
    NS = length(signal);
    FS = 48000;
    Filter = ones(NS,1);
    source_freq= fftshift(fft(signal));
    point1 = fix((NS/2)-(4000*(NS/FS)));      %cutoff freq front 
    point2 = fix((NS/2)+(4000*(NS/FS))); 
    Filter([1:point1 point2+1:end]) = 0;      %Equation of filter 
    f = source_freq.*Filter;
end

% demodulated_filter_f = LP_Filter(demodulated_t);
% demodulated_filter_t = real(ifft(ifftshift(demodulated_filter_f)));
% demodulated_t = modulated_t .* Carriershift;
% demodulated_t = resample(demodulated_t,FS,FS_new);
% demodulated_filter_f = LP_Filter(demodulated_t);
% demodulated_filter_t = real(ifft(ifftshift(demodulated_filter_f)));
% %demodulated_filter_t = diff(demodulated_filter_t);
% demodulated_filter_f= fftshift(fft(demodulated_filter_t));
% % demodulated_filter_t = demodulated_filter_t(100:length(demodulated_filter_t)-200);
% recievedSignal = resample(recievedSignal,FS,FS_new);
