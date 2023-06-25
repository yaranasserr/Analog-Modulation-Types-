
%coherent detector of DSB-SC
% demodulated=SC_t.*Carrier;
% demodulated=demodulated.*2;
% dsbsc=resample(demodulated,FS,FS_new);
% demodulated_F = abs(fftshift(fft(demodulated)));
% Filtered_signal = LP_Filter(dsbsc);
% x = real(ifft(ifftshift(Filtered_signal))); 
% Filtered_signal=abs(Filtered_signal);
% F_SC=linspace(-FS_new/2,FS_new/2,length(demodulated_F));  
% figure();
% subplot(2, 1, 1);         plot(t2,x); 
% title('received signal DSB-SC - Time domain'); 
% subplot(2, 1, 2);        plot(F_SC,abs(Filtered_signal)/length(demodulated)); 
% title('received signal DSB-SC - Frequency domain');
% sound(demodulated,FS);
% pause(8);