%% -----------------Reading the auidios-----------------------%%

[y1,Fs] = audioread('Short_BBCArabic2.wav');

[y2,Fs] = audioread('Short_FM9090.wav');

[y3,Fs] = audioread('Short_QuranPalestine.wav');

[y4,Fs] = audioread('Short_RussianVoice.wav');

[y5,Fs] = audioread('Short_SkyNewsArabia.wav');

[y6,Fs] = audioread('Short_WRNArabic.wav');

%% ------------------Adding the two channels of the signals into one channel---------------%%

x1=y1(:,1)+y1(:,2);

x2=y2(:,1)+y2(:,2);

x3=y3(:,1)+y3(:,2);

x4=y4(:,1)+y4(:,2);

x5=y5(:,1)+y5(:,2);

x6=y6(:,1)+y6(:,2);

%% ------------------ making all signals with the same length then add interpulation -----------------------%%

sig1 = interp(padarray(x1,[896 0],0,'post'),12);

sig2 = interp(padarray(x2,[43904 0],0,'post'),12);

sig3 = interp(padarray(x3,  [2240 0],0,'post'),12);

% sig4 = interp(padarray(x4,[38080 0],0,'post'),12);

% sig5 = interp(padarray(x5,[29568 0],0,'post'),12);

% sig6= interp(x6,10);

%%-------------------- 1st modulated signal----------------------%%

t1=0:1/(12*Fs):(length(sig1)-1)/(12*Fs);

carrier=cos(2*pi*100000*t1);

carrier_trans=transpose(carrier);

modulated_sig1=sig1.*carrier_trans;

modulated_sig1_fd=fft(modulated_sig1);

fr=linspace(-6*Fs,6*Fs,length(modulated_sig1_fd));

modulated_sig1_fd1=fftshift(abs(modulated_sig1_fd));

%%--------------------- 2nd modulated signal---------------------%%

t2=0:1/(12*Fs):(length(sig2)-1)/(12*Fs);

carrier2=cos(2*pi*160000*t2);

carrier_trans2=transpose(carrier2);

modulated_sig2=sig2.*carrier_trans2;

modulated_sig2_fd=fft(modulated_sig2);

fr2=linspace(-6*Fs,6*Fs,length(modulated_sig2_fd));

modulated_sig2_fd1=modulated_sig1_fd1+fftshift(abs(modulated_sig2_fd));

%%--------------------- 3rd modulated signal----------------------%%

t3=0:1/(12*Fs):(length(sig3)-1)/(12*Fs);

carrier3=cos(2*pi*220000*t3);

carrier_trans3=transpose(carrier3);

modulated_sig3=sig3.*carrier_trans3;

modulated_sig3_fd=fft(modulated_sig3);

fr3=linspace(-6*Fs,6*Fs,length(modulated_sig3_fd));

figure

modulated_sig3_fd1=plot(fr3,modulated_sig2_fd1+fftshift(abs(modulated_sig3_fd)));

title('Figure 1: The spectrum of the output of the transmitter')

ylabel('Amplitude')

xlabel('Frequency')

%% ----------------------------spectrums of 6 the signals-----------------------------%%

 subplot(3,2,1); 

 plot(fftshift(abs(fft(x1))))

 title(' Spectrum of sig 1')

  subplot(3,2,2);

  plot(fftshift(abs(fft(x2))))

  title(' Spectrum of sig 2')

  subplot(3,2,3);

  plot(fftshift(abs(fft(x3))))

  title(' Spectrum of sig 3')

  subplot(3,2,4);

  plot(fftshift(abs(fft(x4))))

 title(' Spectrum of sig 4')

  subplot(3,2,5);

  plot(fftshift(abs(fft(x5))))

  title(' Spectrum of sig 5')

  subplot(3,2,6);

  plot(fftshift(abs(fft(x6))))

  title(' Spectrum of sig 6')  

%% ------------------(RF)filter of modulated_signal_1-------------------------------%%  

 multiplexed_time=modulated_sig1+modulated_sig2+modulated_sig3;

 RF = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',70000/(6*Fs),89340/(6*Fs),108500/(6*Fs),130000/(6*Fs),60,1,60);

 Hd = design(RF,'equiripple');

 FILTER1=filter(Hd,multiplexed_time);

 order_filter=order(Hd);

 figure

 plot(fr,fftshift(abs(fft(FILTER1))));

 title('Figure 2: the output of the RF filter 1 (before the mixer)')

ylabel('Amplitude')

xlabel('Frequency')

 %% -----------------(RF)filter of modulated_signal_2-----------------------%%

 RF2 = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',70000/(6*Fs),149000/(6*Fs),170000/(6*Fs),180000/(6*Fs),60,1,60);

 Hd2 = design(RF2,'equiripple');

 FILTER2=filter(Hd2,multiplexed_time);

 order_filter2=order(Hd2);

 figure 

plot(fr2,fftshift(abs(fft(FILTER2))));

 title('Figure 2: the output of the RF filter 2 (before the mixer)')

 ylabel('Amplitude')

xlabel('Frequency')

 %% -----------------(RF)filter of modulated_signal_3------------------------%%

 RF3 = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',70000/(6*Fs),210000/(6*Fs),230000/(6*Fs),240000/(6*Fs),60,1,60);

 Hd3 = design(RF3,'equiripple');

 FILTER3=filter(Hd3,multiplexed_time);

 order_filter3=order(Hd3);

 figure

 plot(fr3,fftshift(abs(fft(FILTER3))));

 title('Figure 2: the output of the RF filter 3 (before the mixer)')

 ylabel('Amplitude')

xlabel('Frequency')

 %% -----------------mix band pass filter of signal_1 with carrier of Wc+WiF----------------%%

 % t1_FILTER=0:1/(12*Fs):(length( multiplexed_time)-1)/(12*Fs); in case of no RF

 t1_FILTER=0:1/(12*Fs):(length(FILTER1)-1)/(12*Fs);

carrier_FILTER=cos(2*pi*130000*t1_FILTER);

carrier_FILTER_trans=transpose(carrier_FILTER);

% modulated_FILTER1= multiplexed_time.*carrier_FILTER1_trans; in case of no RF

modulated_FILTER1=FILTER1.*carrier_FILTER_trans;

modulated_FILTER1_fd=fft(modulated_FILTER1);

f_filter1=linspace(-6*Fs,6*Fs,length(modulated_FILTER1_fd));

figure

plot(f_filter1,fftshift(abs(modulated_FILTER1_fd)));

 title('Figure 3: The output of the mixer for signal 1')

 ylabel('Amplitude')

xlabel('Frequency')

%% ----------------- mix band pass filter of signal_2 with carrier of Wc+WiF----------------%%

  % t2_FILTER=0:1/(12*Fs):(length( multiplexed_time)-1)/(12*Fs); in case of no RF

t2_FILTER=0:1/(12*Fs):(length(FILTER2)-1)/(12*Fs);

carrier_FILTER2=cos(2*pi*190000*t2_FILTER);

carrier_FILTER2_trans=transpose(carrier_FILTER2);

% modulated_FILTER2= multiplexed_time.*carrier_FILTER2_trans; in case of no RF

modulated_FILTER2=FILTER2.*carrier_FILTER2_trans;

modulated_FILTER2_fd=fft(modulated_FILTER2);

f_filter2=linspace(-6*Fs,6*Fs,length(modulated_FILTER2_fd));

figure

plot(f_filter2,fftshift(abs(modulated_FILTER2_fd)));

 title('Figure 3: The output of the mixer for signal 2')

 ylabel('Amplitude')

xlabel('Frequency')

%% ---------------- mix band pass filter of signal_3 with carrier of Wc+WiF-----------------%%

   % t3_FILTER=0:1/(12*Fs):(length( multiplexed_time)-1)/(12*Fs); in case of no RF

t3_FILTER=0:1/(12*Fs):(length(FILTER3)-1)/(12*Fs);

carrier_FILTER3=cos(2*pi*190000*t3_FILTER);

carrier_FILTER3_trans=transpose(carrier_FILTER3);

% modulated_FILTER3= multiplexed_time.*carrier_FILTER3_trans; in case of no RF

modulated_FILTER3=FILTER3.*carrier_FILTER3_trans;

modulated_FILTER3_fd=fft(modulated_FILTER3);

f_filter3=linspace(-6*Fs,6*Fs,length(modulated_FILTER3_fd));

figure

plot(f_filter3,fftshift(abs(modulated_FILTER3_fd)));

 title('Figure 3: The output of the mixer for signal 3')

 ylabel('Amplitude')

xlabel('Frequency')

%% -----------------------(IF) band pass filter2 of signal 1--------------------%% 

 RF_1 = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',1/(6*Fs),10000/(6*Fs),29000/(6*Fs),50000/(6*Fs),60,1,60);

 Hd_1 = design(RF_1,'equiripple');

 FILTER_WIF_1=filter(Hd_1,modulated_FILTER1);

 f_s_filter1=linspace(-6*Fs,6*Fs,length(fft( FILTER_WIF_1)));

  order_filter_1=order(Hd_1);

 figure

 plot(f_s_filter1,fftshift(abs(fft( FILTER_WIF_1))))

 title('Figure 4: Output of the IF filter1')

 ylabel('Amplitude')

xlabel('Frequency')

%% ----------------------(IF) band pass filter2 of signal 2 --------------------%%

 RF_2 = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',1/(6*Fs),10000/(6*Fs),29000/(6*Fs),50000/(6*Fs),60,1,60);

 Hd_2 = design(RF_2,'equiripple');

 FILTER_WIF_2=filter(Hd_2,modulated_FILTER2);

 f_s_filter2=linspace(-6*Fs,6*Fs,length(fft( FILTER_WIF_2)));

   order_filter_2=order(Hd_2);

 figure

 plot(f_s_filter2,fftshift(abs(fft( FILTER_WIF_2))))

 title('Figure 4: Output of the IF filter2')

 ylabel('Amplitude')

xlabel('Frequency')

%% ----------------------(IF) band pass filter2 of signal 3 ---------------------%%

 RF_3 = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',1/(6*Fs),10000/(6*Fs),29000/(6*Fs),50000/(6*Fs),60,1,60);

 Hd_3 = design(RF_3,'equiripple');

 FILTER_WIF_3=filter(Hd_3,modulated_FILTER3);

 f_s_filter3=linspace(-6*Fs,6*Fs,length(fft( FILTER_WIF_3)));

   order_filter_3=order(Hd_3);

 figure

 plot(f_s_filter3,fftshift(abs(fft( FILTER_WIF_3))))

 title('Figure 4: Output of the IF filter3')

 ylabel('Amplitude')

xlabel('Frequency')

 %% -------------------- mix the signal_1 after second band pass filter with cos WIF ----------------------%% 

t1_second_FILTER=0:1/(12*Fs):(length(FILTER_WIF_1)-1)/(12*Fs);

carrier_second_FILTER=cos(2*pi*30000*t1_second_FILTER);

carrier_second_FILTER_trans=transpose(carrier_second_FILTER);

modulated_second_FILTER1=FILTER_WIF_1.*carrier_second_FILTER_trans;

f_second_filter1=linspace(-6*Fs,6*Fs,length(modulated_second_FILTER1));

figure

plot(f_second_filter1,fftshift(abs(fft(modulated_second_FILTER1))))

title('Figure 5: Output of the mixer for signal 1 (before the LPF)')

ylabel('Amplitude')

xlabel('Frequency')

%% --------------------- low pass filter for signal 1 ------------------------%%

d1=fdesign.lowpass('Fp,Fst,Ap,Ast',20000/(6*Fs),25000/(6*Fs),1,60);

Hd_LPF1 = design(d1,'equiripple');

LPF1=filter(Hd_LPF1,modulated_second_FILTER1);

f_LPF1=linspace(-6*Fs,6*Fs,length(modulated_second_FILTER1));

 order_filter_lPF1=order(Hd_LPF1);

figure

plot(f_LPF1,fftshift(abs(fft(LPF1))));

title('Figure 6: Output of the LPF1)') 

ylabel('Amplitude')

xlabel('Frequency')

lpf=downsample(LPF1,12);

sound(lpf,Fs)

%% --------------------- mix the signal_2 after second band pass filter with cos WIF -------------%%

 t2_second_FILTER=0:1/(12*Fs):(length(FILTER_WIF_2)-1)/(12*Fs);

carrier_second_FILTER2=cos(2*pi*30000*t2_second_FILTER);

carrier_second_FILTER_trans2=transpose(carrier_second_FILTER2);

modulated_second_FILTER2=FILTER_WIF_2.*carrier_second_FILTER_trans2;

f_second_filter2=linspace(-6*Fs,6*Fs,length(modulated_second_FILTER2));

figure

plot(f_second_filter2,fftshift(abs(fft(modulated_second_FILTER2))))

title('Figure 5: Output of the mixer for signal 2 (before the LPF)') 

ylabel('Amplitude')

xlabel('Frequency')

%% --------------------- low pass filter for signal 2 ------------------------%%

d2=fdesign.lowpass('Fp,Fst,Ap,Ast',20000/(6*Fs),25000/(6*Fs),1,60);

Hd_LPF2 = design(d2,'equiripple');

LPF2=filter(Hd_LPF2,modulated_second_FILTER2);

f_LPF2=linspace(-6*Fs,6*Fs,length(modulated_second_FILTER2));

 order_filter_lPF2=order(Hd_LPF2);

figure

plot(f_LPF2,fftshift(abs(fft(LPF2))));

title('Figure 6: Output of the LPF2)')

ylabel('Amplitude')

xlabel('Frequency')

lpf2=downsample(LPF2,12);

%sound(lpf2,Fs)

%% -------------------- mix the signal_3 after second band pass filter with cos WIF -------------%%

 t3_second_FILTER=0:1/(12*Fs):(length(FILTER_WIF_3)-1)/(12*Fs);

carrier_second_FILTER3=cos(2*pi*30000*t3_second_FILTER);

carrier_second_FILTER_trans3=transpose(carrier_second_FILTER3);

modulated_second_FILTER3=FILTER_WIF_3.*carrier_second_FILTER_trans3;

f_second_filter3=linspace(-6*Fs,6*Fs,length(modulated_second_FILTER3));

figure

plot(f_second_filter3,fftshift(abs(fft(modulated_second_FILTER3))))

title('Figure 5: Output of the mixer for signal 3 (before the LPF)')

ylabel('Amplitude')

xlabel('Frequency')

%% --------------------- low pass filter for signal 3 ------------------------%%

d3=fdesign.lowpass('Fp,Fst,Ap,Ast',20000/(6*Fs),25000/(6*Fs),1,60);

Hd_LPF3 = design(d3,'equiripple');

LPF3=filter(Hd_LPF3,modulated_second_FILTER3);

f_LPF3=linspace(-6*Fs,6*Fs,length(modulated_second_FILTER3));

order_filter_lPF3=order(Hd_LPF3);

figure

plot(f_LPF3,fftshift(abs(fft(LPF3))));

title('Figure 6: Output of the LPF3)')

ylabel('Amplitude')

xlabel('Frequency')

lpf3=downsample(LPF3,12);

%sound(lpf3,Fs)





