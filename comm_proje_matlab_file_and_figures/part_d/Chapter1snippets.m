% Chapter 1 Code Snippets
clear all;close all;clc;
% Table 1.1, step (a)
% Create an audiorecorder object with 8000 sps and
% a single channel (mono); view its properties:
Fs=8000;
recObj = audiorecorder(Fs, 16, 1);
get(recObj)

% Collect a sample of your speech with a microphone;
% and, plot the signal data:

% Record your voice for 2 seconds. Use display and 
% pause as aids to control the start of the recording.
Trec = 0.5; %2 second record time
disp('Press Enter to start recording.')
pause;%wait for keystroke
disp('Recording.')
recordblocking(recObj, Trec);
disp('End of Recording.')

% Store data in double-precision array.
myRecording = getaudiodata(recObj);

%%%%%%%%%%%%%%%%%%%%%%%%%%% VOICE SIGNAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the waveform.
t = (0:length(myRecording)-1)/Fs; %sample times (sec)
figure(1)
subplot(3,1,1)
plot(t,myRecording);
title('Voice Signal');
xlabel('Time');
ylabel('Amplitude');

%%%%%%%%%%%%%%%%%%%%%%%%%%% SAMPLED SIGNAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (Exsampl e.m) 
% Example of sampl ing , quant ization , and zero- order hold 
%clear;clf; 
td= 0.002;      % o r i g i nal s amp l ing r a t e 5 0 0 H z
t = [0.002:td:8];  % t ime int e rva l o f 1 s e c ond
xsig = getaudiodata(recObj);
Lsig=length(xsig);  % 1 H z + 3 H z s i nu s o i d s

ts=0.02;    %new s amp l i ng r a t e = 5 0 H z .
Nfactor = ts/td;    % s end the s i gn a l through a 1 6 - l eve l uni f o rm quant i z e r


[s_out,sq_out,sqh_out,Delta,SQNR] = sampandquant(xsig,16,td,ts);
% r e c e ive 3 s i gnal s :
% 1 . s amp l e d s i gnal s_o u t
% 2 . s amp l ed and qu an t i z e d s i gn a l s q_out
% 3 . s amp l ed , quant i z ed , and z e r o - order h o l d s i gnal sqh_o u t

% calculate the Fourier transforms
Lfft=2^ceil(log2(Lsig)+1) ;
Fmax=1 / ( 2 * td) ;
Faxis=linspace ( -Fmax, Fmax , Lfft ) ;
Xsig=fftshift ( fft (xsig , Lfft ) ) ;
S_out=fftshift(fft(s_out,Lfft)) ;


% Examples of sampling and reconstruction using
% a ) ideal impulse train through LPF
% b) flat top pulse reconstruction through LPF
% plot the original signal and the sample signals in time
% and frequency domain


subplot (312) ; sfig1a=plot ( t , xsig , ' k ' ) ;
hold on ; sfig1b=plot ( t , s_out ( 1 : Lsig) , ' b ' ) ; hold off ;
set (sfig1a,'Linewidth',2); set ( sfig1b , 'Linewidth' , 2.) ;
xlabel ( ' time ( sec ) ' ) ;
title ( ' Signal and its uniform samples ' ) ;

%%%%%%%%%%%%%%%%%%%%%%%%%% QUANTIZED SIGNAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=16;
L=2^n;

vmax=1;
vmin=-vmax;
del=(vmax-vmin)/L;

part=vmin:del:vmax;
code=vmin-(del/2):del:vmax+(del/2);
[ind,q]=quantiz(myRecording,part,code);

l1=length(ind);
l2=length(q);

for i=1:l1
    if(ind(i)~=0)
        ind(i)=ind(i)-1;
    end
    i=i+1;
end

for i=1:l2
    if(q(i)==vmin-(del/2))
        q(i)=vmin+(del/2);
    end
end
subplot(313)
stem(q);
grid on;
title('Quantized Signal');

figure(2)
%%%%%%%%%%%%%%%%%%%%%% SIGNAL AND ITS 16 LEVEL PCM %%%%%%%%%%%%%%%%%%%%%%%%
td=1/Fs;
xsig=myRecording; %ses sinyali
Lsig=length(xsig);
Lfft=2^ceil(log2(Lsig)+1);
Xsig=fftshift(fft(xsig,Lfft));
Fmax=1/(2*td);
Faxis=linspace(-Fmax,Fmax,Lfft);
ts=0.02 %new sampling rate = 50 Hz
Nfact=ts/td;
%send the signal through a 16-level uniform quantizer
[s_out,sq_out,sqh_out1,Delta,SQNR]= sampandquant(xsig,16,td,ts);
%obtained the PCM signal which is
%   -sampled,quantized,and zero-order hold signal sqh_out
%plot the original signal and the PCM signal in time domain

subplot(3,1,1)
sfig1=plot(t,xsig,'k',t,sqh_out1(1:Lsig),'b');
set(sfig1,'Linewidth',2);
title('Signal and its 16-level PCM signal');
xlabel('time(sec.)')

%%%%%%%%%%%%%%%%%%%%%%%%%%% ENCODED SIGNAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
code=de2bi(ind,'left-msb');
k=1;
for i=1:l1
    for j=1:n
        coded(k)=code(i,j);
        j=j+1;
        k=k+1;
    end
    i=i+1;
end
subplot(3,1,2)
grid on;
stairs(coded);
axis([0 4000 -0.25 1.25]);
title('Encoded Signal')

%%%%%%%%%%%%%%%%%%%%%%%%% DEMODULATED SIGNAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qunt=reshape(coded,n,length(coded)/n);
index=bi2de(qunt','left-msb');
q=del*index+vmin+(del/2);
subplot(3,1,3)
grid on;
plot(q);
title('Demodulated Signal');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(3)
% % p l o t the i deal ly reconstructed s i gnal in t ime
% % and frequency domain
% 
% % c a l culate the recons tructed s i gnal f r om i deal s ampl ing and
% % ideal LPF
% % Maximum LPF bandwidth equal s to BW= f l oor ( ( L f f t /Nfactor ) / 2 ) ;
% 
% BW= 10 ; %Bandwidth i s no larger than l0Hz .
% H_1pf = zeros(1,Lfft) ; H_1pf(Lfft/2 - BW : Lfft / 2 +BW- 1 ) = 1 ; % i deal LPF
% S_recv=Nfactor*S_out.*H_1pf ; % ideal f i l ter ing
% s_recv= real(ifft(fftshift( S_recv ) ) ) ; % recons t ructed £ - domain
% s_recv= s_recv ( 1 : Lsig ) ; % reconstructed t - domain
% 
% % p l o t the i deal ly reconstructed s i gnal in t ime
% % and frequency domain
% subplot(211) ; sfig2a=plot(Faxis , abs(S_recv ) ) ;
% xlabel ( 'frequency ( Hz )' ) ;
% 
% title('Spectrum of ideal filtering (reconstruction)');
% subplot(212) ; sfig2b=plot(t,xsig,'r-.',t,s_recv( 1 : Lsig) ,'b');
% legend ( 'original signal' , 'reconstructed signal');
% xlabel ('time(sec)');
% title('original signal versus ideally reconstructed signal');
% set ( sfig2b , 'Linewidth',2);
% 
% % non - ideal reconstruc t i on
% ZOH=ones ( 1 , Nfactor ) ;
% s_ni =kron( downsample(s_out,Nfactor),ZOH);
% S_ni = fftshift(fft(s_ni,Lfft));
% S_ni=S_ni(:,1);
% S_recv2 = S_ni.*H_1pf ; % ideal f i l tering
% s_recv2 =real(ifft(fftshift(S_recv2))); % rec ons tructed £ - domai n
% s_recv2 = s_recv2(1:Lsig); % recons tructed t - doma in
% % p l o t the i deal ly rec ons tructed s i gnal in t ime
% % and frequency domain
% 
% figure(4)
% subplot(211);sfig3a=plot(t,xsig,'b',t,s_ni(1:Lsig),'r');
% xlabel('time(sec)');
% title('original signal versus flat-top reconstruction');
% subplot(212);sfig3b=plot(t,xsig,'b',t,s_recv2(1:Lsig),'r--');
% legend ( 'original signal','LPF reconstruction');
% xlabel ( 'time(sec)');
% set(sfig3a,'Linewidth',2); set(sfig3b,'Linewidth',2);
% title('original and flat-top reconstructi on a fter LPF');
% 



