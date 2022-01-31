clear all;
close all;
clc;
%Nb is the number of bits to be transmitted
%% Polar
Nb=10000;
% Generate Nb bits randomly
b=rand(1,Nb)>0.5;
%Rb is the bit rate in bits/second
Rb=1000;
fs=10*Rb;
NRZ_out=[];
Vp=5;
 
%Line Coding
for index=1:size(b,2)
 if b(index)==1
 NRZ_out=[NRZ_out [1 1 1 1 1 1 1 1 1 1]*Vp];
 elseif b(index)==0
 NRZ_out=[NRZ_out [1 1 1 1 1 1 1 1 1 1]*(-Vp)];
 end
end
[Pxx f]=pwelch(NRZ_out)
subplot(211)
plot(f/pi*fs/2,db(abs(Pxx))),grid
title('PSD - NRZ - Polar - Rb=1kbps');
xlabel('Frequency (Hz)')
 
 
Nb=10000;
% Generate Nb bits randomly
b=rand(1,Nb)>0.5;
%Rb is the bit rate in bits/second
Rb=3000;
fs=10*Rb;
NRZ_out=[];
Vp=5;
 
%Line Coding
for index=1:size(b,2)
 if b(index)==1
 NRZ_out=[NRZ_out [1 1 1 1 1 1 1 1 1 1]*Vp];
 elseif b(index)==0
 NRZ_out=[NRZ_out [1 1 1 1 1 1 1 1 1 1]*(-Vp)];
 end
end
[Pxx f]=pwelch(NRZ_out)
subplot(212)
plot(f/pi*fs/2,db(abs(Pxx))),grid
title('PSD - NRZ - Polar - Rb=3kbps');
xlabel('Frequency (Hz)')
 
 
 
 


