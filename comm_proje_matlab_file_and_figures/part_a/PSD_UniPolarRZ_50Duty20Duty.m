clear all;
close all;
clc;
%Nb is the number of bits to be transmitted
Nb=10000;
% Generate Nb bits randomly
b=rand(1,Nb)>0.5;
%Rb is the bit rate in bits/second
Rb=1000;
fs=10*Rb;
RZ_out=[];
Vp=5;
 
%% UniPolar %50 1000
Rb=1000;
fs=10*Rb;
RZ_out=[];
Vp=5;
for index=1:size(b,2)
 if b(index)==1
 RZ_out=[RZ_out [1 1 1 1 1 0 0 0 0 0]*Vp];
 elseif b(index)==0
 RZ_out=[RZ_out [1 1 1 1 1 0 0 0 0 0]*0];
 end
end
[Pxx f]=pwelch(RZ_out)
subplot(221)
plot(f/pi*fs/2,db(abs(Pxx))),grid
title('PSD - RZ - UniPolar - %50Duty - Rb=1kbps');
xlabel('Frequency (Hz)')
%% UniPolar %50 3000
Rb=3000;
fs=10*Rb;
RZ_out=[];
Vp=5;
for index=1:size(b,2)
 if b(index)==1
 RZ_out=[RZ_out [1 1 1 1 1 0 0 0 0 0]*Vp];
 elseif b(index)==0
 RZ_out=[RZ_out [1 1 1 1 1 0 0 0 0 0]*0];
 end
end
[Pxx f]=pwelch(RZ_out)
subplot(222)
plot(f/pi*fs/2,db(abs(Pxx))),grid
title('PSD - RZ - UniPolar - %50Duty - Rb=3kbps');
xlabel('Frequency (Hz)')
%% devam
Nb=10000;
% Generate Nb bits randomly
b=rand(1,Nb)>0.5;
%Rb is the bit rate in bits/second
Rb=1000;
fs=10*Rb;
RZ_out=[];
Vp=5;
 
%% UniPolar %20 1000
Rb=1000;
fs=10*Rb;
RZ_out=[];
Vp=5;
for index=1:size(b,2)
 if b(index)==1
 RZ_out=[RZ_out [1 1 0 0 0 0 0 0 0 0]*Vp];
 elseif b(index)==0
 RZ_out=[RZ_out [1 1 0 0 0 0 0 0 0 0]*0];
 end
end
[Pxx f]=pwelch(RZ_out)
subplot(223)
plot(f/pi*fs/2,db(abs(Pxx))),grid
title('PSD - RZ - UniPolar - %20Duty - Rb=1kbps');
xlabel('Frequency (Hz)')
%% UniPolar %20 3000
Rb=3000;
fs=10*Rb;
RZ_out=[];
Vp=5;
for index=1:size(b,2)
 if b(index)==1
 RZ_out=[RZ_out [1 1 0 0 0 0 0 0 0 0]*Vp];
 elseif b(index)==0
 RZ_out=[RZ_out [1 1 0 0 0 0 0 0 0 0]*0];
 end
end
[Pxx f]=pwelch(RZ_out)
subplot(224)
plot(f/pi*fs/2,db(abs(Pxx))),grid
title('PSD - RZ - UniPolar - %20Duty - Rb=3kbps');
xlabel('Frequency (Hz)')
 


