% Author: Patrick Huu Danh Huynh
% Title: BER for FSK modulation in elliptical topology
clc;
clear all;
close all;
%% Initialising assumed values 
num_bit =10^6;
%SNRdB = 0:1:10 ;                                                               %testing for 10 different SNR's
%SNR=10.^(SNRdB/10);

F1= 900*10^6;           
Fc=867*10^6;
T=1/1000; %bit period for bitrate = 1kbps
Fs= 10000;
Ts=1/Fs;
L = T/Ts;

d = randi([0,1],[1,num_bit]);
%% Properties of a ellipse 
eccentricity = 12/13;                                                       %for an ellipse with (x/26)^2+(y/10)^2=1
distance_ER= 48;                                                            %distance between emitter and reader
distance = 2:1:50;                                                          %distance from emitter to tag
y= 48/(eccentricity)-distance;                                              %distance from tag to reader
%% 
min=0;
maximum=2*pi;
startphase= min+rand*(maximum-min);                                         %random start phase

t_tr=y/(3*10^8);

ctr=2*pi*Fc*t_tr*3;


for p=1:length(y)
    mhat(p)=2*10^5/pi.*sqrt(10^-4./(distance(p)*y(p)).^2); 
end

s2= ((mhat.^2)*(L^2))/2;                                                    %non-centrality parameter

M1= mhat.*L/2.*exp(1i*((startphase-2*pi*F1*t_tr)-(-ctr)));

%% BER generation
sig= L/2;
for z= 1:length(distance)
    fun=@(x) (1- igamma(2,x/(2*sig))).*(1/2).*1/sig.*exp((-((x./sig)+(s2(z)./sig))/2)).*((x./s2(z)).^(1/2)).*besseli(1,sqrt(x*s2(z))./sig);
    success(z)= integral(fun,0,10^5);                                      
    theoryBER(z)= 1- success(z);
end

figure
semilogy(distance,theoryBER,'b-','LineWidth',2);
hold on 
grid on
legend1=legend('Theory');  
xlabel('Distance from emitter, metres');
ylabel('Bit Error Rate');
title('BER for FSK modulation bistatic backscatter modulation');