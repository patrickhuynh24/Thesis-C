% Author: Patrick Huu Danh Huynh
% Title: BER for FSK modulation in parallel topology
clc;
clc;
clear all;
close all;
%% Initialising known and random values
num_bit =10^6;
%SNRdB = 0:1:10 ;                                                           %testing for 10 different SNR's
%SNR=10.^(SNRdB/10);

F1= 900*10^6;           
Fc=867*10^6;
T=1/1000;                                                                   %bit period for bitrate = 1kbps
Fs= 10000;
Ts=1/Fs;
L = T/Ts;
h= 0.2;
d = randi([0,1],[1,num_bit]);
x = 1:1:48;
%% Distance forumla / trigonometry
distance= sqrt(x.^2+h^2);
limit= max(x);

min=0;
maximum=2*pi;
startphase= min+rand*(maximum-min);                                         %random start phase
t_tr=distance/(3*10^8);

ctr=2*pi*Fc*t_tr*3;
for m= 1:limit
    SNR(m)= 5*(10^9)/pi*(10^-2*10^-2)/((distance(m)^2)*(sqrt((limit-x(m))^2+h^2))^2);
end

for p=1:length(SNR)
    mhat(p)=sqrt(SNR(p)*2/L); 
end

s2= ((mhat.^2)*(L^2))/2;  %non-centrality parameter

%M0= mhat.*L/2.*exp(-1i*(startphase-2*pi*F0*t_tr)-(-ctr));
M1= mhat.*L/2.*exp(1i*((startphase-2*pi*F1*t_tr)-(-ctr)));

%% BER calculation
sig= L/2;
for z= 1:length(SNR)
    fun=@(x) (1- igamma(2,x/(2*sig))).*(1/2).*1/5.*exp((-((x./sig)+(s2(z)./sig))/2)).*((x./s2(z)).^(1/2)).*besseli(1,sqrt((x/sig).*(s2(z)./sig)));
    success(z)= integral(fun,0,100000);                                      
    theoryBER(z)= 1- success(z);
end

%% Plotting
figure
semilogy(distance,theoryBER,'r-','LineWidth',2);
axis([0 48 10^-2 1])
hold on 
semilogy(distance,BER,'rx','LineWidth',2);
grid on
legend1=legend('FSK');  
xlabel('Distance from emitter, metres');
ylabel('Bit Error Rate');
title('BER for FSK modulation bistatic backscatter modulation');

    