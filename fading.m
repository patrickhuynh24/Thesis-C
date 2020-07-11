clc;
clear all;
close all;
%% Declaring Variables
num_bit =10^6;
SNRdB = 0:1:20 ;                                                               %testing for 10 different SNR's
SNR=10.^(SNRdB/10);
%Fc=867*10^6;
%Pt=13; %13dB
var= 1; %2sigma^2
data = randi([0,1],[1,num_bit]);
L=10;
%% Adjustable Rician a parameters
kdBct= 0;
kdBtr= 0;
kct=10^(kdBct/10);
ktr=10^(kdBtr/10);
%% attenuation paramaters
act= abs(sqrt(kct/(kct+1))+sqrt((var/(kct+1))/2)*(randn(1,1))+1i*randn(1,1));
atr= abs(sqrt(ktr/(ktr+1))+sqrt((var/(ktr+1))/2)*(randn(1,1)+1i*randn(1,1)));
% phase_ct= 2*pi*Fc*24/(3*10^8);
% phase_tr =2*pi*Fc*24/(3*10^8);
phase_ct= unifrnd(0,2*pi);
phase_tr= unifrnd(0,2*pi);
a=act.*atr;
s=10^-9;
T=10; %1kbps
%% Varying carrier power to match SNR
Pc= SNR.*pi^2./(32*T*s^2);
m1=4*(sqrt(2*Pc)*a.*s)./pi;
delta_phi=unifrnd(0,2*pi);
phase_1=phase_ct+phase_tr+delta_phi; 
phi_1=unifrnd(0,2*pi);
phi_0=unifrnd(0,2*pi);
h= m1/2.*L*exp(-1i*phase_1); %compound channel hyperparameter
%e=[exp(1i*phi_0); exp(-1i*phi_0); exp(1i*phi_1); exp(-1i*phi_1)];

%% Noise generated
noise1= (sqrt(L/2))*(randn(1,num_bit) + 1i*randn(1,num_bit));
noise2= (sqrt(L/2))*(randn(1,num_bit) + 1i*randn(1,num_bit));
noise3= (sqrt(L/2))*(randn(1,num_bit) + 1i*randn(1,num_bit));
noise4= (sqrt(L/2))*(randn(1,num_bit) + 1i*randn(1,num_bit));

for k = 1:length(SNRdB) 
    Msend0pos(k)=h(k)*exp(1i*phi_0);
    Msend0neg(k)=h(k)*exp(-1i*phi_0);
    Msend1pos(k)=h(k)*exp(1i*phi_1);
    Msend1neg(k)=h(k)*exp(-1i*phi_1);
    %distribution of sending a 1
    oner0pos(k,:)=noise1;
    oner0neg(k,:)=noise2;
    oner1pos(k,:)=Msend1pos(k)+noise3;
    oner1neg(k,:)=Msend1neg(k)+noise4;
    %distribution of sending a 0;
    zeror0pos(k,:)=Msend0pos(k)+noise1;
    zeror0neg(k,:)=Msend1neg(k)+noise2;
    zeror1pos(k,:)=noise3;
    zeror1neg(k,:)=noise4;
end

%% Composite Hypothesis Testing With Rician Fading (non-coherent)
for i=1:length(SNRdB)
    c=0;
    err=0;
    for n=1:num_bit
       if data(n) == 1
            z1=0;
            z0=0;
            c=c+1;
            z1=(abs(oner1pos(i,n)))^2+(abs(oner1neg(i,n)))^2;
            z0=(abs(oner0pos(i,n)))^2+(abs(oner0neg(i,n)))^2;
        	if z0>z1
                err= err+1;
            end
       else 
            z1=0;
            z0=0;
            c=c+1;
            z1=(abs(zeror1pos(i,n)))^2+(abs(zeror1neg(i,n)))^2;
            z0=(abs(zeror0pos(i,n)))^2+(abs(zeror0neg(i,n)))^2;
        	if z0<z1
                err= err+1;
            end
       end
    end
    errors(i)=err/c;
end
BER= errors;
figure
semilogy(SNRdB,BER,'b-','LineWidth',2);
%axis([0 10 10^-6 1])
% hold on 
% semilogy(SNRdB,BER,'rx','LineWidth',2);
grid on
legend1=legend('Simulation');  
xlabel('SNR, dB');
ylabel('Bit Error Rate');
title('BER for FSK modulation bistatic backscatter modulation');

