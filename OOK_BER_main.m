clc;
clear all;
close all;
%% Initialising
num_bit = 10^7;
SNRdB = -3:1:10 ;                                                          %testing for 10 different SNR's
SNR=10.^(SNRdB/10);
L = 10;
d = randi([0,1],[1,num_bit]);  

noise = (randn(L,num_bit) + 1i*randn(L,num_bit));                          %sigmaSquared/variance of noise = 1
a0=0.1;                                                                     
a1=9*a0;

%% Declaring a and b values of noiseless signal in the channel, then calculating non-centrality parameters
for i=1:length(SNRdB)
    gain(i)=  sqrt(2*SNR(i)/(a1-a0)^2);
    a(i)=a0*gain(i);                                                       %a value  times channel coefficient 
    sa2(i)=L*a(i)^2; 
    b(i)=a1*gain(i);                                                       %sb value 
    sb2(i)= L*b(i)^2; 
    n(i)= ((2/(sqrt(sb2(i))-sqrt(sa2(i))))*(L-0.5)*log(sqrt(sb2(i))/sqrt(sa2(i)))+(sqrt(sb2(i))+sqrt(sa2(i)))/2)^2;  

end

%% Receiving information and decoding
e = 0;
for m = 1:length(SNRdB)
    error=0;
    data =0;
    vcheck = zeros(1,num_bit);
    for k = 1:num_bit                                                      %checking each Bit for each SNR
        if d(k) == 1
            y(:) = 0;
            squares(:) =0;
            sums(:)=0;
            for ii = 1:L
                y(ii) = b(m)+noise(ii,k);
            end
            %Receiver
            squares= (abs(y)).^2;
            sums = sum(squares);
            if sums < n
              error = error +1;
              vcheck(k) = 0;                                               %noise has caused error
            else
              vcheck(k) = 1; 
            end
        elseif d(k) == 0
                                                                           %channel randomises 
            y(:) = 0;
            squares(:) =0;
            sums(:)=0;
            for iii =1:L
                y(iii) = a(m)+noise(iii,k);

            end
            %Receiver 
            squares= (abs(y)).^2;
            sums = sum(squares);
            if sums >= n
                error = error+1;
                vcheck(k) = 1;
            else
                vcheck(k) = 0;
            end
        end
    end
    e(m) = error;
end
%% Simulated BER
BER = e/num_bit;

%% Theoretical BER
theoryBER= 0.5-(0.5*marcumq((sqrt(sb2))/sqrt(0.5),(sqrt(n))/sqrt(0.5),L))+0.5*(marcumq((sqrt(sa2))/sqrt(0.5),(sqrt(n))/sqrt(0.5),L));

%%Plotting
close all
figure
semilogy(SNRdB,theoryBER,'b-','LineWidth',2);
hold on 
semilogy(SNRdB,BER,'rx','LineWidth',2);
axis([-3 8 10^-6 1])
grid on
legend('AWGN-Theory','AWGN-Simulation');
xlabel('SNR, dB');
ylabel('Bit Error Rate');
title('BER for OOK modulation in AWGN channel');