% Author: Patrick Huu Danh Huynh
% Title: OOK modulation for bistatic backscatter communications with BER in
% terms of distance

clc;
clear all;
close all;
%% Initialising
distance = 1:1:48;                                                         %distance from emitter
L = 10;                                                                    %L taps

a0=0.1;                                                                    %a                                                                
a1=9*a0;                                                                   %b
h=1;
%% Declaring a and b values of noiseless signal in the channel, then calculating non-centrality parameters
for i=1:length(distance)
    coeff(i)=  sqrt((0.01*0.01)/((sqrt(distance(i).^2+h^2))^2*(sqrt(48-distance(i))^2+h^2)));
    
    a(i)=a0*coeff(i);                                                      %a value  times channel coefficient 
    sa2(i)=L*(a(i))^2; 
    sa(i) = sqrt(sa2(i));
    
    b(i)=a1*coeff(i);                                                      %sb value 
    sb2(i)= L*(b(i))^2;
    sb(i) = sqrt(sb2(i));
    
    n(i)= ((0.5/(sb(i)-sa(i)))*(L-0.5)*log(sb(i)/sa(i))+(sb(i)+sa(i))/2)^2;  

end
sig = 1/sqrt(2);
%% Theoretical BER
theoryBER= 0.5-0.5*marcumq(sb/sig,sqrt(n)/sig,L)+0.5*(marcumq(sa/sig,sqrt(n)/sig,L));

%% Plotting
semilogy(distance,theoryBER,'bx-','LineWidth',2);
grid on
legend2=legend('FSK','OOK');
xlabel('distance, m');
ylabel('Bit Error Rate');
