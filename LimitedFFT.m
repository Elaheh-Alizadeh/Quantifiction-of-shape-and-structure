function [Parameters,FourierNames,ReconstructedBoundary]=LimitedFFT(Data,N)
Order=N+1;
n=length(Data);% number of the pixels in the boundary
%FS=1;
Y=fft(Data,n);% Run fourier transform
f=1*(0:n/2)/n;
t=(1:n);
N=Order;
An=real(2*((Y(2:N))/n));% Real part of fft is coefficient of cos function 
Bn=imag(2*((Y(2:N))/n));% Imaginary part of fft is coefficient of sin function 
A0=real((Y(1)/n));
%% Keep the first "Order" number of frequencies
Q=zeros(1,length(Y));
Q(1)=Y(1);
Q(2:Order)=2*Y(2:Order);
ReconstructedBoundary=real(ifft(Q));% reconstruct the boundary based on the first "Order" number of frequencies
%%
R2=1-sum((ReconstructedBoundary-Data).^2)/sum((Data-mean(Data)).^2);% Calculate the R squared value 
%%
Parameters=cat(1,1/n,R2,A0,(sqrt(An.^2+Bn.^2))');% We need rotation invariant parameter.
ExtraNames={'w','R-Squared','C0'};
FourierNames={};
for Num=1:Order-1; 
FourierNames=cat(2,FourierNames,{cat(2,'C',num2str(Num))});% this works;
end
FourierNames=cat(2,ExtraNames,FourierNames);
end