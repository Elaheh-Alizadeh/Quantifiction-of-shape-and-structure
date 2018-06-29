function [Roughness_Para,RoughnessNames]=Roughness(y,CartesianOrPolar)

if strcmp(CartesianOrPolar,'Polar')
     Tag1={'_R'};
     Tag2={'_Theta'};

     disp('Polar Zone')
elseif strcmp(CartesianOrPolar,'Cartesian');  
    Tag1={'_X'};
    Tag2={'_Y'};
    disp('Cartesian Zone')
end

[n,nIndeces]=size(y);

if n>2
    error('Input has more than Theta and radius information')
end




Names=[{'MeanAbsolute'}, {'RootMeanSquared'},...
    {'TenPointMeanRoughness'},{'MaximumPeak_ValleyRoughness'},{'MaximumHeight'},{'Skewness'},{'Kurtosis'}];
   
RoughnessNames=cat(2,strcat(Names,Tag1),strcat(Names,Tag2))  ;   

%1) Arithmetic average of absolute values
R_a=mean(abs(y),2);

%2) Root Mean Squared
R_q=rms(y,2);%[1;1];

%3) Ten Point mean Roughness
for i=1:n
    Rp(i,:)=max(y(i,:));
    Rv(i,:)=min(y(i,:));
    Y_Sort(i,:)=sort(y(i,:),2);
    Y_V(i,:)=Y_Sort(i,1:5);
    Y_P(i,:)=Y_Sort(i,length(Y_Sort(i,:))-5:end);
    R_ZIJS(i,:)=mean(Y_V(i,:))+mean(Y_P(i,:));
    %7) Skewness
    R_sk(i,:)=sum(y(i,:).^3)/nIndeces/R_q(i,:)^3;

    %8) Kurtosis
    R_ku(i,:)=sum(y(i,:).^4)/nIndeces/R_q(i,:)^4;

end

%5) Maximum peak to valley roughness
R_y=Rp+Rv;

%6) Maximum Height of the Profile
R_t=Rp-Rv;


Parameters=cat(2,R_a,R_q,R_ZIJS,R_y,R_t,R_sk,R_ku);

Roughness_Para=cat(2,Parameters(1,:),Parameters(2,:));
end 


