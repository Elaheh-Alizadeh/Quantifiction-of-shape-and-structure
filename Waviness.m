function [FourierPara,FourierNames,NonPeriodicBoundary]=Waviness(CellBoundary,L,CartesianOrPolar)
 
 
 
if strcmp(CartesianOrPolar,'Polar')
     Tag1={'_R'};
     Tag2={'_Theta'};
     
    
     %Since theta increases regardless of the oscillations, we need to subtract
     %the phase jump with intercept equal to (2 pi/number of pixels in the boundary).
     nI=size(CellBoundary,2);
     m=2*pi/nI;
     Intercept=2*pi;
     t=(1:nI)';
     Yfixed=-m*t+Intercept;
     CorrectedPhaseJump=unwrap(CellBoundary(2,1:end))';
     Remaining=CorrectedPhaseJump-Yfixed;
     CellBoundary=cat(1,CellBoundary(1,1:end),Remaining');
     disp('Polar Zone')
 
elseif strcmp(CartesianOrPolar,'Cartesian');  
    Tag1={'_X'};
    Tag2={'_Y'};
    disp('Cartesian Zone')
end
 
%Calculate Fourier coefficients for the first component in polar or Cartesian coordinates(X or R)
[A1,Names,yfit1]=LimitedFFT(CellBoundary(1,1:end),L);
 
%Calculate Fourier coefficients for the first component in polar or Cartesian coordinates(Y or Theta)
[A2,~,yfit2]=LimitedFFT(CellBoundary(2,1:end),L);
 
ReconstructedBoundary=cat(1,yfit1,yfit2);
FourierPara=cat(2,A1',A2');
FourierNames=cat(2,strcat(Names,Tag1),strcat(Names,Tag2));
   
    
NonPeriodicBoundary=CellBoundary-ReconstructedBoundary;
 

end

