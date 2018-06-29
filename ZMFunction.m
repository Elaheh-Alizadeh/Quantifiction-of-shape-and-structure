function [ZernikeMoment,ZernikeNames]=ZMFunction(Cell_Mask,nZernikeStart,nZernikeStop)
%
Counter=0;

%%
%Calculate ZM
for n=nZernikeStart:nZernikeStop;
    if rem(n,2)==0
        M=0:2:n;
    elseif rem(n,2)==1
        M=1:2:n;  
    end
    for k=1:length(M);
        Counter=Counter+1; 
        %p = logical(not(Cell_Mask));% for black as cell image
        [~,ZernikeMoment(1,Counter),~ ] = Zernikmoment(Cell_Mask,n,M(k));
        ZernikeNames(1,Counter)={[ num2str(n) '_' num2str(M(k)) ]};    
    end
end
     

        
        
        
        
        
        
        
        
        
      