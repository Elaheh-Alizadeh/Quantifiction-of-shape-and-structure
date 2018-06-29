function [Features,FeatureNames]=CalculateParameters(Cell_Binary,Nuc_Binary,Nuc_Binary2,Cell_Gray,nZernikeStart,nZernikeStop,FSMoments)

%% Textural measures
[Features.BandBased,FeatureNames.BandBased]=Band_Based_Measurements(Cell_Binary, Cell_Gray);

[Features.FractalDimension,FeatureNames.FractalDimension]=GrayScaleFD(Cell_Binary,Cell_Gray);

[Features.GrayScale,FeatureNames.GrayScale]=GLCM_Features4(Cell_Gray);% From Gray Scale Cells     

%% Spreading measure
[Features.Geometric,CellBoundaryRT,CellBoundaryXY,FeatureNames.Geometric]=GeometricPara(Cell_Binary, Nuc_Binary,Nuc_Binary2);  

[Features.ZernikeMoment,FeatureNames.ZernikeMoment]=ZMFunction(Cell_Binary,nZernikeStart,nZernikeStop);

   
%% Irregularity measures
      
[Features.WavinessRT,FeatureNames.WavinessRT,NonPeriodicBoundaryRT]=Waviness(CellBoundaryRT,FSMoments,'Polar')  
[Features.WavinessXY,FeatureNames.WavinessXY,NonPeriodicBoundaryXY]=Waviness(CellBoundaryXY,FSMoments,'Cartesian')  
     
[Features.RoughnessRT,FeatureNames.RoughnessRT]=Roughness(NonPeriodicBoundaryRT,'Polar');
[Features.RoughnessXY,FeatureNames.RoughnessXY]=Roughness(NonPeriodicBoundaryXY,'Cartesian');
     
     
end