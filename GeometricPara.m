%% Geometric Parameter Measurements of Cell Shape (36 Attributes)
%  Note: 
%  modified from Binary_func.m in 
%  Version1Jun2015-29Attributes-CellExtractionPlusGeometryParametes 
%  March 17, 2016 
%  Wenlong Xu 
%  Prasad Group 
%  Colorado State Univ. 
%  Updated on March 28, 2016 
%  ------------------------------------------------------------------------
%% Introduction 
%  The original code includes both image reconstruction from nuclei and the
%  cell shape parameter measurements. The first step was replaced by
%  CellProfiler and the following image processing pipeline. To work better
%  with the new pipeline, the second part was modified in this code. 
%  ------------------------------------------------------------------------
%% Results are Same as Josh's original code
%% New additions to the original 29 attributes 
%  1) Distance between the centroids of cell and nuclei; 
%  2) The angle between the 1) and the orientation of the cell; 
%  3) The ratio between perimeters of convex hull and the cell; 
%  4) Measurements between the points on cell boundary and the cell
%  entroid, including the max, the min, the max/min, the mean and the std. 
%% Inputs 
% WD = 'C:\Users\xuwl\Desktop\CellShapeProject\rawData\DLM8_GDA_TIFF\SepCells\DLM8_GDA_1_Cell11\'; 
% BWCell_File = 'DLM8_GDA_1_Cell11_CellMask.bmp'; 
% BWNuc_File  = 'DLM8_GDA_1_Cell11_NucMask.bmp'; 
% Cell_Mask = imread([WD, BWCell_File]); 
% Nuc_Mask  = imread([WD, BWNuc_File]); 
function [Measurements36,CellBoundryRT,Cell_Radius_xy_Vals,GeometricParaNames]= GeometricPara(Cell_Mask, Nuc_Mask,Nuc_Mask2)%CellBoundry%Cell_Radius_xy_Vals

GeometricParaNames=[{'Cell_Area'} ,{'Cell_Perimeter'} ,{'Cell_Major'} ,{'Cell_Minor'} ,{'Cell_Circularity'} ...
               {'Cell_AspectRatio'} ,{'Cell_Roundness'} ,{'Cell_Solidity'} ,{'Max_Cell_Radius'} ,{'Ratio_Cell_Radius'}...
              { 'Mean_Cell_Radius'} ,{'CV_Cell_Radius'} ,{'Max_Span'} ,{'Hull_Area'} ,{'Hull_Perimeter'}  ...
               {'Ratio_Hull2Cell_Perimeter'} ,{'Hull_Circularity'} ,{'Max_Hull_Radius'} ,{'Ratio_Hull_Radius'} ,{'Mean_Hull_Radius'} ...
               {'CV_Hull_Radius'} ,{'Circle_Diameter'} ,{'Max_Circle2Hull_Radius'} ,{'Ratio_Circle2Hull_Radius'} ,{'Mean_Circle2Hull_Radius'}...
              { 'CV_Circle2Hull_Radius'} ,{'Nuc_Area'} ,{'Nuc_Perimeter'} ,{'Nuc_Major'} ,{'Nuc_Minor'}...
               {'Nuc_Circularity'} ,{'Nuc_AspectRatio'} ,{'Nuc_Roundness'} ,{'Nuc_Solidity'} ,{'Distance_Nuc2Cell_Centroid'}... 
               {'Nuc2Cell_Orientation'}];
Cell_Props = regionprops(Cell_Mask, 'Area', 'Perimeter', 'MajorAxisLength', 'MinorAxisLength', ...
                                    'Solidity', 'ConvexHull', 'ConvexImage', 'BoundingBox', 'PixelList', ...
                                    'Centroid', 'Orientation'); 
Nuc_Props  = regionprops(Nuc_Mask,  'Area', 'Perimeter', 'MajorAxisLength', 'MinorAxisLength', 'Solidity', ...
                                    'Centroid'); 
if length(Cell_Props) > 1 || length(Nuc_Props) > 1 
    warn('Multiple ROIs identified in the input image!\n'); 
end 
%  ------------------------------------------------------------------------
%% 12 Geometric Parameters for the cell from basic measurements in regionprops 
% 1. Cell Area (pixel^2): Area of the cell in pixels
Cell_Area = Cell_Props(1).Area; 
% 2. Cell Perimeter (pixels): Perimeter in pixels of the cell 
Cell_Perimeter = Cell_Props(1).Perimeter;
% 3. Cell Major (pixels): The major axis of an ellipse drawn around the cell 
Cell_Major = Cell_Props(1).MajorAxisLength;
% 4. Cell Minor (pixels): The minor axis of an ellipse drawn around the cell 
Cell_Minor = Cell_Props(1).MinorAxisLength;
% 5. Circularity of Cell (unitless (0-1)): 4pi(Area)/Perimeter^2 
Cell_Circularity = 4*pi*Cell_Area/Cell_Perimeter^2; 
% 6. Aspect Ratio of Cell (unitless): Major/minor axis of ellipse 
Cell_AspectRatio = Cell_Major/Cell_Minor; 
% 7. Roundness of Cell (unitless (0-1)): 4(Area)/(pi(Elliptical Major Axis)^2) 
Cell_Roundness = 4*Cell_Area/(pi*Cell_Major^2); 
% 8. Solidity of Cell (unitless (0-1)): 
% Area of convex hull of cell / Area of the cell 
Cell_Solidity = Cell_Props(1).Solidity; 
% The following are NEWLY added measurements of the cell radius, 
% defined as the distance between the cell centroid and points on the cell
% boundaries 
Cell_Border = bwboundaries(Cell_Mask); 
Cell_Border = Cell_Border{1}; 
Cell_Centroid = Cell_Props(1).Centroid; 
Num_Points_Cell_Border = size(Cell_Border, 1); 
Cell_Radius = sqrt(sum((Cell_Border - repmat(Cell_Centroid, Num_Points_Cell_Border, 1)).^2, 2)); 
Cell_Radius_xy=Cell_Border - repmat(Cell_Centroid, Num_Points_Cell_Border, 1);

Cell_Radius_xy_Vals=Cell_Radius_xy';% To export like Theriots paper    
Theta=atan(Cell_Radius_xy(:,2)./Cell_Radius_xy(:,1));

t=1:length(Theta);
for i=1:length(Theta)
    if Theta(i,1)<0 
        if Cell_Radius_xy(i,1)>0 && Cell_Radius_xy(i,2) <0
            Theta(i,1)=Theta(i,1)+360;
        elseif Cell_Radius_xy(i,1)<0 && Cell_Radius_xy(i,2) >0
            Theta(i,1)=Theta(i,1)+180;
        end
    elseif Theta(i,1)>0 
        if Cell_Radius_xy(i,1)>0 && Cell_Radius_xy(i,2) >0
            Theta(i,1)=Theta(i,1);
        elseif Cell_Radius_xy(i,1)<0 && Cell_Radius_xy(i,2) <0
            Theta(i,1)=Theta(i,1)+180;
        end
    end
end

CellBoundryRT=cat(1,transpose(Cell_Radius),transpose(Theta));%In polar coordinates
% 9. Max radius of the cell (pixels) (new): 
% Maximal Distance from centroid of cell to an exterior point on the cell border  
Max_Cell_Radius = max(Cell_Radius);   
% 10. Ratio of Max/Min Radius of the cell (unitless) (new): 
% Maximum/minimum radius of the cell (see max radius definition above) 
Min_Cell_Radius = min(Cell_Radius); 
Ratio_Cell_Radius = Max_Cell_Radius/Min_Cell_Radius; 
% 11. Mean Radius of the cell (pixels) (new): 
% The mean of all radii drawn from the hull's centroid to an exterior point
Mean_Cell_Radius = mean(Cell_Radius); 
% 12. CV Rad Cell (unitless) (new): 
% The relative variation of radii drawn from the cell’s center to an exterior 
% point. Given by the Standard of Deviation of all Radii divided by the 
% mean of all radii 
STD_Cell_Radius = std(Cell_Radius); 
CV_Cell_Radius = STD_Cell_Radius/Mean_Cell_Radius; 
%  ------------------------------------------------------------------------
%%  9 Geometric Parameters from minimal bounding hull
% Pixels consist of the minimal bounding convex hull: 
Hull_Pixels = Cell_Props(1).ConvexHull; 
% Calculate the distances of all possible pairs of points on the convex
% hull: 
All_Span = pdist(Hull_Pixels); 
% 13. Max span across the hull (pixels): 
% The maximum distance from one point on the convex hull to another 
Max_Span = max(All_Span); 
% The binary image that specifies the convex hull 
Hull_Image = Cell_Props(1).ConvexImage; 
%figure,imshow(Hull_Image)
Hull_Props = regionprops(Hull_Image, 'Area', 'Perimeter', 'Centroid'); 
% 14. Area of the hull (pixel^2): 
% Area of the convex hull (smallest convex polygon that encloses the cell) 
Hull_Area = Hull_Props(1).Area; 
% 15. Perimeter of the hull (pixels): 
% Perimeter of the convex hull 
Hull_Perimeter = Hull_Props(1).Perimeter; 
% 16. Ratio of hull perimeter to cell perimeter (new): 
Ratio_Hull2Cell_Perimeter = Hull_Perimeter/Cell_Perimeter; 
% 17. Circularity of the hull (unitless (0-1)): 4pi(Area)/Perimeter^2  
Hull_Circularity = 4*pi*Hull_Area/Hull_Perimeter^2; 
% All distances from the centroid of Hull to an exterior point on the hull 
Hull_Centroid = Hull_Props.Centroid; 
% This centroid is relative to the bounding box, so we need to add the
% bounding box: 
Cell_BoundingBox = Cell_Props(1).BoundingBox; 
Hull_Centroid = Hull_Centroid + floor(Cell_BoundingBox(1:2)); 
Num_Points_on_Hull = size(Hull_Pixels, 1); 
Hull_Radii = zeros(1, Num_Points_on_Hull); 
for ii = 1:Num_Points_on_Hull 
    Hull_Radii(ii) = sqrt((Hull_Centroid(1)-Hull_Pixels(ii, 1))^2 + (Hull_Centroid(2)-Hull_Pixels(ii, 2))^2); 
end
% 18. Max radius of the hull (pixels): 
% Maximum Distance from centroid of Hull to an exterior point on the hull 
Max_Hull_Radius = max(Hull_Radii); 
% 19. Max/Min Radius of hull (unitless): 
% Maximum/minimum radius of hull (see max radius definition above) 
Min_Hull_Radius = min(Hull_Radii); 
Ratio_Hull_Radius = Max_Hull_Radius/Min_Hull_Radius; 
% 20. Mean Radius of hull (pixels): 
% The mean of all radii drawn from the hull's centroid to an exterior point
Mean_Hull_Radius = mean(Hull_Radii); 
% 21. CV Rad Hull (unitless): 
% The relative variation of radii drawn from the hull’s center to an exterior 
% point. Given by the Standard of Deviation of all Radii divided by the 
% mean of all radii 
STD_Hull_Radius = std(Hull_Radii); 
CV_Hull_Radius = STD_Hull_Radius/Mean_Hull_Radius; 
%  ------------------------------------------------------------------------
%% 5 Geometric Parameters from Minimal Bounding Circle 
Cell_PixelList = Cell_Props(1).PixelList; 
% minboundcircle is an open-source code from John D'Errico 
% The flag in the 3rd input is hullflag to not use its built-in
% functionality to calculate the convex hull. 
[Circle_Center, Circle_Radius] = minboundcircle(Cell_PixelList(:,1),Cell_PixelList(:,2),0); 
% 22. Diameter of bounding circle (pixels): 
% The diameter of the bounding circle drawn around the cell
Circle_Diameter = 2*Circle_Radius; 
% Calculate all radius from the center of the bounding circle to hull: 
Circle2Hull_Radii = zeros(1, Num_Points_on_Hull); 
for ii = 1:Num_Points_on_Hull 
    Circle2Hull_Radii(ii) = sqrt((Circle_Center(1)-Hull_Pixels(ii, 1))^2 + (Circle_Center(2)-Hull_Pixels(ii, 2))^2); 
end
% 23. Max radius from circle to hull (pixels): 
% The maximum distance from the center of the circle to an edge of the convex hull 
Max_Circle2Hull_Radius = max(Circle2Hull_Radii); 

% 24. Max/Min Radius from Circle (unitless): 
% Maximum/minimum radius from the center of the circle to an edge of the convex hull 
Min_Circle2Hull_Radius = min(Circle2Hull_Radii); 
Ratio_Circle2Hull_Radius = Max_Circle2Hull_Radius/Min_Circle2Hull_Radius; 
% 25. Mean radius from circle to hull (pixels): 
% The mean of all radii drawn from the circle's center to the hull 
Mean_Circle2Hull_Radius = mean(Circle2Hull_Radii); 
% 26. CV Radius from circle to hull (unitless): 
% The relative variation of radii drawn from the circle's center to the hull. 
% Given By the Standard of Deviation of all Radii divided by the mean of all radii
Std_Circle2Hull_Radius = std(Circle2Hull_Radii); 
CV_Circle2Hull_Radius = Std_Circle2Hull_Radius/Mean_Circle2Hull_Radius; 
%  ------------------------------------------------------------------------
%% 10 Geometric Parameters for the Nucleus 
% 27. Nucleus Area (pixel^2): Area of the Nucleus in pixels
Nuc_Area = Nuc_Props(1).Area; 
% 28. Nucleus Perimeter (pixels): Perimeter in pixels of the Nucleus 
Nuc_Perimeter = Nuc_Props(1).Perimeter;
% 29. Nucleus Major (pixels): The major axis of an ellipse drawn around the Nucleus 
Nuc_Major = Nuc_Props(1).MajorAxisLength;
% 30. Nucleus Minor (pixels): The major axis of an ellipse drawn around the Nucleus 
Nuc_Minor = Nuc_Props(1).MinorAxisLength;
% 31. Circularity of Nucleus (unitless (0-1)): 4pi(Area)/Perimeter^2 
Nuc_Circularity = 4*pi*Nuc_Area/Nuc_Perimeter^2; 
% 32. Aspect Ratio of Nucleus (unitless): Major/minor axis of ellipse 
Nuc_AspectRatio = Nuc_Major/Nuc_Minor; 
% 33. Roundness of Nucleus (unitless (0-1)): 4(Area)/(pi(Elliptical Major Axis)^2) 
Nuc_Roundness = 4*Nuc_Area/(pi*Nuc_Major^2); 
% 34. Solidity of Nucleus (unitless (0-1)): 
% Area of convex hull of Nucleus / Area of the Nucleus 
Nuc_Solidity = Nuc_Props(1).Solidity; 
% 35. Distance between the centroids of cell and nucleus (pixels) (new):

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Elaheh Made Changes for Nuc

Nuc_Props2= regionprops(Nuc_Mask2,'Centroid'); 
Nuc_Centroid2=Nuc_Props2(1).Centroid;
%Nuc_Centroid = Nuc_Props(1).Centroid; 
Distance_Nuc2Cell_Centroid = sqrt(sum((Cell_Centroid - Nuc_Centroid2).^2)); 
% 36. Angle between the 35) and cell major axis (degree) (new):

Orientation_Nuc2Cell_Centroid = atan((Nuc_Centroid2(2) - Cell_Centroid(2))/(Nuc_Centroid2(1) - Cell_Centroid(1))); 
Orientation_Nuc2Cell_Centroid = Orientation_Nuc2Cell_Centroid * 180/pi; 
Cell_Orientation = Cell_Props(1).Orientation; 
Nuc2Cell_Orientation = abs(Orientation_Nuc2Cell_Centroid - Cell_Orientation); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if Method==1
% out=CellBoundry;
% elseif Method==2
% out=Cell_Radius_xy_Vals;
% end
    
%  ------------------------------------------------------------------------
%% Creat the output vector 
Measurements36 = [Cell_Area, Cell_Perimeter, Cell_Major, Cell_Minor, Cell_Circularity, ...
                 Cell_AspectRatio, Cell_Roundness, Cell_Solidity, Max_Cell_Radius, Ratio_Cell_Radius, ...
                 Mean_Cell_Radius, CV_Cell_Radius, Max_Span, Hull_Area, Hull_Perimeter, ...
                 Ratio_Hull2Cell_Perimeter, Hull_Circularity, Max_Hull_Radius, Ratio_Hull_Radius, Mean_Hull_Radius, ...
                 CV_Hull_Radius, Circle_Diameter, Max_Circle2Hull_Radius, Ratio_Circle2Hull_Radius, Mean_Circle2Hull_Radius, ...
                 CV_Circle2Hull_Radius, Nuc_Area, Nuc_Perimeter, Nuc_Major, Nuc_Minor, ...
                 Nuc_Circularity, Nuc_AspectRatio, Nuc_Roundness, Nuc_Solidity, Distance_Nuc2Cell_Centroid, ...
                 Nuc2Cell_Orientation]; 
%  Names=(['Cell_Area\t Cell_Perimeter\t Cell_Major\t Cell_Minor\t Cell_Circularity\t ' ...
%                'Cell_AspectRatio\t Cell_Roundness\t Cell_Solidity\t Max_Cell_Radius\t Ratio_Cell_Radius\t '...
%                'Mean_Cell_Radius\t CV_Cell_Radius\t Max_Span\t Hull_Area\t Hull_Perimeter\t ' ...
%                'Ratio_Hull2Cell_Perimeter\t Hull_Circularity\t Max_Hull_Radius\t Ratio_Hull_Radius\t Mean_Hull_Radius\t' ...
%                'CV_Hull_Radius\t Circle_Diameter\t Max_Circle2Hull_Radius\t Ratio_Circle2Hull_Radius\t Mean_Circle2Hull_Radius\t '...
%                'CV_Circle2Hull_Radius\t Nuc_Area\t Nuc_Perimeter\t Nuc_Major\t Nuc_Minor\t '...
%                'Nuc_Circularity\t Nuc_AspectRatio\t Nuc_Roundness\t Nuc_Solidity\t Distance_Nuc2Cell_Centroid\t'... 
%                'Nuc2Cell_Orientation','Euler']);
             



end % end of function 