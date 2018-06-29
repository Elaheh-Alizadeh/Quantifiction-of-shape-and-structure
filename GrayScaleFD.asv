function [FD,FDNames]=GrayScaleFD(Cell_Mask,GrayCell_Mask)
FD=zeros(1,4);
% Gray scale images are binerized using four different edge detection methods. 
Methods=[{'Roberts'},{'sobel'},{'zerocross'},{'canny'}];
[nr,nc]=size(GrayCell_Mask);
EdgedImages=zeros(nr,nc,length(Methods));
FDNames=strcat({'FD_'},Methods);

    for j=1:length(Methods)
        EdgedImages(:,:,j)=edge(GrayCell_Mask,Methods{j});
        FD(1,j)=hausDim(EdgedImages(:,:,j));
    end
end
