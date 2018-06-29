function [Z A Phi] = Zernikmoment(p,n,m)
% -------------------------------------------------------------------------
% Copyright C 2012 Amir Tahmasbi
% The University of Texas at Dallas
% a.tahmasbi@utdallas.edu
% http://www.utdallas.edu/~a.tahmasbi/research.html
%
% License Agreement: You are free to use this code in your scientific
%                    research but you should cite the following papers:
%
% [1] - A. Tahmasbi, F. Saki, S. B. Shokouhi, 
%      	Classification of Benign and Malignant Masses Based on Zernike Moments, 
% 	J. Computers in Biology and Medicine, vol. 41, no. 8, pp. 726-735, 2011.
%
% [2] - A. Tahmasbi, F. Saki, H. Aghapanah, S. B. Shokouhi,
%	A Novel Breast Mass Diagnosis System based on Zernike Moments as Shape and Density Descriptors,
%	in Proc. IEEE, 18th Iranian Conf. on Biomedical Engineering (ICBME'2011), 
%	Tehran, Iran, 2011, pp. 100-104.
%
% -------------------------------------------------------------------------
% Function to find the Zernike moments for an N x N binary ROI
%
% [Z A Phi] = Zernikmoment(p,n,m)
% where
%   p = input image N x N matrix (N should be an even number)
%   n = The order of Zernike moment (scalar)
%   m = The repetition number of Zernike moment (scalar)
% and
%   Z = Complex Zernike moment 
%   A = Amplitude of the moment
%   Phi = phase (angle) of the mement (in degrees)
%
% Example: 
%   1- calculate the Zernike moment (n,m) for an oval shape,
%   2- rotate the oval shape around its centeroid,
%   3- calculate the Zernike moment (n,m) again,
%   4- the amplitude of the moment (A) should be the same for both images
%   5- the phase (Phi) should be equal to the angle of rotation

N = size(p,1);
x = 1:N; y = x;% number of columns and rows of a pixelized image
[X,Y] = meshgrid(x,y);
R = sqrt((2.*X-N-1).^2+(2.*Y-N-1).^2)/N;% radius (unit) of maped image of cell to unit circle
Theta = atan2((N-1-2.*Y+2),(2.*X-N+1-2)); % angle  of maped image of cell to unit circle
R = (R<=1).*R;%????????
Rad = radialpoly(R,n,m);    % get the radial polynomial

%A=exp(-1i*m*Theta);
Product = p(x,y).*Rad.*exp(-1i*m*Theta);%p(x,y)=f(c,r) in paper is function of image of the cell
Z = sum(Product(:));        % calculate the moments



cnt = nnz(R)+1;             % counts the number of (nonzero)pixels inside the unit circle. It is actualy lambda. What is one for?  
Z = (n+1)*Z/cnt;            % normalize the amplitude of moments
A =abs(Z);                 % calculate the amplitude of the moment
%A =Z;
Phi = angle(Z)*180/pi;      % calculate the phase of the mement (in degrees)


