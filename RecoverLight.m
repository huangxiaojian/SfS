function [ ] = RecoverLight( )
%RECOVERLIGHT Summary of this function goes here
%   Detailed explanation goes here

%I = imread('RefIntensity.jpg');
%MaskImage = imread('MaskImage.jpg');

albedo = 0.8;
I = importdata('RefImage.mat'); %input data
MaskImage = importdata('RefMaskImage.mat'); %input data
Nref = importdata('RefNormalizedNormalsV.mat'); %input data
%Nref = RefNormalizedNormalsV;
[imageWidth, imageHeight] = size(I);

length = 4;

A = zeros(length, length);
LightDirV = zeros(length, 1);
b = zeros(length, 1);

for k=1:length
    for i = 1:imageWidth
       for j = 1:imageHeight
           if(MaskImage(i,j) > 0.5)
               Nx = Nref(((i-1)*imageHeight+j-1)*3+1);
               Ny = Nref(((i-1)*imageHeight+j-1)*3+2);
               Nz = Nref(((i-1)*imageHeight+j-1)*3+3);
               
               if length == 9
                   Y = [1, Nx, Ny, Nz, Nx*Ny, Nx*Nz, Ny*Nz, Nx*Nx-Ny*Ny, 3*Nz*Nz-1];
               else
                   Y = [1, Nx, Ny, Nz];
               end
               
               b(k) = b(k) + I(i,j)*Y(k);
               %A(k,:) = A(k,:)+ albedo*double(Y(k))*double(Y);
               for ii = 1:length
                   A(k, ii) = A(k, ii) + albedo*Y(k)*Y(ii); 
               end
               %A(k,:) = A(k,:)+ albedo*Y(k)*Y;
           end
       end
    end
end

LightDirV = A\b;

save LightDirV.mat LightDirV;


%Use the LightDirV to calculate the image intensity
ComputedImage = zeros(imageWidth, imageHeight);
for k=1:length
    for i = 1:imageWidth
       for j = 1:imageHeight
           if(MaskImage(i,j) > 0.5)
               Nx = Nref(((i-1)*imageHeight+j-1)*3+1);
               Ny = Nref(((i-1)*imageHeight+j-1)*3+2);
               Nz = Nref(((i-1)*imageHeight+j-1)*3+3);
               %Y = [1, Nx, Ny, Nz, Nx*Ny, Nx*Nz, Ny*Nz, Nx*Nx-Ny*Ny, 3*Nz*Nz-1];
               if length == 9
                   Y = [1, Nx, Ny, Nz, Nx*Ny, Nx*Nz, Ny*Nz, Nx*Nx-Ny*Ny, 3*Nz*Nz-1];
               else
                   Y = [1, Nx, Ny, Nz];
               end
               
               %Y = [1, -p/N, -q/N, 1/N];
               
               ComputedImage(i, j) = max(0, albedo*dot(Y,LightDirV));
           end
       end
    end
end

imshow(ComputedImage);
save ComputedImage.mat ComputedImage;

% deltaImage = zeros(image_width, image_height);
% for i = 1:image_width
%    for j = 1:image_height
%        if(maskImage(i,j) > 0.5)
%            deltaImage(i, j) = image(i, j) - computedImage(i, j);
%        end
%    end
% end
% 
% imshow(deltaImage);

end

