function [  ] = Generate(  )
%GENERATE Summary of this function goes here
%   Detailed explanation goes here

%Intensity = 0.2989 * RM + 0.5870 * GM + 0.1140 * BM; 

image_width = 60*4;
image_height = 60*4;
ball_radius = 28*4;
ball_posX = 30*4;
ball_posY = 30*4;

RefNormalizedNormalsV = zeros(1, image_width*image_height*3);
RefHeightMat = zeros(image_width, image_height);
RefNormalPQV = zeros(1, image_width*image_height*2);
RefMaskImage = zeros(image_width, image_height);
RefLightDir = [-1, 1, 1];

normal_xy = zeros(1, 3);
for i = 1:image_width
    for j = 1:image_height
        if (i-ball_posX)^2+(j-ball_posY)^2 <= ball_radius^2
            normal_xy(1) = i-ball_posX;
            normal_xy(2) = j-ball_posY;
            normal_xy(3) = sqrt(ball_radius^2-normal_xy(1)^2-normal_xy(2)^2);
            
            RefHeightMat(i, j) = normal_xy(3);
            
            if norm(normal_xy)>1e-10
                normal_xy = normal_xy/norm(normal_xy);
            end
            
            RefNormalizedNormalsV(((i-1)*image_height+j-1)*3+1) = normal_xy(1);
            RefNormalizedNormalsV(((i-1)*image_height+j-1)*3+2) = normal_xy(2);
            RefNormalizedNormalsV(((i-1)*image_height+j-1)*3+3) = normal_xy(3);
            
            RefMaskImage(i, j) = 1;
        else
            RefMaskImage(i, j) = 0;
        end
    end
end

RefMaskImage = LabelMaskImage(RefMaskImage);

for i = 1:image_width-1
    for j = 1:image_height-1
        if(RefMaskImage(i, j) > 0.5)
            RefNormalPQV(((i-1)*image_height+j-1)*2+1) = RefHeightMat(i+1, j) - RefHeightMat(i, j);
            RefNormalPQV(((i-1)*image_height+j-1)*2+2) = RefHeightMat(i, j+1) - RefHeightMat(i, j);
        end
    end
end

%save image_normal.mat normal_ori;

albedo = 0.8;
RefImage = zeros(image_width, image_height);
RefNormalizedLightDir = RefLightDir/norm(RefLightDir);
for i = 1:image_width
    for j = 1:image_height
        RefImage(i,j) = max(0, albedo*(RefNormalizedNormalsV(((i-1)*image_height+j-1)*3+1)*RefNormalizedLightDir(1)+RefNormalizedNormalsV(((i-1)*image_height+j-1)*3+2)*RefNormalizedLightDir(2)+RefNormalizedNormalsV(((i-1)*image_height+j-1)*3+3)*RefNormalizedLightDir(3)));
    end;
end;

imshow(RefImage), title('RefImage');

save RefImage.mat RefImage;
save RefMaskImage.mat RefMaskImage;
save RefNormalizedNormalsV.mat RefNormalizedNormalsV;
save RefNormalPQV.mat RefNormalPQV;
save RefHeightMat.mat RefHeightMat;


end

