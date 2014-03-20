function [  ] = PreProcess(  )
%PREPROCESS Summary of this function goes here
%   Detailed explanation goes here

ImageData = imread('face.png');
MaskImage = imread('maskImage.png');
RefHeightMat = load('newHeight.hm');
[imageWidth, imageHeight] = size(RefHeightMat);

RefImage = zeros(imageWidth, imageHeight);
RefNormalizedNormalsV = zeros(1, imageWidth*imageHeight*3);
RefNormalPQV = zeros(1, imageWidth*imageHeight*2); %internal
RefMaskImage = zeros(imageWidth, imageHeight);

MatXYToIndex = zeros(imageWidth, imageHeight);
MatXYToBoundary = zeros(imageWidth, imageHeight);


% mask mapping
maskCount = 0;
for i = 1:imageWidth
    for j = 1:imageHeight
        if MaskImage(i, j, 1) > 100         
            maskCount = maskCount+1;
            MatXYToIndex(i,j) = maskCount;
            IndexToMatX(maskCount) = i;
            IndexToMatY(maskCount) = j;
            RefMaskImage(i,j) = 1;
            r = 0.2989*double(ImageData(i,j,1));
            g = 0.5870*double(ImageData(i,j,2));
            b = 0.1140*double(ImageData(i,j,3));
            intensity = (r+g+b)/256;
            RefImage(i,j) = intensity;
        else
            RefMaskImage(i,j) = 0;
        end
    end
end

% label boundary
RefMaskImage = LabelMaskImage(RefMaskImage);

% boundary mapping
boundaryCount = 0;
for i = 1:imageWidth
    for j = 1:imageHeight
        if(RefMaskImage(i,j) == 2)
            boundaryCount = boundaryCount+1;
            MatXYToBoundary(i,j) = boundaryCount;
            BoundaryIndexX(boundaryCount) = i;
            BoundaryIndexY(boundaryCount) = j;
        end
    end
end

% internal normal
for k = 1:maskCount
    i = IndexToMatX(k);
    j = IndexToMatY(k);
    if RefMaskImage(i,j) < 1.5 %internal
        baseIndex = (i-1)*imageHeight+j-1;
        RefNormalPQV(baseIndex*2+1) = RefHeightMat(i+1, j) - RefHeightMat(i, j);
        RefNormalPQV(baseIndex*2+2) = RefHeightMat(i, j+1) - RefHeightMat(i, j);
        
        normal = [-RefNormalPQV(baseIndex*2+1), -RefNormalPQV(baseIndex*2+2), 1];
        normal = normal/norm(normal);
        
        RefNormalizedNormalsV(baseIndex*3+1) = normal(1);
        RefNormalizedNormalsV(baseIndex*3+2) = normal(2);
        RefNormalizedNormalsV(baseIndex*3+3) = normal(3);
    end
end

% boundary normal
for k = 1:boundaryCount
    i = BoundaryIndexX(k);
    j = BoundaryIndexY(k);
    RidxX = [i-1,i-1,i-1,i,i,i+1,i+1,i+1];
    RidxY = [j-1,j,j+1,j-1,j+1,j-1,j,j+1];
    
    count = 0;
    normal = [0, 0, 0];
    for n = 1:8 % or 4
        x = RidxX(n);
        y = RidxY(n);
        if RefMaskImage(x,y) > 0.5 && RefMaskImage(x,y) < 1.5 %internal
            count = count + 1;
            normal(1) = normal(1) + RefNormalizedNormalsV(((x-1)*imageHeight+y-1)*3+1);
            normal(2) = normal(2) + RefNormalizedNormalsV(((x-1)*imageHeight+y-1)*3+2);
            normal(3) = normal(3) + RefNormalizedNormalsV(((x-1)*imageHeight+y-1)*3+3);
        end
    end
    if count < 1 % count == 0
        fprintf('fuck(%d,%d)\n', i, j);
        LRidxX = [i-2,i-2,i-2,i-2,i-2,i-1,i-1,i,i,i+1,i+1,i+2,i+2,i+2,i+2,i+2];
        LRidxY = [j-2,j-1,j,j+1,j+2,j-2,j+2,j-2,j+2,j-2,j+2,j-2,j-1,j,j+1,j+2];
        for n = 1:16
            x = LRidxX(n);
            y = LRidxY(n);
            if RefMaskImage(x,y) > 0.5 && RefMaskImage(x,y) < 1.5 %internal
                count = count + 1;
                normal(1) = normal(1) + RefNormalizedNormalsV(((x-1)*imageHeight+y-1)*3+1);
                normal(2) = normal(2) + RefNormalizedNormalsV(((x-1)*imageHeight+y-1)*3+2);
                normal(3) = normal(3) + RefNormalizedNormalsV(((x-1)*imageHeight+y-1)*3+3);
            end
        end
        if count == 0
            fprintf('what the hell\n');
        end
    end
    normal = normal / count;
    normal = normal / norm(normal);
    RefNormalizedNormalsV(((i-1)*imageHeight+j-1)*3+1) = normal(1);
    RefNormalizedNormalsV(((i-1)*imageHeight+j-1)*3+2) = normal(2);
    RefNormalizedNormalsV(((i-1)*imageHeight+j-1)*3+3) = normal(3);
end

save MatXYToIndex.mat MatXYToIndex;
save MatXYToBoundary.mat MatXYToBoundary;
save IndexToMatXV.mat IndexToMatX;
save IndexToMatYV.mat IndexToMatY;
save BoundaryIndexXV.mat BoundaryIndexX;
save BoundaryIndexYV.mat BoundaryIndexY;

save RefImage.mat RefImage;
save RefMaskImage.mat RefMaskImage;
save RefNormalizedNormalsV.mat RefNormalizedNormalsV;
save RefNormalPQV.mat RefNormalPQV;
save RefHeightMat.mat RefHeightMat;

% imshow(RefImage);

RefLightDir = [0, 0, 1];
RefNormalizedLightDir = RefLightDir/norm(RefLightDir);
albedo = 0.8;
ComputedImage = zeros(imageWidth, imageHeight);
for i = 1:imageWidth
    for j = 1:imageHeight
        ComputedImage(i,j) = albedo*(RefNormalizedNormalsV(((i-1)*imageHeight+j-1)*3+1)*RefNormalizedLightDir(1)+RefNormalizedNormalsV(((i-1)*imageHeight+j-1)*3+2)*RefNormalizedLightDir(2)+RefNormalizedNormalsV(((i-1)*imageHeight+j-1)*3+3)*RefNormalizedLightDir(3));
    end
end

save ComputedImage.mat ComputedImage;
imshow(ComputedImage);

end