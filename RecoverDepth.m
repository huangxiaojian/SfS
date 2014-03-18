function [  ] = RecoverDepth(  )
%RECOVERDEPTH Summary of this function goes here
%   Detailed explanation goes here


L = importdata('LightDirV.mat'); % result from recover_light step
I = importdata('RefImage.mat'); %input data
MaskImage = importdata('RefMaskImage.mat'); %input data
NormalPQ = importdata('RefNormalPQV.mat');
[image_width, image_height] = size(I);
Nref = zeros(image_width, image_height);
RefZ = importdata('RefHeightMat.mat'); %input data
ComputedImage = importdata('ComputedImage.mat'); % result from recover_light step

albedo = 0.8;

landa = 0.75;

countTotal = 0;%count2+count3
count1 = 0;%internal
count2 = 0;%boundary


%build masked index mapping, including internal and boundary
%Index Mat->index
%Index index->Mat(i,j)
MatXYToIndex = zeros(image_width, image_height);
for i = 1:image_width
    for j = 1:image_height
        if(MaskImage(i, j) > 0.5)
            
            Nref(i, j) = sqrt(NormalPQ(((i-1)*image_height+j-1)*2+1)^2+NormalPQ(((i-1)*image_height+j-1)*2+2)^2+1);
            
            countTotal = countTotal + 1;
            MatXYToIndex(i, j) = countTotal;
            IndexToMatX(countTotal) = i;
            IndexToMatY(countTotal) = j;
        end
    end
end


%build boundary index maaping
for i = 2:image_width
    for j = 2:image_height
        if(MaskImage(i, j) == 2) %boundary
            count2 = count2 + 1;
            BoundaryIndexX(count2) = i;
            BoundaryIndexY(count2) = j;
        end
    end
end

count1 = countTotal - count2; %internal

% for sparse
G = fspecial('gaussian',[3,3],3);
step2GArray = zeros(1, 9);
for k = 1:9
    step2GArray(k) = -G(k);
end
step2GArray(5) = 1 + step2GArray(5);


Ai1 = zeros(1, count1*3);
Aj1 = zeros(1, count1*3);
As1 = zeros(1, count1*3);
B1 = zeros(1, countTotal);
B11 = zeros(1, countTotal);
B12 = zeros(1, countTotal);

indexCount = 0;
for i = 2:image_width
    for j = 2:image_height
        if(MaskImage(i, j) == 1) % internal
            currentIndex = MatXYToIndex(i, j);
            B11(currentIndex) = ComputedImage(i, j); % precise value
            B12(currentIndex) = - albedo * L(1)  - albedo/Nref(i, j) * L(4);       

            B1(currentIndex) = ComputedImage(i, j) - albedo * L(1) - albedo/Nref(i, j) * L(4);       

            indexCount = indexCount + 1;
            Ai1(indexCount) = currentIndex;
            Aj1(indexCount) = MatXYToIndex(i, j);
            As1(indexCount) = -albedo/Nref(i, j) * (L(2)+L(3));%x,y
            indexCount = indexCount + 1;
            Ai1(indexCount) = currentIndex;
            Aj1(indexCount) = MatXYToIndex(i+1, j);
            As1(indexCount) = -albedo/Nref(i, j) * L(2);%x+1,y
            indexCount = indexCount + 1;
            Ai1(indexCount) = currentIndex;
            Aj1(indexCount) = MatXYToIndex(i, j+1);
            As1(indexCount) = -albedo/Nref(i, j) * L(3);%x,y+1
        end
    end
end


A1 = sparse(Ai1, Aj1, As1, countTotal, countTotal);

Ai2 = zeros(1, count1*9);
Aj2 = zeros(1, count1*9);
As2 = zeros(1, count1*9);
B2 = zeros(1, countTotal);

indexCount = 0;
for i = 2:image_width
    for j = 2:image_height
        if(MaskImage(i, j) == 1)
            currentIndex = MatXYToIndex(i, j);
            
            step2IndexArrayX = [i-1, i, i+1, i-1, i, i+1, i-1, i, i+1];
            step2IndexArrayY = [j-1, j-1, j-1, j, j, j, j+1, j+1, j+1];
            
            btmp = RefZ(i, j);
            for k = 1:9
                indexCount = indexCount + 1;
                Ai2(indexCount) = currentIndex;
                Aj2(indexCount) = MatXYToIndex(step2IndexArrayX(k), step2IndexArrayY(k));
                As2(indexCount) = landa * albedo * step2GArray(k);
                
                btmp = btmp - G(k) * RefZ(step2IndexArrayX(k), step2IndexArrayY(k));
            end
            
            B2(currentIndex) = landa * btmp;
         
        end
    end
end

A2 = sparse(Ai2, Aj2, As2, countTotal, countTotal);

count3 = 0;

B3 = zeros(1, countTotal);

indexCount = 0;
breakFlag = 0;
for i = 2:image_width
    for j = 2:image_height
        if(MaskImage(i, j) == 2) %boundary
            currentIndex = MatXYToIndex(i, j);
            
            step2IndexArrayX = [i-1, i, i+1, i-1, i, i+1, i-1, i, i+1];
            step2IndexArrayY = [j-1, j-1, j-1, j, j, j, j+1, j+1, j+1];
            
            count3 = count3 + 1;
            if(count3 == count2)
                breakFlag = 1;
                    break;
            end
            
            btmp = RefZ(i, j);
            for k = 1:9
                if(MaskImage(step2IndexArrayX(k), step2IndexArrayY(k)) > 0.5)
                    indexCount = indexCount + 1;
                    Ai3(indexCount) = currentIndex;
                    Aj3(indexCount) = MatXYToIndex(step2IndexArrayX(k), step2IndexArrayY(k));
                    As3(indexCount) = landa * albedo * step2GArray(k);

                    btmp = btmp - G(k) * RefZ(step2IndexArrayX(k), step2IndexArrayY(k)); 
                end
            end
            
            B3(currentIndex) = landa * btmp;
         
        end
    end
    if(breakFlag == 1)
        break;
    end
end

indexCount = indexCount + 1;
Ai3(indexCount) = currentIndex;
Aj3(indexCount) = currentIndex;
As3(indexCount) = 1;

B3(currentIndex) = 0;

A3 = sparse(Ai3, Aj3, As3, countTotal, countTotal);

A = A1+A2+A3;

B = B11+B12+B2+B3;
B(currentIndex) = 0;

X = A \ B';

% X to HeightMap
HeightImage = zeros(image_width, image_height);
for k = 1:countTotal
    i = IndexToMatX(k);
    j = IndexToMatY(k);
    if(MaskImage(i, j) > 0.5)
        HeightImage(i, j) = X(k);
    else
        HeightImage(i, j) = 0;
    end
end

save HeightImageMat.mat HeightImage;

SaveObjMesh('ResultObj.obj', HeightImage);

mesh(HeightImage);
%figure('HeightImage');

end

