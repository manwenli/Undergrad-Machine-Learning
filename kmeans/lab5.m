clear all
close all
load royce_hall_small
royce = reshape(I,96*144,3);

%%
%1. kmean function
[index,centers] = kmean(royce,10);
testimage = zeros(size(royce));
for i = 1: 10
    j = sum (index == i);
testimage(index==i,:) = repmat(centers(i,:),j,1);
end
testimage1 = reshape(testimage,96,144,3);

figure
imshow(testimage1);

%%
%2.Test the algorithm for different values of K.
for k1 = 2:5
    [index,centers] = kmean(royce,k1);
    segIm = segmentation(royce,centers,index,k1);
    figure
    imshow(segIm)
end


for i = 1:5
    k2 = i*10;
    [index,centers] = kmean(royce,k2);
    segIm = segmentation(royce,centers,index,k2);
    figure
    imshow(segIm)
end
    
%%
%3. 3D Scatterplots of data
[index3,centers3] = kmean(royce,4);

figure
scatter3(royce(:,1),royce(:,2),royce(:,3))
hold on
scatter3(centers3(:,1),centers3(:,2),centers3(:,3),300,'filled')
title('Scatterplot of Data with Centroids, K = 4')
hold off

%%
%4.add zero mean gaussian noise
noisex = sqrt((var(royce(:,1)))/10)*randn(13824,1);
noisey = sqrt((var(royce(:,2)))/10)*randn(13824,1);
noisez = sqrt((var(royce(:,3)))/10)*randn(13824,1);
noise_mat = [noisex,noisey,noisez];
noise_image = royce+noise_mat;

[labels4,centr4] = kmean(noise_image,20);

  segIm4 = segmentation(noise_image,centr4,labels4,20);
  figure
  %imshow(segIm4)
  imshow(reshape(noise_image,96,144,3))

%%
%6.RMSE

for K6 = 1:100
   [idx6,c6]=kmeans(royce,K6); 
   segIm6 = segmentation(royce,c6,idx6,K6);
   segVec6 = reshape(segIm6,96*144,3);
   
   container(:,:,K6) = segVec6;
end

%%
for i = 1:13824
    X = royce(i,:);
    X_h = container(i,:,2);
    RMSE(i) = norm(X-X_h,2)/sqrt(3);
end
RMSE = RMSE.';
imshow(mat2gray(reshape(RMSE,96,144)),[])

%%
%7.
royce_vec = reshape(royce,13824*3,1);
for K7=1:100
    seg_vec = reshape(container(:,:,K7),13824*3,1);
    RMSE_mean(K7) = norm(royce_vec - seg_vec,2)/sqrt(3);
end
figure
plot(RMSE_mean)
title('Mean RMSE')

%%
%8. mean RMSE vs. Compression ratio
for K8 = 1:100
    CR(K8) = (24*K8 + log2(K8)*13824)/(24*13824);
end
line1 = [[ones(36,1)-0.85],[0:35].'];
line2 = [[ones(36,1)-0.9],[0:35].'];
figure
hold on
scatter(CR,RMSE_mean)
plot(line1(:,1),line1(:,2),'black')
plot(line2(:,1),line2(:,2),'black')
title('Mean RMSE vs the Compression Ratio')
hold off

%%
%9.
[index9,centers9] = kmeans(royce,50);

figure
scatter3(royce(:,1),royce(:,2),royce(:,3))
hold on
scatter3(centers9(:,1),centers9(:,2),centers9(:,3),200,'filled')
title('Data with Centroids, K = 50')
hold off