clear all
close all
%1.
load('data_lab3.mat')
for i=1:195
    im_vec(:,i) = reshape(im(:,:,i),101.^2,1);
end

im_cov = cov(im_vec);
im_mean = mean(im_vec);
im_eig = eigs(im_cov,5);
%plot eigenvalues
figure
plot(im_eig,':ok')
title('Eigenvalues in Descending Order')

%eigenvectors, eigenvalues
[PCVec,PCVal]=eigs(im_cov,3);
%pca
PCScore1 = PCVec(:,1).'*((im_vec(1,:)-im_mean)).';
for i = 2:10201
    PCScore1 = [PCScore1;PCVec(:,1).'*((im_vec(i,:)-im_mean)).'];
end

PCScore2 = PCVec(:,2).'*((im_vec(1,:)-im_mean)).';
for i = 2:10201
    PCScore2 = [PCScore2;PCVec(:,2).'*((im_vec(i,:)-im_mean)).'];
end

PCScore3 = PCVec(:,3).'*((im_vec(1,:)-im_mean)).';
for i = 2:10201
    PCScore3 = [PCScore3;PCVec(:,3).'*((im_vec(i,:)-im_mean)).'];
end

%%
figure
scatter(PCScore1,PCScore2)
title('PCA Scores of the first two components')
figure
scatter3(PCScore1.',PCScore2.',PCScore3.')
title('PCA Scores of the first three components');

%%
%plot pixels

%figure
%hold on
%for i = 1:101
%    for j = 1:101
%     plot(squeeze(im(i,j,:)))
%    end
%end

%reshape of PCA, show images of PCA
%for i = 1:3
%    figure
%    imshow(mat2gray(reshape(PCVec(:,i),1,195)),[])
%end

figure
for i = 1:3
    figure
    plot(PCVec(:,i))
end


%display the first 3 eigenvectors
figure
imshow(mat2gray(reshape(PCScore1,101,101)),[]);
figure
imshow(mat2gray(reshape(PCScore2,101,101)),[]);
figure
imshow(mat2gray(reshape(PCScore3,101,101)),[]);

%display with RGB bands
figure
imshow(mat2gray(im(:,:,[30,20,7])),[])


%%
%pca scores for the last few components
[PCVec_all,PCVal_all]=eigs(im_cov,194);
%for i = 192:194
%   figure
%   imshow(mat2gray(reshape(PCVec_all(:,i),1,195)),[])
%end

figure
hold on
for i = 192:194
    plot(PCVec_all(:,i))
end
hold off

%last pca score
PCScore_l1 = PCVec_all(:,194).'*((im_vec(1,:)-im_mean)).';
for i = 2:10201
    PCScore_l1 = [PCScore_l1;PCVec_all(:,194).'*((im_vec(i,:)-im_mean)).'];
end

PCScore_l2 = PCVec_all(:,193).'*((im_vec(1,:)-im_mean)).';
for i = 2:10201
    PCScore_l2 = [PCScore_l2;PCVec_all(:,193).'*((im_vec(i,:)-im_mean)).'];
end

PCScore_l3 = PCVec_all(:,192).'*((im_vec(1,:)-im_mean)).';
for i = 2:10201
    PCScore_l3 = [PCScore_l3;PCVec_all(:,192).'*((im_vec(i,:)-im_mean)).'];
end

figure
imshow(mat2gray(reshape(PCScore_l1,101,101)),[]);
figure
imshow(mat2gray(reshape(PCScore_l2,101,101)),[]);
figure
imshow(mat2gray(reshape(PCScore_l3,101,101)),[]);

%%
%54. Plot eigenvectors


for i = 1:3
    figure
    plot(PCVec(:,i))
end



for i = 1:3
    figure
    plot(S(:,i))
end

%%
% 9. projection of S
S_mean = mean(S.');
S_mean = S_mean.';

for i = 1:3
    coord1_S(i) =(PCVec(:,1).'*(S(:,i)-S_mean));
end

for i = 1:3
    coord2_S(i) = (PCVec(:,2).'*(S(:,i)-S_mean));
end

for i = 1:3
    coord3_S(i) = (PCVec(:,3).'*(S(:,i)-S_mean));
end

figure
hold on
scatter(PCScore1.',PCScore2.')
scatter(coord1_S.',coord2_S.')
title('PCA Scores with Coordinates of Columns of S')
hold off

figure
scatter3(PCScore1.',PCScore2.',PCScore3.') 
hold on
scatter3(coord1_S.',coord2_S.',coord3_S.')
title('PCA Scores with Coordinates of Columns of S')
hold off


%%
%10. Calculate A
X = im_vec.';
A = (inv(S.'* S))*S.'*X;
%plot abundance
for i = 1:3
    figure
    imshow(mat2gray(reshape(A(i,:),101,101)),[])
end

%%
%11.
tic
B_C = [S.'*S, 1*ones(3,1)];
D_E = [ones(1,3),0];
B_C_D_E = [B_C;D_E];

F = S.'*X;
G = [ones(1,10201)];
F_G = [F;G];

A_lambdatrans = inv(B_C_D_E)*F_G;
sol_A = [A_lambdatrans(1,:);A_lambdatrans(2,:);A_lambdatrans(3,:)];
toc
%plot abundance
for i = 1:3
    figure
    imshow(mat2gray(reshape(sol_A(i,:),101,101)),[]) 
end


%%
%12.
% Projection X
for i = 1:2
    X_proj(:,i) =X.'*PCVec(:,i);
end

%Projection of S
for i = 1:2
    S_proj(i,:) = S.'*PCVec(:,i);
end

%calculat new A
tic
Bp_Cp = [S_proj.'*S_proj, -1*ones(3,1)];
Dp_Ep = [ones(1,3),0];
Bp_Cp_Dp_Ep = [Bp_Cp;Dp_Ep];

Fp = S_proj.'*X_proj.';
Gp = [ones(1,10201)];
Fp_Gp = [Fp;Gp];

proj_A_lambdatrans = inv(Bp_Cp_Dp_Ep)*Fp_Gp;
proj_A = [proj_A_lambdatrans(1,:);proj_A_lambdatrans(2,:);proj_A_lambdatrans(3,:)];
toc

%%
%13. RMSE
%unconstrained
for i = 1:10201
    RMSE_unc(i)= (1/sqrt(195))*norm(X(:,i)-S*A(:,i));
end

%constrained
for i = 1:10201
    RMSE_con(i)= (1/sqrt(195))*norm(X(:,i)-S*sol_A(:,i));
end

figure
imshow(mat2gray(reshape(RMSE_unc,101,101)),[]);
figure
imshow(mat2gray(reshape(RMSE_con,101,101)),[]);

%compute mean
RMSE_unc_mean = mean(RMSE_unc);
RMSE_con_mean = mean(RMSE_con);
