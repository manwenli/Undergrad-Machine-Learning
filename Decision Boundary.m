clear all
close all

%%
%1. Generate
mu_1 = [1;2];
mu_2 = [2;4];
sigma_1 = [[0.1,0.05];[0.05,0.2]];
sigma_2 = [[0.2,-0.1];[-0.1,0.3]];

rng(940730)
train_1 = mvnrnd(mu_1,sigma_1,50);
train_2 = mvnrnd(mu_2,sigma_2,50);


%%
%2. Scatterplot
figure
hold on
scatter(train_1(:,1),train_1(:,2),'b')
scatter(train_2(:,1),train_2(:,2),'r')
title('Scatterplot of Training Data')
hold off

%%
%3. Compute W
X = [train_1;train_2];
X_tilde = [ones(100,1),X];
T_1 = [ones(50,1),zeros(50,1)];
T_2 = [zeros(50,1),ones(50,1)];
T = [T_1;T_2];
W_tilde = (inv(X_tilde.' * X_tilde))*X_tilde.'*T;
X_W_tilde = X_tilde * W_tilde;

figure
plot(X_W_tilde)
title('X^~ * W^~')

figure
plot(T)
title('Target Matrix T')

%4. plot decision boundary
W = [W_tilde(2,:); W_tilde(3,:)];
W0 = [W_tilde(1,:)];
x2_1 = (W0(1,2)-W0(1,1)-(W(1,1)-W(1,2))*0.75)/(W(2,1)-W(2,2));
x2_2 = (W0(1,2)-W0(1,1)-(W(1,1)-W(1,2))*2.5)/(W(2,1)-W(2,2));

figure
hold on
scatter(train_1(:,1),train_1(:,2),'r')
scatter(train_2(:,1),train_2(:,2),'b')
plot([0.75,2.5],[x2_1,x2_2],'black','LineWidth',2)
title('Training Sets with Decision Boundary')
hold off

%%
%5.
mu_21 = [2;2];
sig_21 = [[0.2,0.05];[0.05,0.3]];
mu_22 = [2;4];
sig_22 = [[0.4,-0.1];[-0.1,0.3]];
mu_23 = [3;3];
sig_23 = [[0.5,-0.3];[-0.3,0.4]];
rng(1234)
train_21 = mvnrnd(mu_21,sig_21,50);
train_22 = mvnrnd(mu_22,sig_22,50);
train_23 = mvnrnd(mu_23,sig_23,50);

%Compute W
X2 = [train_21;train_22;train_23];
X_tilde2 = [ones(150,1),X2];
T_21 = [ones(50,1),zeros(50,1),zeros(50,1)];
T_22 = [zeros(50,1),ones(50,1),zeros(50,1)];
T_23 = [zeros(50,1),zeros(50,1),ones(50,1)];
T2 = [T_21;T_22;T_23];
W_tilde2 = (inv(X_tilde2.' * X_tilde2))*X_tilde2.'*T2;
X_W_tilde2 = X_tilde2 * W_tilde2;

%Scatterplot
figure
hold on
scatter(train_21(:,1),train_21(:,2),'b')
scatter(train_22(:,1),train_22(:,2),'r')
scatter(train_23(:,1),train_23(:,2))
title('Scatterplot of Training Data')
hold off

%plot decision boundary
%group 1 and 2
W2 = [W_tilde2(2,:); W_tilde2(3,:)];
W02 = [W_tilde2(1,:)];
p_1 = (W02(1,2)-W02(1,1)-(W2(1,1)-W2(1,2))*0.75)/(W2(2,1)-W2(2,2));
p_2 = (W02(1,2)-W02(1,1)-(W2(1,1)-W2(1,2))*2.3)/(W2(2,1)-W2(2,2));

%group 1 and 3
p_3 = (W02(1,3)-W02(1,1)-(W2(1,1)-W2(1,3))*2.3)/(W2(2,1)-W2(2,3));
p_4 = (W02(1,3)-W02(1,1)-(W2(1,1)-W2(1,3))*4)/(W2(2,1)-W2(2,3));

%group 2 and 3
p_5 = (W02(1,3)-W02(1,2)-(W2(1,2)-W2(1,3))*2.3)/(W2(2,2)-W2(2,3));
p_6 = (W02(1,3)-W02(1,2)-(W2(1,2)-W2(1,3))*4)/(W2(2,2)-W2(2,3));

figure
hold on
scatter(train_21(:,1),train_21(:,2),'r')
scatter(train_22(:,1),train_22(:,2),'b')
scatter(train_23(:,1),train_23(:,2));
plot([0.75,2.3],[p_1,p_2],'black','LineWidth',2)
plot([2.3,4],[p_3,p_4],'black','LineWidth',2)
plot([2.3,4],[p_5,p_6],'black','LineWidth',2)



title('Training Sets with Decision Boundary')
hold off

%%
%7. M=100
rng(1234);
train_31 = mvnrnd(mu_21,sig_21,100);
train_32 = mvnrnd(mu_22,sig_22,100);
train_33 = mvnrnd(mu_23,sig_23,100);

X3 = [train_31;train_32;train_33];
X_tilde3 = [ones(300,1),X3];
T_31 = [ones(100,1),zeros(100,1),zeros(100,1)];
T_32 = [zeros(100,1),ones(100,1),zeros(100,1)];
T_33 = [zeros(100,1),zeros(100,1),ones(100,1)];
T3 = [T_31;T_32;T_33];
W_tilde3 = (inv(X_tilde3.' * X_tilde3))*X_tilde3.'*T3;
X_W_tilde3 = X_tilde3 * W_tilde3;

%calculate R
R = 0;
for i = 1:100
    [M,I]=max(X_W_tilde3(i,:));
    if I == 1
        R = R+1;
    end
end

for i = 101:200
    [M,I]=max(X_W_tilde3(i,:));
    if I == 2
        R = R+1;
    end
end

for i = 201:300
    [M,I]=max(X_W_tilde3(i,:));
    if I == 3
        R = R+1;
    end
end

c = 100*R/(3*100);

%%
%8. K-NN

rng(1234);
train_41 = mvnrnd(mu_21,sig_21,100);
train_42 = mvnrnd(mu_22,sig_22,100);
train_43 = mvnrnd(mu_23,sig_23,100);

X4 = [train_41;train_42;train_43];



labels = [ones(1,100),ones(1,100)+1,ones(1,100)+2];

[Class] = cvKnn(X3,X3,labels, 40);

temp = 0;
for i = 1:300
    if Class(i) == labels(i)
        temp = temp+1;
    end
end

