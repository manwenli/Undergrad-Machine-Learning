clear all
close all
load temperature_data

%%
%1. linear fit
phi = [ones(137,1),x];
w = inv((transpose(phi)*phi))*transpose(phi)*t;
er = t - phi*w;
figure
scatter(x,t)
hold on
plot(x,phi*w,'linewidth',2)
legend('data','regression line')
title('Linear Regression Fit')
hold off

figure 
scatter(x,er)
hold on
refline(0,0)
title('Error')
hold off

%%
%2. extrapolate
x_ext = [x.',[2017:2100]].';
x_extra = [2017:2100].';
t_ext = w(1) + w(2)*x_ext;
t_extra = w(1) + w(2)*x_extra;
figure
scatter(x,t)
hold on
scatter(x_extra,t_extra,'filled')
hold on
plot(x_ext,t_ext,'Linewidth',2)
legend('Data','Prediction starting 2017','Regression')
title('Extrapolated Data')
hold off

%3. kernel function

%%
%4.realizations of Gaussian prior
sample_x = linspace(1880,2100,1000);
sample_x = sample_x.';

%compute K
for i = 1:1000
    for j = 1:1000
        K_r(i,j) = compKernel(sample_x(i),sample_x(j),10);
    end
end

sample_realization = mvnrnd(zeros(1000,1),K_r,3);

figure
plot(sample_x,sample_realization,'linewidth',1.5)

title('Realizations of GP Prior')


%%
%5.Mean and variance
%C = K + ?^2 * I_N

%compute K
for i = 1:137
    for j = 1:137
        K(i,j) = compKernel(x(i),x(j),10);
    end
end

C = K + 0.01 * eye(137);

for i = 1:1000
    c(i) = compKernel(sample_x(i),sample_x(i),10)+0.01;
end
c = c.';

for i = 1:1000
    for j = 1:137
        k(j,i) = compKernel(x(j),sample_x(i),10);
    end
end

for i = 1:1000
    mean(i) = k(:,i).'*inv(C)*t;
    sig_sq(i) = c(i)-k(:,i).'*inv(C)*k(:,i);
end

sig_sq = sig_sq.';

for i = 1:1000
    mean_plus(i) = mean(i)+2*sqrt(sig_sq(i));
    mean_minus(i) = mean(i)-2*sqrt(sig_sq(i));
end


figure
scatter(x,t)
hold on
plot(sample_x,mean,'linewidth',2)
hold on
plot(sample_x,mean_plus,'linewidth',2)
hold on
plot(sample_x,mean_minus,'linewidth',2)
title('Predictive Distribution')
legend('Data','Mean,Predictive','Mean+2SD','Mean-2SD')
hold off

%%
%6. realizations
% tao = 10
for i = 1:1000
    for j = 1:1000
        K_1(i,j) = compKernel(sample_x(i),sample_x(j),10);
    end
end

r_1 = mvnrnd(zeros(1000,1),K_1,1);

% tao = 1000
for i = 1:1000
    for j = 1:1000
        K_2(i,j) = compKernel(sample_x(i),sample_x(j),1000);
    end
end

r_2 = mvnrnd(zeros(1000,1),K_2,1);


% tao = 10000
for i = 1:1000
    for j = 1:1000
        K_3(i,j) = compKernel(sample_x(i),sample_x(j),10000);
    end
end

r_3 = mvnrnd(zeros(1000,1),K_3,1);


% tao = 100000
for i = 1:1000
    for j = 1:1000
        K_4(i,j) = compKernel(sample_x(i),sample_x(j),100000);
    end
end

r_4 = mvnrnd(zeros(1000,1),K_4,1);



%%
%6. change values of tao
% tao = 10
for i = 1:137
    for j = 1:137
        K_10(i,j) = compKernel(x(i),x(j),10);
    end
end

sample_10 = mvnrnd(zeros(137,1),K_10,1);

% tao = 1000
for i = 1:137
    for j = 1:137
        K_1000(i,j) = compKernel(x(i),x(j),1000);
    end
end

sample_1000 = mvnrnd(zeros(137,1),K_1000,1);


% tao = 10000
for i = 1:137
    for j = 1:137
        K_10000(i,j) = compKernel(x(i),x(j),10000);
    end
end

sample_10000 = mvnrnd(zeros(137,1),K_10000,1);


% tao = 100000
for i = 1:137
    for j = 1:137
        K_100000(i,j) = compKernel(x(i),x(j),100000);
    end
end

sample_100000 = mvnrnd(zeros(137,1),K_100000,1);

[mymean1,myvar1] = mean_var(K_10, sample_x,10, x, t);
[mymean2,myvar2] = mean_var(K_1000, sample_x,1000, x, t);
[mymean3,myvar3] = mean_var(K_10000, sample_x,10000, x, t);
[mymean4,myvar4] = mean_var(K_100000, sample_x,100000, x, t);



%%
figure
scatter(sample_x,r_1,'filled')
hold on
scatter(sample_x,r_2,'filled')
hold on
scatter(sample_x,r_3,'filled')
hold on
scatter(sample_x,r_4,'filled')
legend('\tau^2 = 10','\tau^2 = 1000','\tau^2 = 10000','\tau^2 = 100000')
title('Realizations of GP Prior with Different \tau^2 values')
hold off

figure
plot(sample_x, mymean1,'linewidth',2)
hold on
plot(sample_x,mymean2,'linewidth',2)
hold on
plot(sample_x,mymean3,'linewidth',2)
hold on
plot(sample_x,mymean4,'linewidth',2)
legend('\tau^2 = 10','\tau^2 = 1000','\tau^2 = 10000','\tau^2 = 100000')
title('Mean of the Predictive Distribution with Different \tau^2 values')
hold off

%%
%7.
mymean3 = mymean3.';
mean_p_sd = mymean3+2*sqrt(myvar3);
mean_m_sd = mymean3-2*sqrt(myvar3);

figure
scatter(x,t)
hold on
plot(sample_x,mymean3,'linewidth',2)
plot(sample_x,mean_p_sd,'linewidth',2)
plot(sample_x,mean_m_sd,'linewidth',2)
legend('original data','mean','mean+2SD','mean-2SD')
title('Mean of the Predictive Distribution with \tau^2 = 10000')
hold off

%%
%8.
sample_new = linspace(1880,2200,1000);
sample_new = sample_new.';
% tao = 10000
for i = 1:137
    for j = 1:137
        Kn_10000(i,j) = compKernel(x(i),x(j),10000);
    end
end

sample_n10000 = mvnrnd(zeros(137,1),Kn_10000,1);


% tao = 100000
for i = 1:137
    for j = 1:137
        Kn_100000(i,j) = compKernel(x(i),x(j),100000);
    end
end

sample_n100000 = mvnrnd(zeros(137,1),Kn_100000,1);

[mymean3_n,mysample3_n] = mean_var(Kn_10000, sample_new,10000, x, t);
[mymean4_n,mysample4_n] = mean_var(Kn_100000, sample_new,100000, x, t);
%%
figure
scatter(x,t)
hold on
plot(sample_new,mymean3_n,'linewidth',2)
hold on
plot(sample_new,mymean4_n,'linewidth',2)
legend('original','\tau^2=10000','\tau^2=100000')
title('Mean of the Predictive Distribution, Until Year 2200')
hold off
%%
figure
scatter(x,t)
hold on
plot(sample_x,mymean3,'linewidth',2)
legend('original','\tau^2=10000')
title('Mean of the Predictive Distribution, Until Year 2100')
hold off

%%
sample_local = linspace(1880,2027,1000);
sample_local = sample_local.';
% tao = 110
for i = 1:137
    for j = 1:137
        K_30(i,j) = compKernel(x(i),x(j),110);
    end
end

sample_30 = mvnrnd(zeros(137,1),K_30,1);


% tao = 500
for i = 1:137
    for j = 1:137
        K_100(i,j) = compKernel(x(i),x(j),700);
    end
end

sample_100 = mvnrnd(zeros(137,1),K_100,1);

[mymean30,mysample30] = mean_var(K_30, sample_local,110, x, t);
[mymean100,mysample100] = mean_var(K_100, sample_local,700, x, t);
%%
figure
scatter(x,t)
hold on
plot(sample_local,mymean30,'linewidth',2)
hold on
plot(sample_local,mymean100,'linewidth',2)
legend('original','\tau^2=110','\tau^2=700')
title('Mean of the Predictive Distribution, Until Year 2027')
hold off
