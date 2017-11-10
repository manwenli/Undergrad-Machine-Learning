function [ mean, sig_sq] = mean_var(K, sample_x, tao_sq,x,t )
C = K + 0.01 * eye(137);

for i = 1:1000
    c(i) = compKernel(sample_x(i),sample_x(i),tao_sq)+0.01;
end
c = c.';

for i = 1:1000
    for j = 1:137
        k(j,i) = compKernel(x(j),sample_x(i),tao_sq);
    end
end

for i = 1:1000
    mean(i) = k(:,i).'*inv(C)*t;
    sig_sq(i) = c(i)-k(:,i).'*inv(C)*k(:,i);
end

sig_sq = sig_sq.';




end

