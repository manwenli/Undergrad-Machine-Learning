function [indices, centroids] = ind_centr(X,numClass )
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here
centroids = zeros(numClass,size(X,2)); 
    randidx = randperm(size(X,1));
    centroids = X(randidx(1:numClass), :);

    
 K = size(centroids, 1);
  indices = zeros(size(X,1), 1);
  m = size(X,1);

  for i=1:m
    k = 1;
    min_dist = sum((X(i,:) - centroids(1,:)) .^ 2);
    for j=2:K
        dist = sum((X(i,:) - centroids(j,:)) .^ 2);
        if(dist < min_dist)
          min_dist = dist;
          k = j;
        end
    end
    indices(i) = k;
  end
  
  [m n] = size(X);
  centroids = zeros(K, n);
  
  for i=1:K
    xi = X(indices==i,:);
    ck = size(xi,1);
    centroids(i, :) = (1/ck) * sum(xi);
  end
end

