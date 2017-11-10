function [label, centroids_value] = lmwkmeans( X,K)
%UNTITLED14 Summary of this function goes here
%   Detailed explanation goes here
max_iterations = 1000;
centroids = initCentroids(X, K);

[indices_old,centroids_old] = ind_centr(X,K )
for i=2:max_iterations
  [indices,centroids] = ind_centr(X,K )
  if (((centroids-centroids_old)/centroids_old)<0.01)
      centroids
      break;
  else
      centroids_old = centroids;
  end

end

label = indices;
centroids_value = centroids;


end

