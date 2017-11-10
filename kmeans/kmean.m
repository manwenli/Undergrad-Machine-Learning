function [label, centroid_value] = kmean( X,K )
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here

max_iterations = 100000;
centroids = CentroidsInitialization(X, K);

indices_old = GroupAssignment(X, centroids);
centroids_old = CentroidsUpdate(X, indices_old, K);

for i=2:max_iterations
  indices = GroupAssignment(X, centroids_old);
  centroids = CentroidsUpdate(X, indices, K);
  if ((abs(norm(centroids)-norm(centroids_old))/norm(centroids_old))<0.01)
      centroids
      break;
  else
      centroids_old = centroids;
      
  end

end

label = indices;
centroid_value = centroids;

end

