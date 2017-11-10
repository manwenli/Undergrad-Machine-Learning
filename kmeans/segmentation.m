function segmentedIm = segmentation(X,Centroids,Index,K)
zero_image = zeros(size(X));
for i = 1: K
    j = sum (Index == i);
zero_image(Index==i,:) = repmat(Centroids(i,:),j,1);
end
segmentedIm = reshape(zero_image,96,144,3);


end

