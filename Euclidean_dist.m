function [outs] = Euclidean_dist(input)
%Generates an array containing the total euclidean distance generted from
%each barcode inreference to each other barcode within each replicate



%input is a raw matrix contiaing the reads of interest..
%rows = replicates, columns = 
%output will be a vector? of the distances
v = 1:size(input,2);
C = nchoosek(v, 2);
outs = zeros(size(input,1),size(C,1));
for i = 1:size(input,1)
    clr = log2(input(i,:)./geomean(input(i,:))); %clr transform
    parfor j = 1:size(C,1)
        outs(i,j) = sqrt((clr(C(j,1)) -  clr(C(j,2)))^2)
    end
end

    
