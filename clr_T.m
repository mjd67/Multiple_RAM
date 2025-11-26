function [clr_mat] = clr_T(mat)
% Replicates must be by column, instances are rows
    clr_mat = zeros(size(mat));
    for i = 1:size(mat,2)
        gm = geomean(mat(:,i));
        clr_mat(:,i) = log2(mat(:,i) / gm);
    end
end