function [ local_summed ] = arfs_local_sum( A, m, n )
% [ local_summed ] = arfs_local_sum( A, m, n )
% Takes the sum of the pixels overlapping at any given point. Assumes that
% the input image (A) is the same size as what is overlapping it.
%
% Used algorithm from "normxcorr2_general" by Dirk Padfield, 
% because it takes 1/3 of the time that the algorithm used by MATLAB does.
%
% Dirk Padfield. "Masked FFT registration". In Proc. Computer
%   Vision and Pattern Recognition, 2010.
%

    s = cumsum(A,1);
    c = [s; repmat(s(end,:),m-1,1) - s(1:end-1,:)];
    s = cumsum(c,2);
    local_summed = [s, repmat(s(:,end),1,n-1) - s(:,1:end-1)];

end

