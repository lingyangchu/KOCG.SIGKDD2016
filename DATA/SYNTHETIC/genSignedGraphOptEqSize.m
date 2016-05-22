function [A, GT_X] = genSignedGraphOptEqSize(NV, NG, sparse_thres)
%% this function generates signed graphs with equal group size
% A: the random permutated matrix of A_org
% NV: number of vertexes
% NG: number of groups
% sparse_thres: the percentage of zero elemetns in the full matrix.
%% code as follow
fprintf('Start GenSignedGraphOPT ... \n');

if NG >= 2
    step = ceil(1.0*NV/NG);
    poses = unique([1:step:NV NV]);
else
    fprintf('Warning: Please generate at least 2 clusters! \n');
    assert(NG >= 2);
end

num_clusters = length(poses)-1;

fprintf('Num of clusters: %d \n', length(poses)-1);

clus_idx = zeros(NV, 1);
for i = 2:length(poses)    
    range = [(poses(i-1)):poses(i)];
    clus_idx(range) = i-1;
end

% build edges
m = ceil(NV^2*(1 - sparse_thres)); % num of non-zero edges to sample
ind = randperm(NV^2, m);
[x, y] = ind2sub([NV, NV], ind);
v = sign(double(clus_idx(x) == clus_idx(y)) - 0.5);
A = sparse(x, y, v, NV, NV);
A = sign(A+A');
A(logical(eye(NV))) = 0;

% avoid non-connected nodes
nc_x = find(sum(abs(A), 2) == 0);
if (~isempty(nc_x))
    nc_y = randperm(NV, length(nc_x))';
    nc_v = ones(length(nc_x), 1)*0.1;
    add_A = sparse(nc_x, nc_y, nc_v, NV, NV);
    add_A = add_A + add_A';
    add_A(logical(eye(NV))) = 0;
    
    % update A to eleminate non-connected nodes
    A = A + add_A;
    A = sign(A);
end

GT_X = recoverX(clus_idx);

fprintf('GenSignedGraphOPT ... Done! \n');






%% ************************************************************************
%% **************************  UTILITY FUNCTIONS  *************************
%% ************************************************************************

function X = recoverX(clus_result)
% recover clus_result into X (i.e., the NV-by-NG dimensional indicator matrix)
%% code as follow

clus_idx = unique(clus_result);

real_NG = length(clus_idx);

NV = length(clus_result);

X = sparse(NV, real_NG);

for i = 1:real_NG
    indi = (clus_result == clus_idx(i));
    
    X(indi, i) = 1;
end