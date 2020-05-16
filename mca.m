function [mcamaps_l3d, pcs_l, mcamaps_r3d, pcs_r, lambdas] = mca(data_l, data_r, varargin)
% [mcamaps_l3d, pcs_l, mcamaps_r3d, pcs_r, lambdas] = mca(data_l, data_r, varargin)
% Version 1.0
% MCA or SVD Toolbox for marine and atmospheric science.
% Function names MCA but not SVD here because the MATLAB build-in "svd" function
%% Input:
%     data_l(lon1,lat1,time)
%     date_r(lon2,lat2,time)
%     The timesteps of two datasets should be the same.
% See also SVD
%%   Syntax
%       [mcamaps_l3d, pcs_l, mcamaps_r3d, pcs_r, lambdas] = mca(data_l, data_r)
% 
%       [mcamaps_l3d, pcs_l, mcamaps_r3d, pcs_r, lambdas] = mca(data_l, data_r, n_mca) 
%           Only calculate the first n_mca mode. 
%% Author:
%	Zelun Wu,
%   Ph.D. student of Physical Oceanography,
%	Xiamen University & University of Delaware
%	zelunwu@stu.xmu.edu.cn, zelunwu@udel.edu
%	15th May, 2020

  
%% Error checks 
narginchk(1,inf) 
%% Input parsing
[N_lon1, N_lat1, N_time1] = size(data_l);
[N_lon2, N_lat2, N_time2] = size(data_r);
assert(N_time1 == N_time2, 'Timesteps of two dataset should be the same');
if nargin > 2
    n_mca = varargin{1};
else
    error('Needs more than two dataset');
end
% core

[data_left_2d, in_notnan_l] = reshape3dto2d(data_l);
[data_right_2d, in_notnan_r] = reshape3dto2d(data_r);

[mcamaps_l, pcs_l, mcamaps_r, pcs_r, lambdas] = svdcore(data_left_2d, data_right_2d, n_mca);

mcamaps_l3d = reshape2dto3d(mcamaps_l, [N_lon1, N_lat1, n_mca], in_notnan_l);
mcamaps_r3d = reshape2dto3d(mcamaps_r, [N_lon2, N_lat2, n_mca], in_notnan_r);

end