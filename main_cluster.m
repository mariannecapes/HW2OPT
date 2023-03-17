
clear all
close all
clc

addpath Data3D
addpath Algorithms
addpath Tools



%% load 3D mesh

name_data = 'dragon3D_true' ;
%name_data = 'bunny3D' ;

load([name_data, '.mat'])

% OUTPUT: 
% t, xbar: 3D mesh 
%          xbar contains the 3D positions of each node
%          t is only used to plot the 3D mesh
% numb_comm: 2 * number of edges in the graph
% degree: diagonal of the degree matrix
%         gives number of neighboors for each node
% neighbors: matrix giving, for each node, all its neighbors
%            to fill the rows, the value (nodes+1) is used
%            Size of the matrix: (nodes) x (max(degree))
% Weights: Weights for the TV norm (incidence matrix)
%          Size of the matrix: (nodes) x (beta_star)
% Ind_current/Ind_neighb: vector giving all the connections of the graph
%              compact version of the incidence matrix
%              Size of the vectors: (beta_star) x 1
%              Corresponds to the linear operato L in the TV
% Ladj:   matrix containing graph connections
%         adjoint operator of L for TV norm on graphs
%         Size of the matrix: (nodes) x (beta_star)

L =@(x) x(Ind_current,:) -x(Ind_neighb,:) ;
Lt =@(y) Ladj*y ;

%% Useful function for quality evaluation

l2_err =@(x) [ sum((xbar(:,1)-x(:,1)).^2)  sum((xbar(:,2)-x(:,2)).^2)  sum((xbar(:,3)-x(:,3)).^2) ] ;
RMSE =@(x) sqrt(sum(l2_err(x)));

%% Create noisy graph

% add Gaussian noise
var_noise = 3e-3 ; %variance of noise
z = xbar + var_noise*randn(nodes,3);
  
disp(['RMSE noisy mesh: ', num2str(RMSE(z))])



%% Algorithm parameters 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TO COMPLETE
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Squared norm of operator L
normL2 = op_norm(L,Lt,size(xbar)) ;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% bounds for the 3D coordinates of the graph
xmin = min(xbar,[],1) - 0.1*abs(min(xbar,[],1)) ;  
xmax = max(xbar,[],1) + 0.1*abs(max(xbar,[],1)) ;

% TV parameters: bun res 2
eta = 0.002 ;

% Maximum number of iterations
NbIt = 100; % To change depending on dataset size
% Display in the algorithm
display =  100;
display_graph =  0;
% Sotpping criteria
Stop_norm = 1e-3 ; % To change depending on dataset size
Stop_crit = 1e-3 ; % To change depending on dataset size

% iniialization
x0 = proj_box(z,xmin,xmax) ;

%% Primal-dual for 3D mesh denoising


disp('******************************************')
disp('Primal-dual algorithm')
[xres,crit,time,rmse,norm_x] = ...
    Primal_Dual_graph(x0, z, Ind_current,Ind_neighb,eta, Ladj, Neighb_mat, xmin,xmax, normL2, NbIt,t, RMSE, Stop_norm, Stop_crit, display, display_graph) ;

%% FB algorithm

disp('******************************************')
disp('FB algorithm')
[xFB, critFB, norm_xFB, rmseFB, timeFB] = ...
    FB(x0, z, Ind_current,Ind_neighb,eta, Ladj, Neighb_mat, xmin,xmax,normL2, NbIt,t, RMSE, Stop_norm, Stop_crit, display, display_graph) ;


%% FISTA

disp('******************************************')
disp('FISTA')
[xFISTA, critFISTA, norm_xFISTA, rmseFISTA, timeFISTA] = ...
    FISTA(x0, z, Ind_current,Ind_neighb,eta, Ladj, Neighb_mat, xmin,xmax,normL2, NbIt,t, RMSE, Stop_norm, Stop_crit, display, display_graph) ;


% %% Primal-dual for 3D mesh denoising (constrained problem)
% 
% % Maximum number of iterations
% l2_bound = 0.4 ;
% disp('******************************************')
% disp('Primal-dual algorithm for constrained problem')
% [xres_c,crit_c, norm_fid_c,time_c,rmse_c,norm_x_c] = ...
%     Primal_Dual_constrained_graph(x0, z, l2_bound, Ind_current,Ind_neighb,eta, Ladj, Neighb_mat, xmin,xmax, normL2, NbIt,t, RMSE, Stop_norm, Stop_crit, display, display_graph) ;
% disp('*******************************************************')
% disp('*******************************************************')
% 

%% Save results

save(['Results/',name_data, '_res.mat'], ...
        't', 'z', 'xbar', ...
        'xres', 'crit', 'time', 'rmse', 'norm_x', ...
        'xFB', 'critFB', 'norm_xFB', 'rmseFB', 'timeFB', ...
        'xFISTA', 'critFISTA', 'norm_xFISTA', 'rmseFISTA', 'timeFISTA', ...
        'xres_c', 'crit_c', 'norm_fid_c', 'norm_x_c', 'rmse_c', 'time_c')


