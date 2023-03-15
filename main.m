
clear all
close all
clc

addpath('Data3D/')
addpath('Algorithms/')
addpath('Tools/')



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

%% Plot original graph and noisy graph

figure
subplot 121
display_3Dmesh(t,xbar)
xlabel('Initial graph')
subplot 122
display_3Dmesh(t,z)
xlabel('Noisy graph')


%% Algorithm parameters 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TO COMPLETE
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Squared norm of operator L
normL2 = op_norm(L, Lt, size(xbar));
%HELLO%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% bounds for the 3D coordinates of the graph
xmin = min(xbar,[],1) - 0.1*abs(min(xbar,[],1)) ;  
xmax = max(xbar,[],1) + 0.1*abs(max(xbar,[],1)) ;

% TV parameters: bun res 2
eta = 0.002 ;

% Maximum number of iterations
NbIt = 2000;
% Sotpping criteria
Stop_norm = 1e-5 ; 
Stop_crit = 1e-5 ;
% Display in the algorithm
display =  1000;
display_graph =  0;

% initialization
x0 = proj_box(z,xmin,xmax) ;

%% Primal-dual for 3D mesh denoising

disp('******************************************')
disp('Primal-dual algorithm')
[xres,crit,rmse,time] = ...
    Primal_Dual_graph(x0, z, Ind_current,eta, L, Lt, Neighb_mat, xmin,xmax, normL2, NbIt,t, RMSE, Stop_norm, Stop_crit, display, display_graph) ;

%% FB algorithm

disp('******************************************')
disp('FB algorithm')
[xFB, critFB, rmseFB, timeFB] = ...
    FB(x0, z, Ind_current,eta, L, Lt, Neighb_mat, xmin,xmax,normL2, NbIt,t, RMSE, Stop_norm, Stop_crit, display, display_graph) ;

%% FISTA

disp('******************************************')
disp('FISTA')
[xFISTA, critFISTA, rmseFISTA, timeFISTA] = ...
    FISTA(x0, z, Ind_current,eta, L, Lt, Neighb_mat, xmin,xmax,normL2, NbIt,t, RMSE, Stop_norm, Stop_crit, display, display_graph) ;


%% Primal-dual for 3D mesh denoising (constrained problem)

% Maximum number of iterations
l2_bound = 0.4 ;
disp('******************************************')
disp('Primal-dual algorithm for constrained problem')
[xres_c,crit_c, norm_fid_c,time_c,rmse_c,norm_x_c] = ...
    Primal_Dual_constrained_graph(x0, z, l2_bound, Ind_current,Ind_neighb,eta, Ladj, Neighb_mat, xmin,xmax, normL2, NbIt,t, RMSE, Stop_norm, Stop_crit, display, display_graph) ;
disp('*******************************************************')
disp('*******************************************************')


 
%% Plot results

figure
subplot 221
display_3Dmesh(t,z)
xlabel(['Noisy, RMSE = ', num2str(RMSE(z))])
subplot 222
display_3Dmesh(t,xres)
xlabel(['Primal dual, RMSE = ', num2str(RMSE(xres))])
subplot 223
display_3Dmesh(t,xFB)
xlabel(['forward-backward, RMSE = ', num2str(RMSE(xFB))])
subplot 224
display_3Dmesh(t,xFISTA)
xlabel(['FISTA, RMSE = ', num2str(RMSE(xFISTA))])

%%
min_crit = min([min(crit), min(critFB), min(critFISTA)]) ;
figure
subplot 121
semilogy(crit-min_crit, 'r')
hold on, semilogy(critFB-min_crit, 'k')
hold on, semilogy(critFISTA-min_crit, 'b')
xlabel('iterations $k$', 'Interpreter', 'latex'), ylabel('$f(x_k) - f(x^\star)$', 'Interpreter', 'latex')
legend('PD', 'FB', 'FISTA')
subplot 122
semilogy(cumsum(time), crit-min_crit, 'r')
hold on, semilogy(cumsum(timeFB), critFB-min_crit, 'k')
hold on, semilogy(cumsum(timeFISTA), critFISTA-min_crit, 'b')
xlim([0,1.5*sum(time)])
xlabel('time (s.)'), ylabel('$f(x_k) - f(x^\star)$', 'Interpreter', 'latex')


figure
subplot 121
plot(rmse, 'r')
hold on, plot(rmseFB, 'k')
hold on, plot(rmseFISTA, 'b')
xlabel('iterations $k$', 'Interpreter', 'latex'), ylabel('RMSE($x_k$)', 'Interpreter', 'latex')
legend('PD', 'FB', 'FISTA')
subplot 122
plot(cumsum(time), rmse, 'r')
hold on, plot(cumsum(timeFB), rmseFB, 'k')
hold on, plot(cumsum(timeFISTA), rmseFISTA, 'b')
xlim([0,1.5*sum(time)])
xlabel('time (s.)'), ylabel('RMSE($x_k$)', 'Interpreter', 'latex')




%% Constrained problem


figure
subplot 131
display_3Dmesh(t,z)
xlabel(['Noisy, RMSE = ', num2str(RMSE(z))])
subplot 132
display_3Dmesh(t,xres_c)
xlabel(['PD const., RMSE = ', num2str(RMSE(xres_c))])
subplot 133
display_3Dmesh(t,xres)
xlabel(['PD, RMSE = ', num2str(RMSE(xres))])


disp_time=2;
min_crit = min(crit_c) ;
figure
subplot 121
semilogy(crit_c-min_crit, 'r')
xlabel('iterations $k$', 'Interpreter', 'latex'), ylabel('$f(x_k) - f(x^\star)$', 'Interpreter', 'latex')
legend('PD - constrained')
subplot 122
semilogy(cumsum(time_c), crit_c-min_crit, 'r')
xlim([0,disp_time*sum(time_c)])
xlabel('time (s.)'), ylabel('$f(x_k) - f(x^\star)$', 'Interpreter', 'latex')
figure
subplot 121
plot(norm_fid_c, 'r')
hold on, plot(l2_bound*ones(size(norm_fid_c)), 'k')
xlabel('iterations $k$', 'Interpreter', 'latex'), ylabel('$\| \Phi x_k - z \|$', 'Interpreter', 'latex')
legend('PD - constrained')
subplot 122
semilogy(cumsum(time_c), norm_fid_c, 'r')
hold on, plot(cumsum(time_c), l2_bound*ones(size(norm_fid_c)), 'k')
xlim([0,disp_time*sum(time_c)])
xlabel('time (s.)'), ylabel('$\| \Phi x_k - z \|$', 'Interpreter', 'latex')


figure
subplot 121
semilogy(norm_x_c, 'r')
hold on, semilogy(norm_x, 'k')
xlabel('iterations $k$', 'Interpreter', 'latex'), ylabel('$\|x_{k+1} - x_k\|$', 'Interpreter', 'latex')
legend('PD - const', 'PD')
subplot 122
semilogy(cumsum(time_c), norm_x_c, 'r')
hold on, semilogy(cumsum(time), norm_x, 'k')
xlim([0,disp_time*sum(time_c)])
xlabel('time (s.)'), ylabel('$\|x_{k+1} - x_k\|$', 'Interpreter', 'latex')


figure
subplot 121
plot(rmse_c, 'r')
hold on, plot(rmse, 'k')
xlabel('iterations $k$', 'Interpreter', 'latex'), ylabel('RMSE($x_k$)', 'Interpreter', 'latex')
legend('PD - const', 'PD')
subplot 122
plot(cumsum(time_c), rmse_c, 'r')
hold on, plot(cumsum(time), rmse, 'k')
xlim([0,disp_time*sum(time_c)])
xlabel('time (s.)'), ylabel('RMSE($x_k$)', 'Interpreter', 'latex')

