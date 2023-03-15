

name_data = 'dragon3D_true' ;
load(['Results/', name_data, '_res.mat'])

l2_err =@(x) [ sum((xbar(:,1)-x(:,1)).^2)  sum((xbar(:,2)-x(:,2)).^2)  sum((xbar(:,3)-x(:,3)).^2) ] ;
RMSE =@(x) sqrt(sum(l2_err(x)));

%% Plot original graph and noisy graph

figure
subplot 121
display_3Dmesh(t,xbar)
xlabel('Initial graph')
subplot 122
display_3Dmesh(t,z)
xlabel('Noisy graph')

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
disp_time=2;
min_crit = 0.99*min([min(crit), min(critFB), min(critFISTA)]) ;
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
xlim([0,disp_time*sum(time)])
xlabel('time (s.)'), ylabel('$f(x_k) - f(x^\star)$', 'Interpreter', 'latex')


figure
subplot 121
semilogy(norm_x(2:end), 'r')
hold on, semilogy(norm_xFB(2:end), 'k')
hold on, semilogy(norm_xFISTA(2:end), 'b')
xlabel('iterations $k$', 'Interpreter', 'latex'), ylabel('$\|x_{k+1} - x_k\|$', 'Interpreter', 'latex')
legend('PD', 'FB', 'FISTA')
subplot 122
semilogy(cumsum(time), norm_x, 'r')
hold on, semilogy(cumsum(timeFB), norm_xFB, 'k')
hold on, semilogy(cumsum(timeFISTA), norm_xFISTA, 'b')
xlim([0,disp_time*sum(time)])
xlabel('time (s.)'), ylabel('$\|x_{k+1} - x_k\|$', 'Interpreter', 'latex')


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
xlim([0,disp_time*sum(time)])
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

