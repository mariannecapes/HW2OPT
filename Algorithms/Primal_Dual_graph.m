function [x,crit,rmse,time, norm_x] = Primal_Dual_graph...
    (x0, z, Ind_current,eta, L, Lt, Neighb_mat, xmin,xmax,normL2, NbIt,tri, RMSE, Stop_norm, Stop_crit, display, display_graph)

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

disp('******************************************')
disp('Primal Dual Algorithm for 3D mesh denoising')
disp('******************************************')


% --------------------------------------------------------------
% dimensions
nb_nodes = length(z) ;
nb_edges = length(Ind_current)/2 ;
% Initialisation 
x = x0 ;
u = zeros(2*nb_edges,3);
Lx = L(x) ;
% --------------------------------------------------------------

% Define objective function
f =@(x) Obj_function( x,z,Neighb_mat, L, eta ) ;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TO COMPLETE
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define gradient operator
gradh =@(x) norm(x-z);
gradh = real(x -z); %TBC
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define proximity operator for TV regularisation
proxg =@(x,T) prox_tv_graph( x , Neighb_mat , Ind_current , T ) ;  % prox TV


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TO COMPLETE
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step-sizes
gamma = 10 ;
sigma = gamma/normL2 ; % DO NOT CHANGE - optimised for bunny example
tau = ... ;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('******************************************')
disp(['step size primal: tau = ',num2str(tau)])
disp(['step size dual: sigma = ',num2str(sigma)])
disp('******************************************')


% Define variables
crit = zeros(NbIt+1, 1) ;
crit(1) = f(x) ;
norm_x = zeros(NbIt, 1) ;
norm_x(1) = norm(x(:)-z(:))/norm(x(:)) ;
rmse = zeros(NbIt+1, 1) ;
rmse(1) = RMSE(x) ;
time = zeros(NbIt,1) ;

display_values(0 , rmse(1), crit(1),norm_x(1) )

if(display_graph>0)
figure(100)
subplot 121
display_3Dmesh(tri,z)
title('3D mesh with noise added to it')
subplot 122
display_3Dmesh(tri,x)
title(['Iteration ',num2str(0)])
end



for it = 1:NbIt
    
    start_t = tic;
    xold = x ;
    Lxold = Lx ;
    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TO COMPLETE
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --------------------------------------------------------------
    % primal update: gradients and 3D box constraint
    % --------------------------------------------------------------
    x_int = xold - tau * (gradh(xold) + Lt(u));
    x = proj_box(x_int, xmin, xmax);
    % --------------------------------------------------------------
    
    % --------------------------------------------------------------
    % dual update : prox TV
    % --------------------------------------------------------------
    Lx = 2*Lx - Lxold;
    u_int = u+sigma*Lx;
    u = u_int - sigma*proxg(u_int, 1/sigma);
   
    ...
    % --------------------------------------------------------------
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % --------------------------------------------------------------
    % Save information
    time(it+1) = toc(start_t);
    
    crit(it+1) = f( x ) ;
    norm_x(it+1) = norm(x(:)-xold(:))/norm(x(:)) ;
    rmse(it+1) = RMSE(x);
    % --------------------------------------------------------------
    
    % --------------------------------------------------------------
    % Display if needed
    if(mod(it,display_graph)==0)
        figure(100)
        subplot 133
        display_3Dmesh(tri,x)
        title(['Iteration ',num2str(it)])
    end
    if(mod(it,display)==0)
        display_values(it , rmse(it+1),crit(it+1),norm_x(it) ) ;
        pause(1)
    end
    % --------------------------------------------------------------
    
    % --------------------------------------------------------------
    % stopping criteria
    if it >3 ...
            && norm_x(it) < Stop_norm ...
            && abs(crit(it+1) - crit(it)) < Stop_crit * crit(it+1)
        break;
    end
    % --------------------------------------------------------------
    
end

crit = crit(1:it+1) ;
norm_x = norm_x(1:it+1) ;
rmse = rmse(1:it+1) ;
time = time(1:it+1) ;

% --------------------------------------------------------------
% Final display
disp('-------------------------------------------')
disp(['TOTAL NUMBER OF ITERATIONS : ',num2str(it)]);
disp(['TOTAL TIME : ',num2str(sum(time))]);
disp(['FINAL RMSE : ',num2str(rmse(it+1))])
disp(['Minimum value : ',num2str(crit(it+1))])
disp('-------------------------------------------')
% --------------------------------------------------------------


end


% *************************************************************************
% FONCTIONS
% *************************************************************************

function display_values(n , RMSE,crit, normx )
disp('-------------------------------------------')
disp(['Iteration : ',num2str(n)])
disp('-------------------------------------------')
disp(['Crit = ',num2str(crit)])
disp(['norm(x-xold) = ',num2str(normx)])
disp('-------------------------------------------')
disp(['RMSE     : ',num2str(RMSE)])
disp('-------------------------------------------')
disp(' ')
end

 
function crit = Obj_function( x,z,Neighb_mat,L, eta )

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TO COMPLETE
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid =   ;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lx = L(x) ;
s = sqrt(Neighb_mat*(Lx.^2));
tv = eta*sum(s(:));

crit = fid + tv ;

end

function p  = prox_tv_graph(u,Neighb_mat,Ind_current,T)

p = zeros(size(u)) ;
s = sqrt(Neighb_mat*u.^2) ; % l2-norm per node n: || (u_i)_{i in V_n} ||
s = s(Ind_current,:) ;

for i = 1:3
    ind = (s(:,i) > T) ;
    p(ind, i) = (1-T./ s(ind, i)) .* u(ind, i) ;
end

end



