function [x,crit, norm_fid,time,rmse,norm_x] = Primal_Dual_constrained_graph...
    (x0, z, l2_bound, Ind_current,Ind_neighb,eta, L, Lt, Neighb_mat, xmin,xmax,normL2, NbIt,tri, RMSE, Stop_norm, Stop_crit, display, display_graph)

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

disp('******************************************')
disp('Primal Dual Algorithm for 3D mesh denoising')
disp('Constrained problem')
disp('******************************************')


% --------------------------------------------------------------
% dimensions
nb_nodes = length(z) ;
nb_edges = length(Ind_current)/2 ;
% Initialisation
x = x0 ;
u = zeros(2*nb_edges,3);
Lx = L(x) ;
v = zeros(size(x));
% --------------------------------------------------------------

% Define objective function
f =@(x) Obj_function( x,Neighb_mat,Ind_current,Ind_neighb, eta ) ;
norm_res =@(x) Norm_fidelity(x,z) ; % norm of the residual for constrained problem

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TO COMPLETE
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define projection on l2 ball, center 0, radius l2_bound
proj_l2 =@(x) ... ;
    % Define projection on l2 ball, center z, radius l2_bound
proj_l2_c =@(x) ... ;
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define proximity operator for TV regularisation
proxg =@(x,T) prox_tv_graph( x , Neighb_mat , Ind_current , T ) ;  % prox TV



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TO COMPLETE
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step-sizes
gamma = 20 ; % DO NOT CHANGE - optimised for bunny example
sigma = ... ;
    tau = ... ;
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('******************************************')
disp(['step size primal: tau = ',num2str(tau)])
disp(['step size dual: sigma = ',num2str(sigma)])
disp('******************************************')


crit = zeros(NbIt+1, 1) ;
crit(1) = f(x) ;
norm_fid = zeros(NbIt+1, 1) ;
norm_fid(1) = norm_res(x) ;
norm_x = zeros(NbIt, 1) ;
crit_it = zeros(NbIt, 1) ;
norm_x(1) = norm(x(:)-z(:))/norm(x(:)) ;
rmse = zeros(NbIt+1, 1) ;
rmse(1) = RMSE(x) ;
time = zeros(NbIt,1) ;

display_values(0 , rmse(1), crit(1), norm_fid(1), l2_bound,norm_x(1), 0 )

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
    ...
    % --------------------------------------------------------------
    
    % --------------------------------------------------------------
    % dual update : prox TV
    % --------------------------------------------------------------
    Lx = ... ;
    ...
    % --------------------------------------------------------------

    % --------------------------------------------------------------
    % dual update : data fidelity
    % --------------------------------------------------------------
    ...
    % --------------------------------------------------------------
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % --------------------------------------------------------------
    % Save information
    time(it+1) = toc(start_t);
    
    crit(it+1) = f( x ) ;
    norm_fid(it+1) = norm_res(x) ;
    crit_it(it) = abs(crit(it+1)-crit(it))/abs(crit(it+1)) ;
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
        display_values(it , rmse(it+1),crit(it+1), norm_fid(it+1), l2_bound,norm_x(it), crit_it(it) ) ;
        pause(1)
    end
    % --------------------------------------------------------------
    
    % --------------------------------------------------------------
    % stopping criteria
    if it >3 ...
            && norm_x(it) < Stop_norm ...
            && crit_it(it) < Stop_crit ...
            && norm_fid(it) < 1.1 * l2_bound
        break;
    end
    % --------------------------------------------------------------
    
end

crit = crit(1:it+1) ;
norm_fid = norm_fid(1:it+1) ;
norm_x = norm_x(1:it+1) ;
rmse = rmse(1:it+1) ;
time = time(1:it+1) ;

% --------------------------------------------------------------
% Final display
disp('-------------------------------------------')
disp(['FINAL RMSE : ',num2str(rmse(it+1))])
disp(['Minimum value : ',num2str(crit(it+1))])
disp(['TOTAL TIME : ',num2str(sum(time))]);
disp(['TOTAL NUMBER OF ITERATIONS : ',num2str(it)]);
disp('-------------------------------------------')
% --------------------------------------------------------------


end


% *************************************************************************
% FONCTIONS
% *************************************************************************

function display_values(n , RMSE,crit, norm_fid, l2_bound, normx, crit_it )
disp('-------------------------------------------')
disp(['Iteration : ',num2str(n)])
disp('-------------------------------------------')
disp(['Crit = ',num2str(crit)])
disp(['Norm residual = ',num2str(norm_fid)])
disp([' vs. l2 bound = ', num2str(l2_bound)])
disp(['norm(crit - critold) = ',num2str(crit_it)])
disp(['norm(x-xold) = ',num2str(normx)])
disp('-------------------------------------------')
disp(['RMSE     : ',num2str(RMSE)])
disp('-------------------------------------------')
disp(' ')
end


function crit = Obj_function( x,Neighb_mat,Ind_current,Ind_neighb, eta )

u = x(Ind_current,:) - x(Ind_neighb,:) ;
s = sqrt(Neighb_mat*(u.^2));
tv = eta*sum(s(:));

crit = tv ;

end


function norm_fid = Norm_fidelity( x,z )

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TO COMPLETE
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
norm_fid = ...  ; % || x - z ||
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end




