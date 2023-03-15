function [x, crit, rmse, time, norm_x] = FB...
    (x0, z, Ind_current,eta, L, Lt, Neighb_mat, xmin,xmax,normL2, NbIt,tri, RMSE, Stop_norm, Stop_crit, display, display_graph)


% --------------------------------------------------------------
% Initialisation 
x = x0 ;
% --------------------------------------------------------------

% Define objective function
f =@(x) Obj_function( x,z,Neighb_mat,L, eta ) ;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TO COMPLETE
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define step-size
gamma = ... ;

% Define gradient operator
gradh =@(x) ... ;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define proximity operator
prox =@(x, T) prox_TV_graph_ind(x, xmin, xmax, Neighb_mat,Ind_current,L, Lt,T, eta, normL2) ;

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


% Iterations
for it = 1:NbIt
    start = tic;
    xold = x ;
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TO COMPLETE
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % iterates
    ...
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    time(it+1) = toc(start) ;
    
    crit(it+1) = f(x) ;
    norm_x(it+1) = norm(x(:)-xold(:))/norm(x(:)) ;
    rmse(it+1) = RMSE(x) ;
    
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
    
    if it >3 ...
            && abs(crit(it+1)-crit(it))/abs(crit(it+1)) < Stop_crit ...
            && norm_x(it) < Stop_norm
        break
    end
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
fid = ...  ;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lx = L(x) ;
s = sqrt(Neighb_mat*(Lx.^2));
tv = eta*sum(s(:));

crit = fid + tv ;

end
