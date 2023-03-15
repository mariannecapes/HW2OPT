function p = prox_TV_graph_ind(z, xmin, xmax, Neighb_mat,Ind_current,L, Lt,T, eta, normL2)
% Compute the proximity operator of TV + ind on [xmin, xmax] for a 3D mesh

gamma = 1.9 * min( 1, 1/normL2) ;
nb_edges = length(Ind_current)/2 ;
nb_nodes = length(z) ;

u1 = zeros(2*nb_edges,3);
u2 = zeros(nb_nodes,3);

NbIt = 300 ;
stop_norm = 1e-4 ;
stop_crit = 1e-4 ;

x = z ;

for it = 1:NbIt
    xold = x;

    x = z - 0.5 * ( Lt(u1) + u2 ) ;
    
    Lx = L(x) ;
    u1_ = u1 + gamma * Lx ;
    p1 = prox_tv_graph(u1_/gamma,Neighb_mat,Ind_current,2*T*eta/gamma) ;
    u1 = u1_ - gamma * p1 ;
    
    u2_ = u2 + gamma * x ;
    p2 = proj_box(u2_/gamma,xmin,xmax) ;
    u2 = u2_ - gamma * p2 ;
    
    crit(it) = comp_crit(x,z,Neighb_mat,L, T, eta ) ;
    norm_x(it) = norm(x(:)-xold(:)) ;
    
    
    if mod(it,500)==0
        disp('*********************')
        disp(['it = ', num2str(it)])
        disp(['crit = ', num2str(crit(it))])
        disp(['norm(x-xold) = ', num2str(norm_x(it))])
        figure(10)
        subplot 121, semilogy(crit), ylabel('crit')
        subplot 122, semilogy(norm_x), ylabel('norm( x-xold)')
        pause(0.1)
    end
    
    if it>10 && ...
            abs(crit(it)-crit(it-1)) < stop_crit * crit(it) &&  ...
            norm_x(it) < stop_norm * norm(x(:))
        break
    end
    
end
p=x ;

end

function crit = comp_crit(x,z,Neighb_mat,L, T, eta )

u = L(x) ;
s = sqrt(Neighb_mat*(u.^2));
tv = T*eta*sum(s(:));

crit = 0.5* sum(abs(x(:)-z(:)).^2) + tv ;

end

function p  = prox_tv_graph(u,Neighb_mat,Ind_current,T)

p = zeros(size(u)) ;
s = sqrt(Neighb_mat*u.^2) ;
s = s(Ind_current,:) ;

for i = 1:3
    ind = (s(:,i) > T) ;
    p(ind, i) = (1-T./ s(ind, i)) .* u(ind, i) ;
end

end


