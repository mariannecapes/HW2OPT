function p = proj_box(x,xmin,xmax)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TO COMPLETE
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p(:,1) = min(max(xmin, x(:,1)), xmax); % projection dim 1
p(:,2) = min(max(xmin, x(:,2)), xmax); % projection dim 2
p(:,3) = min(max(xmin, x(:,2)), xmax); % projection dim 3
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end