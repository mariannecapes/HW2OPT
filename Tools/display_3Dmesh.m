function display_3Dmesh( tri, x )
trisurf(tri,x(:,1),x(:,2),x(:,3))
axis equal, %axis equal sets the aspect ratio so that the data units are the same in every direction. The aspect ratio of the x-, y-, and z-axis is adjusted automatically according to the range of data units in the x, y, and z directions.
axis tight, %axis tight sets the axis limits to the range of the data
view(10,80) 
pause(0.1)
end
