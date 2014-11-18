figure
[x y z] = sphere(20); 
s = surface(x,y,z,'facecolor','interp','cdata',z);
set(s,'edgecolor','black','facealpha','flat','alphadata',x.*z);
alpha('scaled'); 
axis equal
box on
grid on
campos([2 13 10]); 
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Sphere with Alpha Data');
plot2svg('sphere.svg');
