%For the default condition
%run ./boundarySampling [number of points]>output.txt
x = linspace(-10,10);
y = linspace(-10,10);
[X,Y] = meshgrid(x,y);
Z=-1*X.^2 - 16*Y.^2 + 25;
v=[0,0];
A = dlmread('output.txt');
figure()
scatter(A(:,1),A(:,2),'filled');
hold on; contour(X,Y,Z,v,'red','LineWidth',2);
title("Sampling the Boundary of a Spectrahedron");
legend('Sampled Points','Exact Boundary');