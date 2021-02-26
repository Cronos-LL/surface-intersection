clear all 
close all
clc

%% données de départ, centre de l'ellipsoïde, matrice de covariance, et rotation pour aligner sur l'axe y

c = [-15944.907660346551, -16133.554742185399, 2642181.8304501893];

span = 6000;

angle_sat = atan2(c(2), c(1));
rotation = rot_z(angle_sat);

c = (rotation*c')';

cov = [762292.80181018, 324361.96788065, -206982.2744863; 324361.96788065, 261521.23380526, -118602.16131351;...
    -206982.2744863, -118602.16131351, 276185.96438456];

cov = rotation*cov*rotation';
cov = inv(cov);

%% calcul de f1 et f2, fonctions représentants les surfaces du cône et de l'ellispoïde

x = linspace(c(1) - 6000, c(1) + 6000, 50); 
y = linspace(c(2) - 6000, c(2) + 6000, 50); 
z = linspace(c(3) - 6000, c(3) + 6000, 50); 

[x3, y3, z3] = meshgrid(x, y, z);

f1 = (x3-c(1)).^2*cov(1,1) + (y3-c(2)).^2*cov(2,2) + (z3-c(3)).^2*cov(3,3) + 2*(x3-c(1)).*(y3-c(2))*cov(1,2) ...
    + 2*(x3-c(1)).*(z3-c(3))*cov(1,3) + 2*(y3-c(2)).*(z3-c(3))*cov(2,3) - 6^2;

k = tan(0.4*pi/180);

f2 = y3.^2 + x3.^2 - z3.^2*k^2;

%% détermination de l'intersection entre l'ellispoïde et le cône en utilisant la fonction contour

tic

n1 = 200;
n2 = 20;

phi = linspace(-pi/2, pi/2, n1);

x2 = linspace(c(3)-span, c(3)+span, n2)*k.*cos(phi)';
y2 = linspace(c(3)-span, c(3)+span, n2)*k.*sin(phi)';

z2 = sqrt(x2.^2 + y2.^2)/k;


f3 = (x2 - c(1)).^2*cov(1,1) + (y2 - c(2)).^2*cov(2,2) + (z2 - c(3)).^2*cov(3,3)...
    + 2*(x2 - c(1)).*(y2 - c(2))*cov(1,2) + 2*(x2 - c(1)).*(z2 - c(3))*cov(1,3)...
    + 2*(y2 - c(2)).*(z2 - c(3))*cov(2,3) - 6^2;

C = contour(x2, y2, f3, [0 0]);

xL = C(1, 2:end);
yL = C(2, 2:end);

zL = sqrt(xL.^2 + yL.^2)/k;

toc

%% plot en 3d des deux surfaces et de la courbe noire représentant l'intersection

figure(2);
hold on
plot3(xL, yL, zL, 'k-', 'Linewidth', 1);

patch(isosurface(x3, y3, z3, f1, 0), 'FaceColor', [0.5 1.0 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
patch(isosurface(x3, y3, z3, f2, 0), 'FaceColor', [1.0 0.5 0.0], 'EdgeColor', 'none', 'FaceAlpha', 0.6);

view(3); camlight; axis vis3d;
xlabel('Axe x')
ylabel('Axe y')
zlabel('Axe z')
legend('Intersection', 'Ellipsoïde', 'Cône')
title("Intersection obtenue entre l'ellipsoïde associé à l'incertitude en position du satellite et le cône représentatif du lobe principal de l'antenne")
hold off

%% fonctions

function matrice = rot_z(angle)

    matrice = [[cos(angle), sin(angle), 0]; [-sin(angle), cos(angle), 0]; [0, 0, 1]];
    
end






