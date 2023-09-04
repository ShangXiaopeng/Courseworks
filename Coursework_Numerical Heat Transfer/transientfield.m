
clear all
close all
clc

global eps Nmax

eps = 1e-5;
Nmax = 1000;
x0 = 0;
xM = 2;
y0 = 0;
yN = 1;

T0 = 300;

deltax = 0.05;
deltay = 0.05;
delta = 0.05;
deltat = 0.01;

x = [x0:deltax:xM]'; 
y = [y0:deltay:yN]';

M = length(x) - 1; 
N = length(y) - 1;

[X,Y] = meshgrid(x,y);

% T_ana = T0 + exp(X) .* sin(pi*Y);

% internal nodes
X = X(2:N,2:M); 
Y = Y(2:N,2:M);

% source term f(x,y)
% f = (pi^2 - 1) * exp(X) .* sin(pi*Y);
f = zeros(N-1,M-1);

T = T0 * ones(N+1, M+1);
T_save = zeros(11, N+1, M+1);
T_save(1,:,:) = T;

% boundary conditions
T(:,1) = T0 + sin(pi*y);
T(:,end) = T0 + exp(2) * sin(pi * y);
T(1,:) = T0;
T(end,:) = T0;

T_jacob = T;
T_jacob0 = T_jacob;

% T_gauss_seidel = T;
% T_gauss_seidel0 = T_gauss_seidel;

% T_sor = T;
% T_sor0 = T_sor;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

p = -1 / deltax^2;
q = -1 / deltay^2;
r = 2 * (1/deltax^2 + 1/deltay^2) + 1/deltat;

% f(:,1) = f(:,1) - q * T(2:N,1);
% f(:,end) = f(:,end) - q * T(2:N,end);
% f(1,:) = f(1,:) - p * T(1,2:M);
% f(end,:) = f(end,:) - p * T(end,2:M);

count = 1;
error = 0;

for k = 1:100
    
    f = (pi^2 - 1) * exp(X) .* sin(pi*Y) + T_jacob0(2:N,2:M) / deltat;
    f(:,1) = f(:,1) - q * T_jacob(2:N,1);
    f(:,end) = f(:,end) - q * T_jacob(2:N,end);
    f(1,:) = f(1,:) - p * T_jacob(1,2:M);
    f(end,:) = f(end,:) - p * T_jacob(end,2:M);

    e = ones(M-1,1);
%   C = 1 / delta^2 * spdiags([-e 4*e -e],[-1 0 1], M-1, M-1);
    C = spdiags([p*e r*e p*e],[-1 0 1], M-1, M-1);
    D = q * eye(M-1);
    e = ones(M-1,1);
    A = kron(eye(N-1),C) + kron(spdiags([e e],[-1 1],N-1,N-1),D);

    A = full(A);
    f = f'; 
    b = f(:);

    X_jacob = jacobi(A,b);
%   X_gauss_seidel = gauss_seidel(A,b);
%   X_sor = sor(A,b);

    T_jacob(2:N,2:M) = reshape(X_jacob,M-1,N-1)';
    error = max(abs(T_jacob - T_jacob0)./T_jacob0);
    
    
    if(error < eps)
        break;
    end
    
    if mod(k,5) == 0
        count = count + 1;
        T_save(count,:,:) = T_jacob;
    end
    
    T_jacob0 = T_jacob;
end

count = count + 1;
T_save(count,:,:) = T_jacob;
 





T_plot = ones(N+1, M+1);
for ii = 1: N+1
    for jj = 1:M+1
        T_plot(ii,jj) = T_save(10,ii,jj);
    end
end

figure(1),
set(gcf,'Units','centimeters','Position',[1 2 17.5 15]);
set(gca,'Position',[0.175 0.17 0.775 0.78])

contourf(T_plot)
xlabel x, ylabel y, zlabel 'T (Jacob, K)'
colorbar('SouthOutside')




% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % post process % % % % % % % % % % % % % % % % % 
% figure(1), spy(A,5)
% figure(1), 
% set(gcf,'Units','centimeters','Position',[1 2 17.5 15]);
% set(gca,'Position',[0.175 0.17 0.775 0.78])
% 
% subplot(3,2,1), mesh(x,y,T_jacob)
% xlabel x, ylabel y, zlabel 'T (Jacob, K)'
% subplot(3,2,2), mesh(x,y,abs((T_jacob - T_ana)./T_ana*100))
% axis([0 2 0 1 0 0.006])
% xlabel x, ylabel y, zlabel 'Error (%)'
% ztick=0:0.003:0.006;
% set(gca,'ZTick',ztick,'ZTickLabel',ztick)
% 
% subplot(3,2,3), mesh(x,y,T_gauss_seidel)
% xlabel x, ylabel y, zlabel 'T(Gauss Seidel, K)'
% subplot(3,2,4), mesh(x,y,abs((T_gauss_seidel - T_ana)./T_ana*100))
% axis([0 2 0 1 0 0.006])
% xlabel x, ylabel y, zlabel 'Error (%)'
% ztick=0:0.003:0.006;
% set(gca,'ZTick',ztick,'ZTickLabel',ztick)
% 
% subplot(3,2,5), mesh(x,y,T_sor)
% xlabel x, ylabel y, zlabel 'T(SOR, K)'
% subplot(3,2,6), mesh(x,y,abs((T_sor - T_ana)./T_ana*100))
% axis([0 2 0 1 0 0.006])
% xlabel x, ylabel y, zlabel 'Error (%)'
% ztick=0:0.003:0.006;
% set(gca,'ZTick',ztick,'ZTickLabel',ztick)

% temperature field
% figure(2),
% set(gcf,'Units','centimeters','Position',[1 2 17.5 15]);
% set(gca,'Position',[0.175 0.17 0.775 0.78])
% 
% subplot(2,2,1), contourf(T_ana)
% xlabel x, ylabel y, zlabel 'T (Analytical, K)'
% 
% subplot(2,2,2), contourf(T_jacob)
% xlabel x, ylabel y, zlabel 'T (Jacob, K)'
% 
% subplot(2,2,3), contourf(T_gauss_seidel)
% xlabel x, ylabel y, zlabel 'T (Gauss Seidel, K)'
% 
% subplot(2,2,4), contourf(T_sor)
% xlabel x, ylabel y, zlabel 'T (SOR, K)'
% 
% colorbar('SouthOutside')
% % colorbar('Ticks',[300,302,304,306,308],...
% %          'TickLabels',{'300','302','304','306','308'})


