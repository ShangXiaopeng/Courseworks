clear all
clc

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

Re = 10;
L = 1;
nu = 1e-6;
U = Re * nu / L;
M = 100; %
dl = L/M; % delta h
dt = 1e-5;                      
psi = zeros(M + 1,M + 1);
omega = zeros(M+1,M+1);
psi_save = zeros(25, M+1, M+1);
psi_save(1,:,:) = psi;

eps = 1e-2;
count = 1;


for k = 1:4000
    

    for inner_iter = 1:2
        
        
    err = 0;

    omega(2:M,1) = -2*(psi(2:M,2) - psi(2:M,1)) / dl^2;
    omega(2:M,M+1) = -2*(psi(2:M,M) - psi(2:M,M+1)) / dl^2;

    omega(1,2:M) = -2*(psi(2,2:M) - psi(1,2:M))/dl^2;
    omega(M+1,2:M) = -2*(psi(M,2:M) - psi(M+1,2:M) + 1 * dl)/dl^2;
   
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

    for i = 2:M
        
        for j = 2:M
        
%             psi(i,j) = 1/4 * (omega(i,j) * dl^2 + psi(i+1,j) + psi(i-1,j) + psi(i,j+1) + psi(i,j-1));

            u(i,j) = (psi(i,j+1) - psi(i,j-1)) / (2*dl);
            v(i,j) = -((psi(i+1,j) - psi(i-1,j)) / (2*dl));

            delta_omega = dt*( -1/2/dl*(u(i,j)*(omega(i+1,j) - omega(i-1,j)) + v(i,j)*(omega(i,j+1) - omega(i,j-1))) ...
            + 1/Re/ dl^2*(omega(i+1,j) + omega(i-1,j) + omega(i,j+1) + omega(i,j-1) - 4*omega(i,j)) );
        
       
            omega(i,j) = omega(i,j) + delta_omega;
           
            psi(i,j) = 1/4 * (omega(i,j) * dl^2 + psi(i+1,j) + psi(i-1,j) + psi(i,j+1) + psi(i,j-1));
            
            
            err = max(abs(delta_omega/omega(i,j)),err);
          
        end
    end
    end
    
    if (mod(k,200) == 0)
        count = count + 1;
        k   
        err
        psi_save(count,:,:) = psi; %
        pause(0.5)
    end
    
    if err < eps
        break;
    end
end

count = count + 1;
psi_save(count,:,:) = psi; %








psi_plot = ones(M+1, M+1);
for ii = 1: M+1
    for jj = 1:M+1
        psi_plot(ii,jj) = psi_save(26,ii,jj);
    end
end

figure(1),
set(gcf,'Units','centimeters','Position',[1 2 17.5 15]);
set(gca,'Position',[0.175 0.17 0.775 0.78])

contour(psi,25)
xlabel x, ylabel y, zlabel 'T (Jacob, K)'
axis([0 100.1 0 100.1])
% colorbar('SouthOutside')





        
% contour(psi)
% xlabel x, ylabel y

% contour(psi,100)
% xlabel x, ylabel y


