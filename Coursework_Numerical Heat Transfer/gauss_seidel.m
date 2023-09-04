function xnew = gauss_seidel(A,b)

global eps Nmax


n = length(b);

xnew = zeros(n,1);
xold = zeros(n,1);
count = 0;

while 1
    
    for i = 1:n
       sum = 0; 
       for j = 1:n
           xj = 0; % by default, i.e., i == j, xj = 0
           if j < i
               xj = xnew(j);
           elseif j > i
               xj = xold(j);
%            else
%                xj = 0;
           end
           sum =  sum + A(i,j) * xj;
       end
       xnew(i) = (-sum + b(i)) / A(i,i);
    end
    
    if max(abs(xnew - xold)) < eps
        
%            x = A * xnew;
        disp('converged!');
        disp(xnew);
        break;
    else
       disp(count);
%      disp(xnew);
    end
       
    if(count == Nmax)
        
        disp('max iteration cyles exceeded!');
        disp(xnew);
        break;
    end
    count = count + 1;
    xold = xnew;
end

end