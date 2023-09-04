function xnew = jacobi(A,b)

global eps Nmax


n = length(b);

xnew = zeros(n,1);
xold = zeros(n,1);
count = 0;

while 1
    
    for i = 1:n
       sum = 0; 
       for j = 1:n
           if j ~= i
               sum =  sum + A(i,j) * xold(j);
           end
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