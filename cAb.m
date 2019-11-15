function alpha=cAb(c,A,b)
    x = gaussianelim(A,b);
    alpha = c*x;
end

function x = gaussianelim(A,b);
    [row,col]=size(A);
    n = row;
    x = zeros(n,1);
    for k=1:n-1
        %do partial pivoting
        for i=k+1:n 
            if A(i,k) > A(k,k)
                A([i, k], :) = A([k, i], :);
                b([i, k], :) = b([k, i], :);
            end
        end
        for i=k+1:n
            xMultiplier = A(i,k)/A(k,k);
            for j=k+1:n
                A(i,j) = A(i,j)-xMultiplier*A(k,j);
            end
            b(i) = b(i)-xMultiplier*b(k);
        end
    end
    % backsubstitution:
    x(n) = b(n)/A(n,n);
    for i=n-1:-1:1
        summation = b(i);
        for j=i+1:n
            summation = summation-A(i,j)*x(j);
        end
        x(i) = summation/A(i,i);
    end
end