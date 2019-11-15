function x = gaussianelim(A,B);
    [row,col]=size(A);
    [rowB,colB]=size(B);
    n = row;
    x = zeros(n,colB);
    for k=1:n-1
        %do partial pivoting
        for i=k+1:n 
            if A(i,k) > A(k,k)
                A([i, k], :) = A([k, i], :);
                B([i, k], :) = B([k, i], :);
            end
        end
        for i=k+1:n
            xMultiplier = A(i,k)/A(k,k);
            for j=k+1:n
                A(i,j) = A(i,j)-xMultiplier*A(k,j);
            end
            B(i,:) = B(i,:)-xMultiplier*B(k,:);
        end
    end
    % backsubstitution:
    for k=1:colB
        x(n,k) = B(n,k)/A(n,n);
        for i=n-1:-1:1
            summation = B(i,k);
            for j=i+1:n
                summation = summation-A(i,j)*x(j,k);
            end
            x(i,k) = summation/A(i,i);
        end
    end
end