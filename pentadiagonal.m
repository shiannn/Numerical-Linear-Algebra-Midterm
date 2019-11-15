function x = penta(A,myalpha,f)
    %[A,myalpha,f]=getMatrix(m);
    m = size(A,1);
    
    xtemp = zeros(m+4,1);
    xtemp(1) = myalpha(1);
    xtemp(2) = myalpha(2);
    xtemp(m+3) = myalpha(3);
    xtemp(m+4) = myalpha(4);
    f = f - A*xtemp;
    
    [Arow Acol] = size(A);
    Anew = A(:,3:Acol-2);
    %x = Anew\f;
    x = gaussianelim(Anew,f);
    x = [myalpha(1:2);x;myalpha(3:4)];
end

function [A,myalpha,f]=getMatrix(m)
    coeff = rand(m,5);
    f = rand(m,1);
    myalpha = rand(4,1);
    A = zeros(m,m+4);
    for i=1:m 
        A(i,i:i+4) = coeff(i,:);
    end
    %A
    %x
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
end