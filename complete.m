function x = gaussianelimComplete(A,b);
    [row,col]=size(A);
    m = row;
    n = col;
    x = zeros(n,1);
    index = [1:n]';
    for k=1:m-1
        %do complete pivoting
        maxRow = -1;
        maxCol = -1;
        maxElement = -1;
        for i=k:m 
            for j=k:n
                if abs(A(i,j)) > maxElement
                    %row change
                    maxRow = i;
                    %column change
                    maxCol = j;
                    %Element
                    maxElement = abs(A(i,j));
                end
            end
        end
        A([maxRow,k],:) = A([k,maxRow],:);
        b([maxRow,k],:) = b([k,maxRow],:);
        A(:,[maxCol,k]) = A(:,[k,maxCol]);
        index([maxCol,k],:) = index([k,maxCol],:);
        
        for i=k+1:m
            xMultiplier = A(i,k)/A(k,k);
            A(i,k:n) = A(i,k:n) - xMultiplier*A(k,k:n);
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
    %At = A;
    %bt = b;
    toSort = [x,index];
    AfterSort = sortrows(toSort,2);
    x = AfterSort(:,1);
end