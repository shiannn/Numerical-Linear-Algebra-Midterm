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