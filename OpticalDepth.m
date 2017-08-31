function f = OpticalDepth(alpham,rangebin)
% OpticalDepth(alpham) integrates absorption coefficient alpham (2D matrix)
% from initial height to height r, such that f is a function of height too.

numrow = size(alpham,1);            % Number of rows in alpha matrix
optdepth = [];

for i=1:numrow
    if i==1
        optdepth(i,:) = alpham(1,:)*rangebin;
    else 
        alphai = alpham(1:i,:)*rangebin;     % Selects range of rows from from row 1 to row i
        optdepth(i,:) = trapz(alphai);
    end
end
f = optdepth;
end