function [qpts] = Quadrature_2D_Quadilateral_element(n)
if nargin==0
    n=3;
end
q=quadrature(n);
qpts=zeros(3,n^2);
for i=1:n
        qpts(1,1+n*(i-1):n+n*(i-1))=q(1,i);
        qpts(2:3,1+n*(i-1):n+n*(i-1))=q(1:2,1:end);
        qpts(3,1+n*(i-1):n+n*(i-1))=qpts(3,1+n*(i-1):n+n*(i-1)).*q(2,i);
        
end
end
function [qpts] = quadrature(n)
%   QUADRATURE
%     quadrature(n) returns a quadrature table for a rule with n
%     integration points.  The first row of the table gives the quadrature
%     point location, and the second gives the quadrature weights.

    u = 1:n-1;
    u = u./sqrt(4*u.^2 - 1);

    A = zeros(n);
    A(2:n+1:n*(n-1)) = u;
    A(n+1:n+1:n^2-1) = u;

    [v, x] = eig(A);
    [x, k] = sort(diag(x));    
    qpts = [x'; 2*v(1,k).^2];
end

