

function [uv] = floater97 (X,T)
L = adjacent(T);
L=L|L';
bounds = boundaryloop(T);
assert(length(bounds)>5,'Small Boundary');
t = linspace(0,2*pi,length(bounds)+1);
Bv = [cos(t') sin(t')];
Bv(end,:)=[];

L=diag(sparse(sum(L,2)))-L;
rhs = zeros(size(X,1),2);
rhs(bounds,:)=Bv;
L(bounds,:)=0;
L(bounds,bounds)=eye(numel(bounds));
V = (L\rhs)';
uv = V'*sqrt(1/pi);
end
