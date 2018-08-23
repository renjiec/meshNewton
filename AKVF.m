function [val, Ki, Kj, K ] = AKVF( X,T,A2, imX)
 nV = size(X,1);
 nt = size(T,1);

u = T(:,1);
v = T(:,2);
w = T(:,3);
x1 = real(imX(u)); x2 = real(imX(v)); x3 = real(imX(w));
y1 = imag(imX(u)); y2 = imag(imX(v)); y3 = imag(imX(w));
x12 = x1 - x2;
y12 = y1 - y2;
x23 = x2 - x3;
y23 = y2 - y3;
x31 = x3 - x1;
y31 = y3 - y1;

area = abs ( y12.*x31 - x12.*y31 );
area = (area.* area)./(A2/2);
area = sqrt(0.5 * area);
x12 = x12./area;
y12 = y12./area;
x23 = x23./area;
y23 = y23./area;
x31 = x31./area;
y31 = y31./area;

val = [x23.*x23 + 2*y23.*y23; 2*x23.*x23 + y23.*y23; -x23.*y23;...
       x31.*x31 + 2*y31.*y31; 2*x31.*x31 + y31.*y31; -x31.*y31;...
       x12.*x12 + 2*y12.*y12; 2*x12.*x12 + y12.*y12; -x12.*y12;...
       x31.*x23 + 2*y23.*y31; 2*x31.*x23 + y23.*y31; -x31.*y23; -x23.*y31;...
       x12.*x23 + 2*y12.*y23; 2*x12.*x23 + y12.*y23; -x12.*y23; -x23.*y12;...
       x12.*x31 + 2*y12.*y31; 2*x12.*x31 + y12.*y31; -x12.*y31; -x31.*y12;...
       -x23.*y23;...
       -x31.*y31;...
       -x12.*y12;...
       x31.*x23 + 2*y23.*y31; 2*x31.*x23 + y23.*y31; -x31.*y23; -x23.*y31;...
       x12.*x23 + 2*y12.*y23; 2*x12.*x23 + y12.*y23; -x12.*y23; -x23.*y12;...
       x12.*x31 + 2*y12.*y31; 2*x12.*x31 + y12.*y31; -x12.*y31; -x31.*y12];

    Ki = [u;u+nV;u+nV; v;v+nV;v+nV; w;w+nV;w+nV; u; u+nV; u+nV; u;u; u+nV; u+nV; u;v; v+nV; v+nV; v;    u; v;w; v; v+nV; v; v+nV;w; w+nV; w; w+nV;w; w+nV; w; w+nV];
    Kj = [u;u+nV;u; v;v+nV;v; w;w+nV;w; v; v+nV; v; v+nV;w; w+nV; w; w+nV;w; w+nV; w; w+nV;    u+nV; v+nV; w+nV; u; u+nV; u+nV; u;u; u+nV; u+nV; u;v; v+nV; v+nV; v];
    if(nargout>1)
      K = sparse(Ki, Kj, val);
    end

end

