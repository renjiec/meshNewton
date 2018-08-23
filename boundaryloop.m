function [ ind ] = boundaryloop( t )
  i = reshape(t, numel(t),1);
  j = reshape(t(:,[2 3 1]),numel(t),1);
  A = sparse(i,j,1);
  B=A-A';
  [x,y] = find(B==1);
  if(isempty(x))
      ind = []; return;
  end
  ind = zeros(1,1);
  
  ind(1) = y(1);
  while(true)
      n = y(find(x==ind(end)));
      if(n~=ind(1))
          ind(end+1)=n;
      else
          break;
      end
  end

end

