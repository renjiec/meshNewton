function [z, allStats,preData] = deformMeshNew(x, t, P2PID, P2PCP, z, nIter, lambda, energy_type, e_param, methodName, preData, trimeshH)
  fC2R = @(x) [real(x) imag(x)];
  fR2C = @(x) complex(x(:,1), x(:,2));
  Gthresh = 1e-3;
  Ethresh = 1e-17;
  isDeform = ~isempty(P2PID);
  if strcmp(energy_type,'SARAP')
    Ethresh = 1e-8;
  end
  if isDeform
    Gthresh = 1e-10;
    Ethresh = 1e-6;
  end
  nf = size(t, 1);
  nv = size(x, 1);
if strcmp(methodName,'CM')
 assert( ~(strcmp(energy_type,'BCONF') || strcmp(energy_type,'BARAP')),'Unsupported energy'); 
end

  faceElens = sqrt( meshFaceEdgeLen2s(x, t) );
  faceAngles = meshAnglesFromFaceEdgeLen2(faceElens.^2);
  flatXPerFace = [zeros(nf,1) faceElens(:,3) faceElens(:,2).*exp(1i*faceAngles(:,1))];


  %% initilization
  Areas = signedAreas(x, t);
  Tarea = sum(Areas);

  isometric_energyies = { 'ISO', 'ExpSD', 'AMIPS', 'SARAP', 'HOOK', 'ARAP', 'BARAP', 'BCONF'};
  hessian_projections = { 'NP', 'KP-Newton(ours)', 'FP-Newton', 'FP6', 'CM' };


  energy_param = e_param;
  hession_proj = methodName;

  mexoption = struct('energy_type', find(strcmpi(energy_type, isometric_energyies), 1) - 1, ...
                     'hessian_projection', find(strcmpi(hession_proj, hessian_projections), 1) - 1, ...
                     'energy_param', energy_param, 'verbose', 0);

  D2 = -1i/4*(flatXPerFace(:,[2 3 1])-flatXPerFace(:,[3 1 2]))./Areas;
  D = sparse(repmat(1:nf,3,1)', t, D2);
  D2t = D2.';
  DirDiff = (1i*([conj(flatXPerFace(:,3)) flatXPerFace(:,1) flatXPerFace(:,2)] - [flatXPerFace(:,2) conj(flatXPerFace(:,3)) flatXPerFace(:,1)])./repmat(Areas*2,1,3)).';
  fDeformEnergy = @(z) meshIsometricEnergyC(D*conj(z), D*z, D2t, Areas, mexoption); 
  if isDeform
    fP2PEnergy = @(z) norm(z(P2PID)-P2PCP)^2 * lambda;
  else
    fP2PEnergy = @(z) 0;
  end

  Iscale = 1e-10;
  t2 = [t t+nv]';
  %% initialization, get sparse matrix pattern
  if strcmp(methodName,'AKVF')
    [~, Mi, Mj, ~] = AKVF( x,t,Areas, z);
    Tt = uint64((t).');
  else
    [xmesh, ymesh] = meshgrid(1:6, 1:6);
    Mi = t2(xmesh(:), :);
    Mj = t2(ymesh(:), :);
  end
    H = sparse(Mi,Mj, Mi);
    [Hi, Hj] = find(H);
    Hi = uint64(Hi);
    Hj = uint64(Hj);
    nzidx = ij2nzIdxs(H, uint64(Mi), uint64(Mj));
    % nonzero indices of the matrix
    Hnonzeros0 = zeros(nnz(H),1);
    idxDiagH = ij2nzIdxs(H, uint64(1:nv*2), uint64(1:nv*2));
    Hnonzeros0(idxDiagH) = Iscale*2;
  scale = 1;
  if(~isDeform)
    x = abs(D*conj(z)).^2;
    y = abs(D*z).^2;
    scales = [];
    if strcmp('ISO',energy_type) &&  strcmp('KP-Newton(ours)',methodName)
      scales = sort((((x-y).^-3).*(x+3*y)).^(1/4));
    end
    if strcmp('BCONF',energy_type) &&  strcmp('KP-Newton(ours)',methodName)
      scales = sort((x-y).^-0.5);
    end
    if strcmp('BARAP',energy_type) &&  strcmp('KP-Newton(ours)',methodName)
      scales = zeros(nf,1);
      for f = 1:nf
       mscale = eig([0 0 0 (e_param*(x(f)-y(f))^-2)/(2+e_param);1 0 0 0;0 1 0 0; 0 0 1 2/(2+e_param) * sqrt(x(f))^-1]);
       scales(f) = max(mscale.*(abs(imag(mscale))<1e-30));
      end
      scales = sort(scales);
    end
    if ~isempty(scales)
     b = [1 -1]; diff = filter(b,1,scales);
     scale = scales(find((diff(ceil(nf/2):end)>0.1 ),1,'first')+floor(nf/2)-1);
     if isempty(scale)
       scale = max(scales);
     end
     z = z* scale;
    end
  end


  %% main loop
  g2GIdx = uint64(t2);
  allStats = zeros(nIter+1, 8); % statistics
  en = fDeformEnergy(z) +fP2PEnergy(z);
  allStats(1, [5 7 8]) = [0 0  en/Tarea];
  for it=1:nIter
    tt = tic;
    NE = 0;
    fz = D*conj(z);
    gz = D*z;
    if strcmp(methodName,'AKVF')
      [e, g] = meshIsometricEnergyC(fz, gz, D2t, Areas, mexoption);
      hs = AKVFC(Areas*2,z.',Tt);
    else
      [e, g, hs,NE] = meshIsometricEnergyC(fz, gz, D2t, Areas, mexoption);
    end
    G = myaccumarray(g2GIdx, g, zeros(nv*2,1));
    Hnonzeros = myaccumarray( nzidx, hs, Hnonzeros0 );
    if isDeform
      G([P2PID P2PID+nv]) = G([P2PID P2PID+nv]) + (fC2R(z(P2PID))-fC2R(P2PCP))*lambda;
      Hnonzeros(idxDiagH([P2PID P2PID+nv])) =  Hnonzeros(idxDiagH([P2PID P2PID+nv])) + lambda;
    end

    s = 1./sqrt(Hnonzeros(idxDiagH));
    Hnonzeros2 = jacobiC(Hnonzeros, s, Hi, Hj);
    Hmat = sparse(double(Hi),double(Hj),Hnonzeros2+Hnonzeros0);
    [dz] = Hmat\(-s.*G);
    dz = s.*dz;

      
    
    dz = fR2C( reshape(dz, [], 2) );

    %% orientation preservation
    dfz = D*conj(dz);
    dgz = D*dz;
    ls_t = min( maxtForPositiveArea(fz, gz, dfz, dgz)*0.9, 1 );

    %% line search energy decreasing
    fMyFun = @(t) fDeformEnergy( dz*t + z ) + fP2PEnergy( dz*t + z );
    normdz = norm(dz);

    dgdotfz = dot( G, [real(dz); imag(dz)] );
    
    ls_alpha = 0.2; ls_beta = 0.5;
    fQPEstim = @(t) en+ls_alpha*t*dgdotfz;

    e_new = fMyFun(ls_t);
    while ls_t*normdz>1e-12 && e_new > fQPEstim(ls_t)
        ls_t = ls_t*ls_beta;
        e_new = fMyFun(ls_t);
    end
    ediff = en-e_new;
    en = e_new;
    
    %% update
    z = dz*ls_t + z;
    allStats(it+1, [3 4 5 7 8]) = [NE norm(G) (toc(tt)*1000) 0  en/Tarea];
    if (norm(G) <=Gthresh) || (ediff/Tarea<Ethresh)
      break; 
    end
    if(nargin>=12)
      title(strcat(methodName,' Iteration: ', num2str(it),' Energy: ', num2str(en,'%10.5e\n' )));
      set(trimeshH,'vertices',[real(z) imag(z) zeros(size(z))]);drawnow; %2D PARAM
  end
  end


end