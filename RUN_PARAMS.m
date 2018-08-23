
%input obj file path
inDir = 'venus_parametrized_nonlin.obj';
iEnergy = 1;
energies = {'ISO','SARAP','BCONF','BARAP'};
[X,T,~,~] = readObjCut(inDir,false);

%--------------------------
params = [1,1,0.1,0.1];
methods = {'KP-Newton(ours)','FP-Newton','AKVF', 'CM'};
mid = [ sum(minmax(X(:,1)'))/2 sum(minmax(X(:,2)'))/2 sum(minmax(X(:,3)'))/2];   X = X- repmat(mid,size(X,1),1);
FN = cross( X(T(:,2),:)-X(T(:,1),:), X(T(:,3),:)-X(T(:,1),:) );X = X/sqrt(sum(sqrt(sum(FN.^2,2))/2));
initX = floater97(X,T);initX = complex(initX(:,1), initX(:,2));
triF = figure;
set(triF, 'Position', [100 100 480*3 480]);
subplot(1,4,4);
trimesh(T,X(:,1),X(:,2),X(:,3),'edgecolor', [.2 .2 .2]);set(gcf,'color','w');view(2);
subplot(1,4,3);
trimeshH = trimesh(T,X(:,1),X(:,2),X(:,3),'edgecolor', [.2 .2 .2]);set(gcf,'color','w');grid off; axis equal; 
view(2);drawnow;
for iMethod = 1:length(methods)
  subplot(1,4,3);
  if (iEnergy == 3 | iEnergy == 4 ) & strcmp(methods{iMethod},'CM')
    continue;
  end
  [~,statsAll] = deformMesh(X, T, [], [], initX, 1000, 1000, energies{iEnergy}, params(iEnergy), methods{iMethod},[],trimeshH);
  minE = min(statsAll(statsAll(:,end)>0,end))-1e-4;
  subplot(1,4,1);
  loglog(cumsum(statsAll(:,5)')/1000+1,(statsAll(:,end))'-minE,'-x','markers',4);
  hold on;
  ylabel(energies{iEnergy});
  xlabel('Time(s)');
  legend(methods);
  subplot(1,4,2);
  loglog((statsAll(:,end))'-minE,'-x','markers',4);
  hold on;
  legend(methods);
  xlabel('Iteration');
end
