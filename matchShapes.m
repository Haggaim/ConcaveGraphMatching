function energy = matchShapes(params)
seed = getoptions(params,'seed',1);
rng(seed)
shape1name = getoptions(params,'shape1name','');
shape2name = getoptions(params,'shape2name','');
distanceType = getoptions(params,'distanceType','geodesic');
doSaveEnergy = getoptions(params,'doSaveEnergy',0);
n = getoptions(params,'n',30);
un = getoptions(params,'un',100);
params.dosave = getoptions(params,'dosave',0);
dovis = getoptions(params,'dovis',1);
posDefRadial30Const = getoptions(params,'posDefRadial30Const',1);

% preprocess shape 1
[V1,F1] = read_mesh(shape1name);
if size(V1,1)<size(V1,2), V1=V1';F1=F1';end
V1 = V1(:,1:3);
if ~isempty(F1)
    V1 = V1/sqrt(CORR_calculate_area(F1,V1));
    adj1 = triangulation2adjacency(F1,V1);
end
% preprocess shape 2
[V2,F2] = read_mesh(shape2name);
if size(V2,1)<size(V2,2), V2=V2';F2=F2';end
V2 = V2(:,1:3);
if ~isempty(F1)
    V2 = V2/sqrt(CORR_calculate_area(F2,V2));
    adj2 = triangulation2adjacency(F2,V2);
end


switch distanceType
    case 'geodesicRadial31'
        AFF1 = graphallshortestpaths(adj1);
        AFF2 = graphallshortestpaths(adj2);
        AFF1 = radial31(AFF1,posDefRadial30Const);
        params.AFF1 = AFF1;
        AFF2 = radial31(AFF2,posDefRadial30Const);
        params.AFF2 = AFF2;
    case 'geodesicRadial30'
        AFF1 = graphallshortestpaths(adj1);
        AFF2 = graphallshortestpaths(adj2);
        AFF1 = radial30(AFF1,posDefRadial30Const);
        params.AFF1 = AFF1;
        AFF2 = radial30(AFF2,posDefRadial30Const);
        params.AFF2 = AFF2;
        
    case 'euclidean'
        % calc distances
        AFF1 = squareform(pdist(V1));
        AFF2 = squareform(pdist(V2));
    case 'geodesic'
        AFF1 = graphallshortestpaths(adj1);
        AFF2 = graphallshortestpaths(adj2);
        
    case 'biharmonic'
        L1 = cotmatrix(V1,F1);
        M1 = massmatrix(V1,F1);
        L2 = cotmatrix(V2,F2);
        M2 = massmatrix(V2,F2);
        
        [ U1, E1 ] = eigs(L1,M1, 500,'sm');
        [ U2, E2 ] = eigs(L2,M2, 500,'sm');
        %throw away first e.v.
        U1 = U1(:,2:end);
        E1 = E1(2:end,2:end);
        U2 = U2(:,2:end);
        E2 = E2(2:end,2:end);
        
        % scale by sqrt e.v.
        U1 = U1*diag(1./diag((abs(E1))));
        U2 = U2*diag(1./diag((abs(E2))));
        
        
        % calc distances
        AFF1 = squareform(pdist(U1));
        AFF2 = squareform(pdist(U2));
        
    case 'gaussiangeodesic'
        AFF1 = exp(-graphallshortestpaths(adj1).^2/params.beta);
        AFF2 = exp(-graphallshortestpaths(adj2).^2/params.beta);
end
params.AFF1 = AFF1;
params.AFF2 = AFF2;

idx1 = chooseFarthestPoints(squareform(pdist(V1)),n);
idx2 = chooseFarthestPoints(squareform(pdist(V2)),n);

params.idx1 = idx1;
params.idx2 = idx2;


if strcmp('geodesic',distanceType) || isempty(F1)
    params.geodist1 = AFF1;
    params.geodist2 = AFF2;
else
    adj1 = triangulation2adjacency(F1,V1);
    params.geodist1 = graphallshortestpaths(adj1);
    adj2 = triangulation2adjacency(F2,V2);
    params.geodist2 = graphallshortestpaths(adj2);
end


%% project with our functional
X_proj_new = frankWolfeSamplingSolver(V1,F1,idx1,AFF1,V2,F2,idx2,AFF2,params);

%% visualize
[~,mappedIdx2] = max(X_proj_new,[],2);
mappedIdx2 = idx2(mappedIdx2);
if dovis
    figure, imagesc(X_proj_new), title('Assignment matrix')
    PlotResultAfterLocalMinimization(V1',F1',V2',F2',idx1,mappedIdx2,'source','target')
    savefig('coarse_a_',params)
    PlotResultAfterLocalMinimization(V2',F2',V1',F1',mappedIdx2,idx1,'target','source')
    savefig('coarse_b_',params)
end

%% upsample

uidx1 = chooseFarthestPoints(squareform(pdist(V1)),un);
uidx2 = chooseFarthestPoints(squareform(pdist(V2)),un);
% generate initialization assignment
[~,NNidx1] = min(params.geodist1(uidx1,idx1),[],2);
[~,NNidx2] = min(params.geodist2(uidx2,mappedIdx2));
X0 = zeros(un,un);
X0(sub2ind([un un],1:un,NNidx2(NNidx1)))=1;
% params.optimalTranslation = (un/(1000*n))*optimalTranslation;

params.Xinit = X0;
uX_proj_new= frankWolfeSolver(AFF1(uidx1,uidx1),AFF2(uidx2,uidx2),params);


%% visualize upsampled solution
if dovis
    figure, imagesc(uX_proj_new),title('Assignment matrix: after upsampling')
    [~,umappedIdx2] = max(uX_proj_new,[],2);
    umappedIdx2 = uidx2(umappedIdx2);
    
    PlotResultAfterLocalMinimization(V1',F1',V2',F2',uidx1,umappedIdx2,'source','target')
    savefig('upsampled_a_',params);
    PlotResultAfterLocalMinimization(V2',F2',V1',F1',umappedIdx2,uidx1,'target','source')
    savefig('upsampled_b_',params);
end

% calc energy
X = uX_proj_new;
energy = trace(ones(un)*X*AFF2(uidx2,uidx2).^2*X')-2*trace(AFF1(uidx1,uidx1)*X*AFF2(uidx2,uidx2)*X')+trace(AFF1(uidx1,uidx1).^2*X*ones(un)*X');



if params.dosave
    [~,name1,~]=fileparts(params.shape1name);
    [~,name2,~]=fileparts(params.shape2name);
    if size(uidx1,1)>size(uidx1,2) uidx1 = uidx1';end
    if size(umappedIdx2,1)>size(umappedIdx2,2) umappedIdx2 = umappedIdx2';end
    
    dlmwrite(fullfile(params.outdir,['coarse_idx_', name1 '_' name2 '.txt']),[idx1',mappedIdx2']);
    dlmwrite(fullfile(params.outdir,['upsampled_idx_', name1 '_' name2 '.txt']),[uidx1',umappedIdx2']);
end

if doSaveEnergy
    [~,name1,~]=fileparts(params.shape1name);
    [~,name2,~]=fileparts(params.shape2name);
    save(fullfile(params.outdir,['energy_', name1 '_' name2 '.mat']),'energy')
end
end

function savefig(prefix,params)
if params.dosave
    [~,name1,~]=fileparts(params.shape1name);
    [~,name2,~]=fileparts(params.shape2name);
    saveas(gcf,fullfile(params.outdir,[prefix, name1 '_' name2 '.fig']))
    
end
end

