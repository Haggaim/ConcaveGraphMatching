function matchShapesOneSided(params)
%% preprocess
shape1name = getoptions(params,'shape1name','');
shape2name = getoptions(params,'shape2name','');
k = getoptions(params,'k',50);
n = getoptions(params,'n',100);
uk = getoptions(params,'uk',500);
un = getoptions(params,'un',500);
seed = getoptions(params,'seed',30);
outdir = getoptions(params,'outdir','');
dosave = getoptions(params,'dosave',0);
params.dosave = getoptions(params,'dosave',0);

rng(seed)

% preprocess shape 1
[V1,F1] = read_mesh(shape1name);
if size(V1,1)<size(V1,2), V1=V1';F1=F1';end
V1 = V1/sqrt(CORR_calculate_area(F1,V1));
adj1 = triangulation2adjacency(F1,V1);

% preprocess shape 2
[V2,F2] = read_mesh(shape2name);
if size(V2,1)<size(V2,2), V2=V2';F2=F2';end
V2 = V2/sqrt(CORR_calculate_area(F2,V2));
adj2 = triangulation2adjacency(F2,V2);


dist1 = graphallshortestpaths(adj1);
dist2 = graphallshortestpaths(adj2);

params.dist1 = dist1;
params.dist2 = dist2;

idx1 = chooseFarthestPoints(squareform(pdist(V1)),k);
idx2 = chooseFarthestPoints(squareform(pdist(V2)),n);

%% project with our functional
[X_proj_new,Efinal] = frankWolfeSamplingSolverOneSided(V1,F1,idx1,dist1,V2,F2,idx2,dist2,params);

%% visualize
figure, imagesc(X_proj_new),title('Assignment matrix')  
[~,mappedIdx2] = max(X_proj_new,[],2);
mappedIdx2 = idx2(mappedIdx2);
PlotResultAfterLocalMinimization(V1',F1',V2',F2',idx1,mappedIdx2,'source','target')
savefig('coarse_a_',params)
PlotResultAfterLocalMinimization(V2',F2',V1',F1',mappedIdx2,idx1,'target','source')
savefig('coarse_b_',params)


%% upsample
if ~isfield(params,'gt1')
    uidx1 = chooseFarthestPoints(squareform(pdist(V1)),uk);
    uidx2 = chooseFarthestPoints(squareform(pdist(V2)),un);
else
    uidx1 = addFartherstPoints(squareform(pdist(V1)),uk,params.gt1);
    uidx2 = chooseFarthestPoints(squareform(pdist(V2)),un)';
end
% generate initialization assignment
[~,NNidx1] = min(dist1(uidx1,idx1),[],2);
[~,NNidx2] = min(dist2(uidx2,mappedIdx2));
X0 = zeros(uk,un);
X0(sub2ind([uk un],1:uk,NNidx2(NNidx1)))=1;

params.Xinit = X0;

[uX_proj_new,Efinal]= frankWolfeSolverSmartOnesided(dist1(uidx1,uidx1),dist2(uidx2,uidx2),params);
%% visualize upsampled solution
figure, imagesc(uX_proj_new),title('Assignment matrix: after upsampling')  
[~,umappedIdx2] = max(uX_proj_new,[],2);
umappedIdx2 = uidx2(umappedIdx2);
PlotResultAfterLocalMinimization(V1',F1',V2',F2',uidx1,umappedIdx2,'source','target')
savefig('upsampled_a_',params);
PlotResultAfterLocalMinimization(V2',F2',V1',F1',umappedIdx2,uidx1,'target','source')
savefig('upsampled_b_',params);

if dosave
    [~,name1,~]=fileparts(params.shape1name);
    [~,name2,~]=fileparts(params.shape2name);
    if size(uidx1,1)>size(uidx1,2) uidx1 = uidx1';end
    if size(umappedIdx2,1)>size(umappedIdx2,2) umappedIdx2 = umappedIdx2';end
    dlmwrite(fullfile(outdir,['coarse_idx_', name1 '_' name2 '.txt']),[idx1',mappedIdx2']);
    dlmwrite(fullfile(outdir,['upsampled_idx_', name1 '_' name2 '.txt']),[uidx1',umappedIdx2']);
    params.dist = [];
    params.dist1 = [];
    params.dist2 = [];
    save(fullfile(params.outdir,['res_', name1 '_' name2 '.mat']),'Efinal','params')
end
end

function savefig(prefix,params)
if params.dosave
    [~,name1,~]=fileparts(params.shape1name);
    [~,name2,~]=fileparts(params.shape2name);
    saveas(gcf,fullfile(params.outdir,[prefix, name1 '_' name2 '.fig']))
    saveas(gcf,fullfile(params.outdir,[prefix, name1 '_' name2 '.png']))
end
end