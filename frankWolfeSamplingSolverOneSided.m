function  [X_proj_new,Efinal] = frankWolfeSamplingSolverOneSided(V1,F1,idx1,dist1,V2,F2,idx2,dist2,params)
numSamples = getoptions(params,'numSamples',5);
k = numel(idx1);
n = numel(idx2);

% find feature points and sample triplets
params.dist = params.dist1;
samples1 = getPointSamplesByFPS(V1,F1,params);
params.dist = params.dist2;
samples2 = getPointSamplesByFPS(V2,F2,params);
params.dist = [];
% find NN in the set of chosen points
newSamples1 = zeros(size(samples1));
for ii = 1:numel(samples1)
    [~,newSamples1(ii)] = min(dist1(idx1,samples1(ii)));
end
newSamples1 = idx1(newSamples1);

newSamples2 = zeros(size(samples2));
for ii = 1:numel(samples2)
    [~,newSamples2(ii)] = min(dist2(idx2,samples2(ii)));
end
newSamples2 = idx2(newSamples2);

% calculate initializations
Xinits = {};
randSampleIdx1 = randi(size(newSamples1,1),[numSamples 1]);
randSampleIdx2 = randi(size(newSamples2,1),[numSamples 1]);

for ii = 1:numSamples
    % generate initialization assignment
    [~,NNidx1] = min(dist1(idx1,newSamples1(randSampleIdx1(ii),:)),[],2);
    [~,NNidx2] = min(dist2(idx2,newSamples2(randSampleIdx2(ii),:)));
    X0 = zeros(k,n);
    X0(sub2ind([k n],1:k,NNidx2(NNidx1)))=1;
    Xinits = [Xinits {X0}];    
end

dist1reduced = dist1(idx1,idx1);
dist2reduced = dist2(idx2,idx2);


% calc translation vector
params.X_trans = -2*sum(dist1reduced,2)*sum(dist2reduced,1)+...
    sum(dist1reduced.^2,2)*n*ones(1,n)+...
    ones(k,1)*k*sum(dist2reduced.^2,1);
params.row_max = max(params.X_trans,[],2);
params.X_trans = repmat(params.row_max,[1,n]);


% run optimization
objs = [];
X_proj_new = cell(1);
for ii = 1:numel(Xinits)
    if mod(ii,50)==0 fprintf('running on Xinit number %d\n',ii); end
    params.Xinit = Xinits{ii};
    [X_proj_new{ii},objs(ii)]= frankWolfeSolverSmartOnesided(dist1reduced,dist2reduced,params);
end

% take best map
[Efinal,minidx] = min(objs);
X_proj_new = X_proj_new{minidx};
disp('-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------')
fprintf('Frank-Wolfe sampler-solver: Optimization has finished. %d initializations were tried, Efinal = %.3f \n',numel(objs),Efinal);
disp('-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------')
end
