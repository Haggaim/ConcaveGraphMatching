function  [X_proj_new,Efinal] = frankWolfeSamplingSolver(V1,F1,idx1,AFF1,V2,F2,idx2,AFF2,params)
numSamples = getoptions(params,'numSamples',10);
geodist1 = getoptions(params,'geodist1',[]);
geodist2 = getoptions(params,'geodist2',[]);

n = numel(idx1);

% find feature points and sample triplets
params.dist = params.geodist1;
samples1 = getPointSamplesByFPS(V1,F1,params);
params.dist = params.geodist2;
samples2 = getPointSamplesByFPS(V2,F2,params);
params.dist = [];

% find NN in the set of chosen points
newSamples1 = zeros(size(samples1));
for ii = 1:numel(samples1)
    [~,newSamples1(ii)] = min(geodist1(idx1,samples1(ii)));
end
newSamples1 = idx1(newSamples1);
newSamples2 = zeros(size(samples2));
for ii = 1:numel(samples2)
    [~,newSamples2(ii)] = min(geodist2(idx2,samples2(ii)));
end
newSamples2 = idx2(newSamples2);

% calculate initializations
Xinits = {};

if size(newSamples1,1)<size(newSamples1,2) newSamples1=newSamples1';end
if size(newSamples2,1)<size(newSamples2,2) newSamples2=newSamples2';end

randSampleIdx1 = randi(size(newSamples1,1),[numSamples 1]);
randSampleIdx2 = randi(size(newSamples2,1),[numSamples 1]);


for ii = 1:numSamples
    % generate initialization assignment
    [~,NNidx1] = min(geodist1(idx1,newSamples1(randSampleIdx1(ii),:)),[],2);
    [~,NNidx2] = min(geodist2(idx2,newSamples2(randSampleIdx2(ii),:)));
    X0 = zeros(n);
    targetIdx = NNidx2(NNidx1);
    if size(targetIdx,1)>size(targetIdx,2) targetIdx = targetIdx';end
    X0(sub2ind([n n],1:n,targetIdx))=1;
    Xinits = [Xinits {X0}];
end

aff1reduced = AFF1(idx1,idx1);
aff2reduced = AFF2(idx2,idx2);

% run optimization
objs = [];
X_proj_new = cell(1,numel(Xinits));
for ii = 1:numel(Xinits)
    if mod(ii,50)==0 fprintf('running on Xinit number %d\n',ii); end
    params.Xinit = Xinits{ii};
   [X_proj_new{ii},objs(ii)] = frankWolfeSolver(aff1reduced,aff2reduced,params);
end
% take map with lowest energy
[Efinal,minidx] = min(objs);
X_proj_new = X_proj_new{minidx};
disp('-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------')
fprintf('Frank-Wolfe sampler-solver: Optimization has finished. %d initializations were tried , Efinal = %.3f \n',numel(objs),Efinal);
disp('-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------')
end
