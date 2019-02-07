function [samples] = getPointSamplesByFPS(V,F,params)
% V,F mesh representation
params.null = [];
n = getoptions(params, 'numFPSPoints', 5);
numCoarseMathces = getoptions(params,'numCoarseMathces',0);
dist = getoptions(params,'dist',[]);
idx = chooseFarthestPoints(dist,n);

% choose all triplets
samplesIdx = combinator(numel(idx),numCoarseMathces,'p');
samples = idx(samplesIdx);

end