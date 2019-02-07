function testGeodesicMatching(params)

params.null = [];
params.distanceType = getoptions(params,'distanceType','geodesic'); 
params.farthestType = 'euclidean';
params.initType = 'centroid';
params.doLineSearch = 0;
params.n = getoptions(params,'n',300);
params.numCoarseMathces = getoptions(params,'numCoarseMathces',3);
params.numAGDMaxPoints = getoptions(params,'numFPSPoints',8);
params.idx1 = getoptions(params,'idx1',14);
params.idx2 = getoptions(params,'idx2',15);
params.un = getoptions(params,'un',1500);
params.shape1name = '14.off';
params.shape2name = '15.off';
params.numSamples = getoptions(params,'numSamples',100);
params.beta=0.01;

matchShapes(params);
end


