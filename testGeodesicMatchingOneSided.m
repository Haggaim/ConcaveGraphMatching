function testGeodesicMatchingOneSided(params)
params.null = [];
params.k = getoptions(params,'k',100);
params.n = getoptions(params,'n',100);
params.numCoarseMathces = getoptions(params,'numCoarseMathces',3);
params.numFPSPoints = getoptions(params,'numFPSPoints',8);
params.idx1 = getoptions(params,'idx1',14);
params.idx2 = getoptions(params,'idx2',15);
params.uk = getoptions(params,'uk',500);
params.un = getoptions(params,'un',500);
params.shape1name = '14.off';
params.shape2name = '15.off';
params.doLineSearch = getoptions(params,'doLineSearch',0);
params.numSamples = getoptions(params,'numSamples',100);
matchShapesOneSided(params)
end