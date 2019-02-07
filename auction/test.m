function test()
%% DEMO
	% compile mex file
	mex -largeArrayDims auctionAlgorithmSparseMex.cpp -lut
	
	% create sample data
	N = 100;
	
	A = rand(N,N);
	
	% create sparse matrix, since the mex implementation uses the Matlab
	% sparse matrix data structure
	A = sparse(A);

	% scale A such that round(Ascaled) has sufficient accuracy
	scalingFactor = 10^6;
	Ascaled = A*scalingFactor;
	
	% solve assignment problem
	tic
	[assignments,P] = sparseAssignmentProblemAuctionAlgorithm(Ascaled);
	toc
end