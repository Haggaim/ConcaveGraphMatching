function  plotMeshAndPoints( V, F, params )
%===============================================================
% module:
% ------
% plotMeshAndPoints.m
%
% paper:
% -------
% Point registration via efficient convex relaxation.
% Haggai Maron, Nadav Dym, Itay Kezurer, Shahar Kovalsky,Yaron Lipman
%
% Description:
% -----------
% Plots a given mesh and points on it according to a predefined color map
%===============================================================

%--------------------------------------------
% Initialization
%--------------------------------------------
params.null = [];
n = getoptions(params,'n',1);
verInd = getoptions(params,'verInd',1:n);
buttonDownFuncHandle = getoptions(params,'funcHandle',@vertexCallback);
scattColor = getoptions(params,'scattColor','c');
scattSize = getoptions(params,'scattSize',500);
meshIdx = getoptions(params,'meshIdx',1);
labels = getoptions(params,'labels',{});

% check size of verInd compared to n
if n ~= size( verInd, 2 )
    %     warning('n = %d, size(verInd) = %d. Taking n to be %d.', n, size( verInd, 2 ), size( verInd, 2 ));
    n = size( verInd, 2 );
end
%============================================

%--------------------------------------------
% Take out NaNs
%--------------------------------------------
nanInd = isnan(verInd);
%============================================

%--------------------------------------------
% Extract verteces and faces information
%--------------------------------------------
axis equal; axis off;
if isempty(labels)
    labels = cellstr( num2str([1:n]') );
end

if ~isempty(F)
    patch('vertices',V','faces',F','FaceColor',[0.6 0.6 0.6],'EdgeColor','none','FaceAlpha',1,'ButtonDownFcn',{buttonDownFuncHandle,meshIdx});
else
    scatter3(V(1,:), V(2,:), V(3,:),10,[0.6 0.6 0.6]);
end
hold on
scatter3(V(1,verInd(~nanInd)), V(2,verInd(~nanInd)), V(3,verInd(~nanInd)),scattSize,scattColor,'filled');

axis equal, axis off
addRot3D;
%============================================

end




