%% This script loads h5 files and processes the data
clc,clear;
close all;


%% Load .h5 file from directory
dir = "/home/shengduo/pylith-developer/build/debug/pylith-nonRegSlipLawWithVaryingB/examples/bar_shearwave/quad4/output/fault/";
As = [0.008, 0.01, 0.012, 0.014, 0.016];
Bs = [0.008, 0.01, 0.012, 0.014, 0.016];

% As = [0.011, 0.015];
% Bs = [0.011, 0.015];

nOfFiles = length(As) * length(Bs);
saveDir = "/home/shengduo/InverseProblems/GPRWorkingField/data/";
saveFileName = "trainGrid.mat";

% Get sizes
filename = "A" + string(As(1)) + "_B" + string(Bs(1)) + "-fault.h5";
faultFileName = dir + filename;
t = h5read(faultFileName, '/time');
t = reshape(t, [1, size(t, 3)]);
XYZs = h5read(faultFileName, '/geometry/vertices');

SlipRate = h5read(faultFileName, '/vertex_fields/slip_rate');
nOfNodes = size(SlipRate, 2);
nOfTSteps = length(t);

% Allocate memory for data matrices 
ts = zeros(nOfFiles, nOfTSteps);
Vxs = zeros(nOfFiles, nOfNodes, nOfTSteps);
Vys = zeros(nOfFiles, nOfNodes, nOfTSteps);
Us = zeros(nOfFiles, 2);

fileNo = 1;
for i = 1:1:length(As)
    for j = 1:1:length(Bs)
        filename = "A" + string(As(i)) + "_B" + string(Bs(j)) + "-fault.h5";
        faultFileName = dir + filename;
        t = h5read(faultFileName, '/time');
        t = reshape(t, [1, size(t, 3)]);
        SlipRate = h5read(faultFileName, '/vertex_fields/slip_rate');
        nOfNodes = size(SlipRate, 2);
        nOfTSteps = length(t);
        Xs = h5read(faultFileName, '/geometry/vertices');
        Vx = reshape(SlipRate(1, :, :), [nOfNodes, nOfTSteps]);
        Vy = reshape(SlipRate(2, :, :), [nOfNodes, nOfTSteps]);
        
        ts(fileNo, :) = t;
        Us(fileNo, :) = [As(i), Bs(j)];
        Vxs(fileNo, :, :) = Vx;
        Vys(fileNo, :, :) = Vy;
        fileNo = fileNo + 1;
    end
end

%% Do moment-time integral for the points
momentTerms = 4;
observations = zeros(nOfFiles, momentTerms * nOfNodes);
for m = 1 : 1 : momentTerms
    for fileNo = 1 : 1 : nOfFiles
        observations(fileNo, (m - 1) * nOfNodes + 1 : m * nOfNodes) = ...
            trapz(t, power(reshape(Vxs(fileNo, :, :), [nOfNodes, nOfTSteps])', m));
    end
end

% Save Nsave Vsave
save(saveDir + saveFileName, 'Us', 'observations');