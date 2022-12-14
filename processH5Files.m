%% This script loads h5 files and processes the data
clc,clear;
close all;


%% Load .h5 file from directory
dir = "/home/shengduo/pylith-developer/build/debug/pylith-nonRegSlipLawWithVaryingB/examples/bar_shearwave/quad4/output/fault/";
FourierTerms = 16;

% As = [0.008, 0.01, 0.012, 0.014, 0.016];
% Bs = [0.008, 0.01, 0.012, 0.014, 0.016];

% As = [0.002, 0.004, 0.006, 0.008, 0.01];
% Bs = [0.006 0.008 0.01 0.012 0.014 0.016];

As = [0.006 0.0064 0.0068 0.0072 0.0076 0.008];
Bs = [0.012 0.0124 0.0128 0.0132 0.0136 0.014];

% As = [0.011, 0.015];
% Bs = [0.011, 0.015];

% As = [0.003, 0.007];
% Bs = [0.009, 0.013];

nOfFiles = length(As) * length(Bs);
saveDir = "/home/shengduo/InverseProblems/GPRWorkingField/data/";
saveFileName = "trainGrid1102_1.mat";

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

%% Calculate Fourier coefficients (cos(k pi/T t)) for the point histories
observations = zeros(nOfFiles, FourierTerms * nOfNodes);
for m = 1 : 1 : FourierTerms
    for fileNo = 1 : 1 : nOfFiles
        nOfTSteps = length(ts(fileNo, :));
        nOfNodes = size(Vxs(fileNo, :, :), 2);
        T = max(ts(fileNo, :));
        kPiTt = (m - 1) * pi / T * ts(fileNo, :);
        cosines = cos(kPiTt);
        observations(fileNo, (m - 1) * nOfNodes + 1 : m * nOfNodes) = ...
            trapz(ts(fileNo, :), (reshape(Vxs(fileNo, :, :), [nOfNodes, nOfTSteps]) .* cosines)');
    end
end

% Save Nsave Vsave
save(saveDir + saveFileName, 'Us', 'observations');