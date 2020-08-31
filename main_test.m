% This is an attempt to write an optimisation code for a Sandwich Beam
% Author: Robert Jonsson
% Date: 2020-08-25

% Clear workspace
clear all
close all
clc

% Add paths to solvers
addpath('.') %path to material database [you could pick a different material database]
addpath('/Users/robertkth/Documents/GitHub/GCMMA-MMA-code-1.5') %path to GCMMA MATLAB functions
addpath('/Users/robertkth/Documents/GitHub/beamEB') %path to constraint solver [you could add a different solver]
addpath('/Users/robertkth/Documents/GitHub/Timoshenko') %path to constraint solver [you could add a different solver]
addpath('/Users/robertkth/Documents/GitHub/DriveCycles') %path to drive cycles, for use phase.
import ecoOptimize.*

% Import the material data to be worked with
global materialsData
materialsData = importdata('materialData.mat');

% Define the model to be worked with
global model
model = initmodelSWBeam; %The initial model is defined.
model=blendMaterials(model,materialsData);
model=updateDependentVars(model);
figure(1), clf, dispModel(model,0)

% Define optimisation parameters
xval=[model.B(1) model.H(1) model.H(2) model.H(3)]';
xnam={'B(1)' 'H(1)' 'H(2)' 'H(3)'}';
xmin=[0.05 0.001 0.01 0.001]';
xmax=[0.3 0.1 0.2 0.1]';
maxiter=20;

% initiate GCMMA
gcmma=GCMMA.init(@ecoOptimize.optFuncs,xval,xnam,xmin,xmax);

% run GCMMA
disp(['Optimizing for: ',model.objfunc])
gcmma.displive=1;
%figure(2), clf, gcmma.plotlive=1;
gcmma.maxoutit=20;
[gcmma,xval]=GCMMA.run(gcmma);
[f0val,fval]=optFuncs(xval,xnam,false);

% view results
figure(2), clf, GCMMA.plotIter(gcmma)
figure(1), dispModel(model,1)
mass=computeMass(model)
LCE=computeLCE(model)
LCCO2=computeLCCO2(model)
LCCost=computeLCCost(model)

%% Functions

function model = initmodelSWBeam %This is used to set up the initial model
model.objfunc='LCE'; %Objective function
model.fmax=[5e-3]; %Maximum deflection
model.driveDistTotal=100e3; %Driving distance of vehicle
model.solver='beamTSAna'; %Solver used
model.loadcase='cantilever_pt'; %Load case (pt = point load)
model.P=-1e4; %Applied load
model.xP=0.5; %Spanwise fraction of beam where load is applied
model.L=1; %Beam length
model.xsection='sandwich'; %Section type
model.B=0.3; %Cross-section breadth
model.H=[0.05 0.05 0.05]; %Beam thicknesses [lower face, core, upper face].
model.material={'CFRP' 'PVC' 'CFRP'}; %Beam materials 
model.alpha=[1 1 1]; %Material mixtures.
model.drivecycle = 'WLTP';
end

function dispModel(model,fill)
N=max([numel(model.B) numel(model.H)]);
B=repmat(model.B,1,N-numel(model.B)+1);
H=repmat(model.H,1,N-numel(model.H)+1);
H0=([0 cumsum(H(1:end-1))]-sum(H)/2);
B0=-B/2;
C=colormap;
for i=1:N
  R=rectangle('Position',[B0(i) H0(i) B(i) H(i)]);
  if fill
    set(R,'Facecolor',C(50*(i-1)+1,:))
  end
end
axis([B0(i)-B(i)/10 B0(i)+B(i)+B(i)/10 H0(i)-H(i)/10 H0(i)+H(i)+H(i)/10])
axis auto, axis equal
xlabel('b [m]'), ylabel('h [m]')
end


