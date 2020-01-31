% This file is part of ecoOptimize, a code to optimize a design model for 
% minimum eco impacts subject to functional requirements.
% 
% Copyright (C) 2020 Ciarán O'Reilly <ciaran@kth.se>
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
% 
% main.m

restart=1;

if restart
  
  clear all
  clear global
  addpath('.') %path to material database [you could pick a different material database]
  addpath('../beamEB') %path to constraint solver [you could add a different solver]
  
  %% load material database
  global materialsData
  materialsData=importdata('materialData.mat');
    
  %% initiate model of the panel
  global model material
  model=initModelBeam;
  material=ecoOptimizeFuncs.blendMaterials(model,materialsData);
  model=ecoOptimizeFuncs.updateMaterialProps(model,material);
  model=ecoOptimizeFuncs.updateDependentVars(model);
  figure(1), clf, dispModel(model,0)
  
  %% set optimisation params
  xval=[model.B(1) model.B(2) model.B(3) model.H(1) model.H(2) model.H(3)]';
  xnam={'B(1)' 'B(2)' 'B(3)' 'H(1)' 'H(2)' 'H(3)'}';
  xmin=[0.05 0.005 0.05 0.005 0.005 0.005]';
  xmax=[0.2 0.2 0.2 0.2 0.2 0.2]';
  maxiter=20;
  
  %% initiate GCMMA
  gcmma=GCMMAFuncs.init(xval,xnam,xmin,xmax,maxiter);

end

%% run GCMMA
disp(['Optimizing for: ',model.objfunc])
figure(2)
[gcmma,xval]=GCMMAFuncs.run(gcmma,xval,xnam,xmin,xmax,true);
[f0val,fval]=ecoOptimizeFuncs.optFuncs(xval,xnam,false);

%% view results
figure(2), clf, GCMMAFuncs.plotIter(gcmma)
figure(1), dispModel(model,1)
mass=ecoOptimizeFuncs.computeMass(model)
LCE=ecoOptimizeFuncs.computeLCE(model)
LCCO2=ecoOptimizeFuncs.computeLCCO2(model)
LCCost=ecoOptimizeFuncs.computeLCCost(model)

%%
function model=initModelBeam
  model.objfunc='LCE';
  model.fmax=[1e-3];
  model.fscale=[1e3];
  model.driveDistTotal=1e5;
  model.solver='beamEBComsol';
  model.loadcase='simple_pt';
  model.P=-1e4;
  model.xP=0.5;
  model.L=1;
  model.xsection='layered';
  model.B=[0.15 0.15 0.15];
  model.H=[0.05 0.05 0.05];
  model.material={'Steel' 'Al-alloys' 'Steel'};
  model.alpha=[1 1 1];
end

%%
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