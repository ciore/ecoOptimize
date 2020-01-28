% This file is part of ecoOptimize, a code to optimize a design model for 
% minimum eco impacts subject to functional requirements.
% 
% Copyright (C) 2018 Ciarán O'Reilly <ciaran@kth.se>
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
  
  %% select a material
  addpath('.')
  global materialsData
  materialsData=load('materialData.mat');
  materialsData=materialsData.data;
  
  %% initiate model of the panel
  global model material
  model.loadcase='simple_pt';
  model.F=-1e4;
  model.xF=0.5;
  model.L=1;
  model.xsection='layered';
  model.B=[0.15 0.15 0.15];
  model.H=[0.05 0.05 0.05];
  model.material={'Steel' 'Steel' 'Steel'};
  model.alpha=[1 1 1];
  material=ecoOptimizeFuncs.blendMaterials(model,materialsData);
  model=ecoOptimizeFuncs.updateMaterialProps(model,material);
  model=ecoOptimizeFuncs.updateDependentVars(model);
  figure(1), clf, subplot(1,2,1), ecoOptimizeFuncs.dispModel(model)

  
  %% set optimisation params
  xval=[model.B(1) model.B(2) model.B(3) model.H(1) model.H(2) model.H(3)]';
  xnam={'B(1)' 'B(2)' 'B(3)' 'H(1)' 'H(2)' 'H(3)'};
  xmin=[0.05 0.01 0.05 0.01 0.01 0.01]';
  xmax=[0.2 0.2 0.2 0.2 0.2 0.2]';
  maxiter=10;
  
  %% initiate GCMMA
  gcmma=GCMMAFuncs.init(xval,xnam,xmin,xmax,maxiter);

end

%% run GCMMA
figure(2)
[gcmma,xval]=GCMMAFuncs.run(gcmma,xval,xnam,xmin,xmax,true);
[f0val,fval]=ecoOptimizeFuncs.optFunctions(xval,xnam,false);

%% plot result
figure(1), subplot(1,2,2), ecoOptimizeFuncs.dispModel(model)
figure(2), clf, GCMMAFuncs.plotIter(gcmma)

%% check results
mass=ecoOptimizeFuncs.computeMass(model)
LCE=ecoOptimizeFuncs.computeLCE(model)
beam=computeEulerBernoulli(model);
figure(3), clf, plot(beam.x,beam.w), xlabel('x [m]'), xlabel('w [m]')
fval=max(abs(beam.w))
% comsol=runCOMSOLBeam(model);
% v=mpheval(comsol,'v','edim',1,'dataset','dset1');
% figure(3), clf, plot(beam.x,beam.w,v.p(1,:),v.d1,'*'), xlabel('x [m]'), xlabel('w [m]')
