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
% along with this program.  If not, see <https://www.gnu.org/licenses/>.% 
% 
% main.m

restart=1;

if restart
  
  clear all
  clear global
  addpath('.') %path to material database
  addpath('../beamEB') %path to constraint solver
  
  %% select a material
  global materialsData
  materialsData=importdata('materialData.mat');
  
  %% initiate model of the panel
  global model material
  model.loadcase='simple_pt';
  model.F=-1e4;
  model.xF=1;
  model.L=2;
  model.xsection='layered';
  model.B=1;
  model.H=[0.05 0.05 0.05];
  model.material={'GFRP' 'PUR' 'GFRP';'CFRP' 'PVC' 'CFRP'};
  model.alpha=[0.3 0.4 0.5];
  material=ecoOptimizeFuncs.blendMaterials(model,materialsData);
  model=ecoOptimizeFuncs.updateMaterialProps(model,material);
  model=ecoOptimizeFuncs.updateDependentVars(model);
  figure(1), clf, subplot(2,1,1), ecoOptimizeFuncs.dispModel(model)
  
  %% set optimisation params
  xval=[model.H(1) model.H(2) model.H(3) model.alpha(1,1) model.alpha(1,2) model.alpha(1,3)]';
  xnam={'H(1)' 'H(2)' 'H(3)' 'alpha(1,1)' 'alpha(1,2)' 'alpha(1,3)'};
  xmin=[0.001 0.001 0.001 0 0 0]';
  xmax=[0.2 0.2 0.2 1 1 1]';
  maxiter=30;
  
  %% initiate GCMMA
  gcmma=GCMMAFuncs.init(xval,xnam,xmin,xmax,maxiter);

end

%% run GCMMA
figure(2)
[gcmma,xval]=GCMMAFuncs.run(gcmma,xval,xnam,xmin,xmax,true);
[f0val,fval]=ecoOptimizeFuncs.optFuncs(xval,xnam,false);

%% plot results
figure(1), subplot(2,1,2), ecoOptimizeFuncs.dispModel(model)
figure(2), clf, GCMMAFuncs.plotIter(gcmma)

%% check results
mass=ecoOptimizeFuncs.computeMass(model)
LCE=ecoOptimizeFuncs.computeLCE(model)
beam=computeEulerBernoulli(model);
fval=max(abs(beam.w))
figure(3), clf, plot(beam.x,beam.w), xlabel('l [m]'), ylabel('w [m]')
% comsol=runCOMSOLBeam(model);
% v=mpheval(comsol,'v','edim',1,'dataset','dset1');
% figure(3), clf, plot(beam.x,beam.w,v.p(1,:),v.d1,'*'), xlabel('x [m]'), xlabel('w [m]')
