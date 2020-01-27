% This file is part of LEnOP, a code to optimize a design model for 
% minimum life cycle energy subject to functional requirements.
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
% main_beam_gcmma.m

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
  model.xF=1;
  model.L=2;
  model.xsection='layered';
  model.B=1;
  model.H=[0.05 0.05 0.05];
  model.material={'CFRP' 'PUR' 'CFRP'};
  model.alpha=[1 1 1];
  material=LEnOpFunctions.blendMaterials(model,materialsData);
  model=LEnOpFunctions.updateMaterialProps(model,material);
  model=LEnOpFunctions.updateDependentVars(model);
  
  %% set optimisation params
  xval=[model.H(1) model.H(2) model.H(3)]';
  xnam={'H(1)' 'H(2)' 'H(3)'};
  xmin=[0.001 0.001 0.001]';
  xmax=[0.2 0.2 0.2]';
  maxiter=10;
  
  %% initiate GCMMA
  gcmma=GCMMAFunctions.init(xval,xnam,xmin,xmax,maxiter);

end

%% run GCMMA
[gcmma,xval]=GCMMAFunctions.run(gcmma,xval,xnam,xmin,xmax,true);
[f0val,fval]=LEnOpFunctions.optFunctions(xval,xnam,false);

%% plot convergence
figure(1), clf, GCMMAFunctions.plotIter(gcmma)

%% check results
mass=LEnOpFunctions.computeMass(model)
LCE=LEnOpFunctions.computeLCE(model)
beam=computeEulerBernoulli(model);
fval=max(abs(beam.w))
comsol=runCOMSOLBeam(model);
v=mpheval(comsol,'v','edim',1,'dataset','dset1');
figure(2), clf, plot(beam.x,beam.w,v.p(1,:),v.d1,'*')
