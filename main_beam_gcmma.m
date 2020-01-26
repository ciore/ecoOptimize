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
  materialsData=load('materialData2.mat');
  materialsData=materialsData.data;
  global material
  matInd=[1];
  fnames=fieldnames(materialsData);
  for i=1:numel(fnames)
    eval(['material.',fnames{i},'=materialsData.',fnames{i},'(matInd);'])
  end
  
  %% initiate model of the panel
  global model
  model.F=-1e4;
  model.xF=0.5;
  model.L=2;
  model.B=[0.1 0.1 0.1];
  model.H=[0.05 0.05 0.05];
  model.El=[2.1e11 2.1e11 2.1e11];
  model.nul=[0.285 0.285 0.285];
  model.rhol=[7.8e3 7.8e3 7.8e3];
  model.xsection='layered';
  model.loadcase='simple_pt';
  model=LEnOpFunctions.updateDependentVars(model);

%   %% map materials to model
%   model.E=material.youngsModulus{1};
%   model.rho=material.density{1};
%   model.nu=material.poissonsRatio{1};
%   model=LEnOpFunctions.updateDependentVars(model);
  
  %% set optimisation params
  xval=[model.B(1) model.B(2) model.B(3) model.H(1) model.H(2) model.H(3)]';
  xnam={'B(1)' 'B(2)' 'B(3)' 'H(1)' 'H(2)' 'H(3)'};
  xmin=[0.05 0.01 0.05 0.01 0.01 0.01]';
  xmax=[0.2 0.2 0.2 0.2 0.2 0.2]';
  maxiter=10;
  
  %% initiate GCMMA
  gcmma=GCMMAFunctions.init(xval,xnam,xmin,xmax,maxiter);

end

%% run GCMMA
[gcmma,xval]=GCMMAFunctions.run(gcmma,xval,xnam,xmin,xmax,true);
[f0val,fval]=LEnOpFunctions.optFunctions(xval,xnam,false);

%% plot convergence
figure(1), clf, GCMMAFunctions.plotIter(gcmma)

% %% check results
% mass=LEnOpFunctions.computeMass(model)
% LCE=LEnOpFunctions.computeLCE(model,material)
% beam=computeEulerBernoulli(model);
% fval=max(abs(beam.w))
% comsol=runCOMSOLBeam(model);
% v=mpheval(comsol,'v','edim',1,'dataset','dset1');
% figure(2), clf, plot(beam.x,beam.w,v.p(1,:),v.d1,'*')
