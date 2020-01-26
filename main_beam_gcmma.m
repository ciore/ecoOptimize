% main_beam_gcmma.m

restart = 1;

if restart
  
  clear all
  clear global
  
  %% select a material
  addpath('.')
  materialsData=load('materialData.mat');
  materialsData=materialsData.data;
  global material
  material=materialsData(1);
  
  %% initiate model of the panel
  global model
  model=initModel;
%   model.E=material.youngsModulus{1};
%   model.rho=material.density{1};
%   model.nu=material.poissonsRatio{1};
  
  %% set optimisation params
  xval=[model.B(1) model.B(2) model.B(3) model.H(1) model.H(2) model.H(3)]';
  xnam={'B(1)' 'B(2)' 'B(3)' 'H(1)' 'H(2)' 'H(3)'};
  xmin=[0.05 0.01 0.05 0.01 0.01 0.01]';
  xmax=[0.2 0.2 0.2 0.2 0.2 0.2]';
  maxiter=10;
  
  %% initiate GCMMA
  gcmma=initGCMMA(xval,xnam,xmin,xmax,maxiter);

end

%% run GCMMA
[gcmma,xval]=runGCMMA(gcmma,xval,xnam,xmin,xmax,true);
[f0val,fval]=optFunctions(xval,xnam,false);

%% plot convergence
figure(1), clf, plotIter(gcmma)

%% check results
mass=computeMass(model)
LCE=computeLCE(model,material)
beam=computeEulerBernoulli(model);
fval=max(abs(beam.w))
comsol=runCOMSOLBeam(model);
v=mpheval(comsol,'v','edim',1,'dataset','dset1');
figure(2), clf, plot(beam.x,beam.w,v.p(1,:),v.d1,'*')





%% FUNCTIONS

%%
function model=initModel
  global model
  model.F=-1e4;
  model.xF=0.5;
  model.L=2;  
  model.B=[0.1 0.1 0.1];
  model.H=[0.01 0.08 0.01];
  model.El=[2.1e11 2.1e11 2.1e11];
  model.nul=[0.285 0.285 0.285];
  model.rhol=[7.8e3 7.8e3 7.8e3];
  model.xsection='layered';
  model.loadcase='simple_pt';
  model=updateDependentVars(model);
end

%% 
function model=updateDependentVars(model)
  switch model.xsection
    case 'layered'
      I0=model.B.*model.H.^3/12;
      d=([0 cumsum(model.H(1:end-1))]+model.H/2-sum(model.H)/2);
      A=model.B.*model.H;
      I=I0+A.*d.^2;
      model.E=sum(model.El.*I)/sum(I);
      model.nu=sum(model.nul.*A)/sum(A);      
      model.rho=sum(model.rhol.*A)/sum(A);
      model.A=sum(A);
      model.I=sum(I);
  end
end

%% 
function mass=computeMass(model)
  mass=model.A*model.L*model.rho;
end

%% 
function LCE=computeLCE(model,material)
  mass=computeMass(model);
  %% production phase
  Ep_kg=material.productionEnergy{1}; %[J/kg]
  Ep=Ep_kg*mass;
  %% use phase
  driveDistTotal=1e5; %[km] %total driving distance
  usemodel='simple';
  switch usemodel
    case 'simple'
      %simple model based on fuel efficiency
      Eu_km_kg=10.8/100*33.7e6/1400; %[J/km/kg]
      Eu=Eu_km_kg*driveDistTotal*mass; %[J]
    case 'physicsbased'
      %more advanced model based on drive cycle energy (Koffler (2010))
      crr=0.01*0.85; %[-] %coefficient of effective rolling resistance
      CR=11013; %[m] drive cycle rolling resistance constant
      CA=1227; %[m^2/s^2] drive cycle acceleration constant
      diffEff=0.42; %[-] differential efficiency
      Eu_km_kg=(9.81*crr*CR+CA)/CR*1e3/diffEff; %[J/km/kg]
      Eu=Eu_km_kg*driveDistTotal*mass; %[J]
  end
  %% end-of-life phase
  Ee_kg=material.eoLPotentialEnergy{1}; %[J/kg]
  Ee=Ee_kg*mass; %[J]
  %% total
  LCE=1e-9*(Ep+Eu+Ee);
end

%%
function fval=computeConstraints(model,update)
  addpath('../beamEB')
  method='analytical';
  model=updateDependentVars(model);
  switch method
    case 'analytical'
      beam=computeEulerBernoulli(model);
      fval(1)=max(abs(beam.w));
    case 'comsol'  
      if nargin<2
        update=0;
      end
      if model.H<1e-6
        fval=1e20;
      else
        if update
          comsol=runCOMSOLBeam(model,1);
        else
          comsol=runCOMSOLBeam(model,0);
        end
        fval(1)=max(abs(mpheval(comsol,'v','edim',1,'dataset','dset1','dataonly','on')));
    %     v=mpheval(comsol,'freq','edim',1,'dataset','dset2');
    %     fval(2)=-v.d1(1,1);
      end
    otherwise 
      print('Computional method for constraints not defined')
  end
end

%%
function [f0val,fval,df0dx,dfdx] = optFunctions(x,xnam,grads)
%  This file calculates function values and gradients
  global model material
  fmax=[1e-3]';% -60]';
  scale=[1e3]';% 1]';
  for i=1:numel(x)
    eval(['model.',xnam{i},'=x(i);'])
  end
  model=updateDependentVars(model);
  f0val = computeLCE(model,material);
  fval  = scale.*[computeConstraints(model,1)'-fmax];
  if grads
    dx=1e-4;
    df0dx=[];
    dfdx=[];
    for i=1:numel(x)
      eval(['model.',xnam{i},'=x(i)+dx;'])
      model=updateDependentVars(model);
      df0dx = [df0dx; (computeLCE(model,material)-f0val)/dx]; 
      dfdx  = [dfdx, scale.*[computeConstraints(model,1)'-(fval./scale+fmax)]/dx];
      eval(['model.',xnam{i},'=x(i);'])
      model=updateDependentVars(model);
    end
  end
end

%%
function gcmma = initGCMMA(xval,xnam,xmin,xmax,maxiter)
% parameters and the starting point are defined
  [f0val,fval,df0dx,dfdx] = optFunctions(xval,xnam,true);
  m=size(fval,1);
  n=size(xval,1);
  epsimin = 0.0000001;
  xold1   = xval;
  xold2   = xval;
  low     = xmin;
  upp     = xmax;
  c       = 1000*ones(m,1);
  d       = 1*ones(m,1);
  a0      = 1;
  a       = 0*ones(m,1);
  raa0    = 0.01;
  raa     = 0.01*ones(m,1);
  raa0eps = 0.000001;
  raaeps  = 0.000001*ones(m,1);
  kkttol  = 0;
  innerit = 0;
  outeriter = 0;
  maxoutit  = maxiter;
  iout = [outeriter innerit];
  xout = [xval'];
  fout = [f0val];
  cout = [fval'];
  for name=who'
    if not(ismember(name,{'xval','xnam','xmin','xmax'}))
      eval(['gcmma.',name{:},'=',name{:},';'])
    end
  end
end

%%
function [gcmma,xval] = runGCMMA(gcmma,xval,xnam,xmin,xmax,plotlive)
  % This is a slightly modified version of GCMMA-MMA-code-1.5/gctoymain.m
  %
  addpath('../GCMMA-MMA-code-1.5')
  % 
  names=fieldnames(gcmma);
  for i=1:numel(names)
    eval([names{i},'=gcmma.',names{i},';'])
  end
  %%%% begin optimisation
  %%%% If outeriter=0, the user should now calculate function values
  %%%% and gradients of the objective- and constraint functions at xval.
  %%%% The results should be put in f0val, df0dx, fval and dfdx:
  if outeriter < 0.5
    [f0val,fval,df0dx,dfdx] = optFunctions(xval,xnam,true);
    innerit=0;
    iout = [outeriter innerit];
    xout = [xval'];
    fout = [f0val];
    cout = [fval'];
  end
  %
  %%%% The outer iterations start:
  kktnorm = kkttol+10;
  outit = 0;
  while kktnorm > kkttol & outit < maxoutit
    outit   = outit+1;
    outeriter = outeriter+1;
  %%%% The parameters low, upp, raa0 and raa are calculated:
    [low,upp,raa0,raa] = ...
    asymp(outeriter,n,xval,xold1,xold2,xmin,xmax,low,upp, ...
          raa0,raa,raa0eps,raaeps,df0dx,dfdx);
  %%%% The GCMMA subproblem is solved at the point xval:
    [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,f0app,fapp] = ...
    gcmmasub(m,n,outeriter,epsimin,xval,xmin,xmax,low,upp, ...
             raa0,raa,f0val,df0dx,fval,dfdx,a0,a,c,d);
  %%%% The user should now calculate function values (no gradients)
  %%%% of the objective- and constraint functions at the point xmma
  %%%% ( = the optimal solution of the subproblem).
  %%%% The results should be put in f0valnew and fvalnew.
    [f0valnew,fvalnew] = optFunctions(xmma,xnam,false);
  %%%% It is checked if the approximations are conservative:
    [conserv] = concheck(m,epsimin,f0app,f0valnew,fapp,fvalnew);
  %%%% While the approximations are non-conservative (conserv=0),
  %%%% repeated inner iterations are made:
    innerit=0;
    if conserv == 0
      while conserv == 0 & innerit <= 15
        innerit = innerit+1;
  %%%% New values on the parameters raa0 and raa are calculated:
        [raa0,raa] = ...
        raaupdate(xmma,xval,xmin,xmax,low,upp,f0valnew,fvalnew, ...
                  f0app,fapp,raa0,raa,raa0eps,raaeps,epsimin);
  %%%% The GCMMA subproblem is solved with these new raa0 and raa:
        [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,f0app,fapp] = ...
        gcmmasub(m,n,outeriter,epsimin,xval,xmin,xmax,low,upp, ...
                 raa0,raa,f0val,df0dx,fval,dfdx,a0,a,c,d);
  %%%% The user should now calculate function values (no gradients)
  %%%% of the objective- and constraint functions at the point xmma
  %%%% ( = the optimal solution of the subproblem).
  %%%% The results should be put in f0valnew and fvalnew:
        [f0valnew,fvalnew] = optFunctions(xmma,xnam,false);
  %%%% It is checked if the approximations have become conservative:
        [conserv] = concheck(m,epsimin,f0app,f0valnew,fapp,fvalnew);
      end
    end
  %%%% No more inner iterations. Some vectors are updated:
    xold2 = xold1;
    xold1 = xval;
    xval  = xmma;
  %%%% The user should now calculate function values and gradients
  %%%% of the objective- and constraint functions at xval.
  %%%% The results should be put in f0val, df0dx, fval and dfdx:
    [f0val,fval,df0dx,dfdx] = optFunctions(xval,xnam,true);
  %%%% The residual vector of the KKT conditions is calculated:
    [residu,kktnorm,residumax] = ...
    kktcheck(m,n,xmma,ymma,zmma,lam,xsi,eta,mu,zet,s, ...
             xmin,xmax,df0dx,fval,dfdx,a0,a,c,d);
    iout = [iout; outeriter innerit];
    xout = [xout; xval'];
    fout = [fout; f0val];
    cout = [cout; fval'];
    for i=1:numel(names)
      eval(['gcmma.',names{i},'=',names{i},';'])
    end
    dispIter(gcmma)
    if plotlive
      plotIter(gcmma),drawnow
    end
  end
end

%%
function plotIter(gcmma)
  subplot(3,2,1)
  plot(gcmma.iout(:,1),gcmma.fout)
  xlabel('iter'),ylabel('f0val')
  subplot(3,2,2)
  semilogy(gcmma.iout(2:end,1),abs(gcmma.fout(2:end,:)-gcmma.fout(1:end-1,:))./abs(gcmma.fout(1,:)))
  xlabel('iter'),ylabel('df0val/di')
  subplot(3,2,3)
  plot(gcmma.iout(:,1),gcmma.xout)
  xlabel('iter'),ylabel('xval')
  subplot(3,2,4)
  semilogy(gcmma.iout(2:end,1),abs(gcmma.xout(2:end,:)-gcmma.xout(1:end-1,:))./abs(gcmma.xout(1,:)))
  xlabel('iter'),ylabel('dxval/di')
  subplot(3,2,5)
  plot(gcmma.iout(:,1),gcmma.cout)
  xlabel('outit'),ylabel('fval')
  subplot(3,2,6)
  semilogy(gcmma.iout(2:end,1),abs(gcmma.cout(2:end,:)-gcmma.cout(1:end-1,:))./abs(gcmma.cout(1,:)))
  xlabel('iter'),ylabel('dfval/di')
end

%%
function dispIter(gcmma)
  if gcmma.iout(end,1)==1
    fprintf(['GCMMA iterations:\n'])
    fprintf(['   iter   nini   f0val',repmat(' ',1,12-5),'xval',repmat(' ',1,gcmma.n*12-4),'fval\n'])
  end
  format=['%7i%7i',repmat('%12.3e',1,1+gcmma.n+gcmma.m),'\n'];
  fprintf(format,[gcmma.iout(end,:) gcmma.fout(end,:) gcmma.xout(end,:) gcmma.cout(end,:)])
end

