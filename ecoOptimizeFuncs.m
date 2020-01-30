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

classdef ecoOptimizeFuncs
  methods(Static)
    
    %%
    function material=blendMaterials(model,materialsData)
      N=max(sum(not(cellfun(@isempty,model.material)),1));
      if N>1
        model.alpha(N,:)=1-sum(model.alpha(1:N-1,:),1);
      end
      fnames=fieldnames(materialsData);
      fnames=fnames(not(ismember(fnames,'info')));
      for j=1:numel(fnames)
        eval(['material.',fnames{j},'=zeros(size(model.alpha));'])
        for i=1:numel(model.alpha)
          if not(cellfun(@isempty,model.material(i)))
            eval(['material.',fnames{j},'(i)=materialsData(ismember([materialsData.info],model.material(i))).',fnames{j},'{1};'])
          end
        end
        eval(['material.',fnames{j},'=sum(material.',fnames{j},'.*model.alpha,1);'])
      end
    end

    %% 
    function model=updateMaterialProps(model,material)
      model.Ei=material.youngsModulus;
      model.nui=material.poissonsRatio;
      model.rhoi=material.density;
      model.EProdi=material.productionEnergy;
      model.EDispi=material.disposalEnergy;
      model.EEoLi=material.eolEnergy;
      model.CO2Prodi=material.productionCO2;
      model.CO2Dispi=material.diposalCO2;
      model.CO2EoLi=material.eolCO2;
      model.Costi=material.price;
    end
    
    %%
    function model=updateDependentVars(model)
    % this functions must be update for different model paramaterizations
      switch model.xsection
        case 'layered'
          I0=model.B.*model.H.^3/12;
          d=([0 cumsum(model.H(1:end-1))]+model.H/2-sum(model.H)/2);
          A=model.B.*model.H;
          I=I0+A.*d.^2;
          model.Ai=A;
          model.Ii=I;
          model.E=sum(model.Ei.*I)/sum(I);
          model.nu=sum(model.nui.*A)/sum(A);!
          model.rho=sum(model.rhoi.*A)/sum(A);
          model.A=sum(A);
          model.I=sum(I);
          model.EProd=sum(model.EProdi.*A)/sum(A);
          model.EDisp=sum(model.EDispi.*A)/sum(A);
          model.EEoL=sum(model.EEoLi.*A)/sum(A);
          model.CO2Prod=sum(model.CO2Prodi.*A)/sum(A);
          model.CO2Disp=sum(model.CO2Dispi.*A)/sum(A);
          model.CO2EoL=sum(model.CO2EoLi.*A)/sum(A);
          model.Cost=sum(model.Costi.*A)/sum(A);
      end
    end
    
    %%
    function mass=computeMass(model)
      mass=model.A*model.L*model.rho;
    end
    
    %%
    function LCE=computeLCE(model)
      mass=ecoOptimizeFuncs.computeMass(model);
      % production phase
      Ep_kg=model.EProd; %[J/kg]
      Ep=Ep_kg*mass;
      % use phase
      usemodel='physicsbased';
      switch usemodel
        case 'simple'
          %simple model based on fuel efficiency
          Eu_km_kg=10.8/100*33.7e6/1400; %[J/km/kg]
          Eu=Eu_km_kg*model.driveDistTotal*mass; %[J]
        case 'physicsbased'
          %more advanced model based on drive cycle energy (O'Reilly (2016))
          crr=0.01*0.85; %[-] %coefficient of effective rolling resistance
          CR=11013; %[m] drive cycle rolling resistance constant
          CA=1227; %[m^2/s^2] drive cycle acceleration constant
          diffEff=0.42; %[-] differential efficiency (petrol)
          Eu_km_kg=(9.81*crr*CR+CA)/CR*1e3/diffEff; %[J/km/kg]
          Eu=Eu_km_kg*model.driveDistTotal*mass; %[J]
      end
      % disposal phase
      Ed_kg=model.EDisp; %[J/kg]
      Ed=Ed_kg*mass; %[J]
      % end-of-life phase
      Ee_kg=model.EEoL; %[J/kg]
      Ee=Ee_kg*mass; %[J]
      % total
      LCE=[Ep Eu Ed Ee]*1e-9;
    end
    
    %%
    function LCCO2=computeLCCO2(model)
      mass=ecoOptimizeFuncs.computeMass(model);
      % production phase
      CO2p_kg=model.CO2Prod; %[kg/kg]
      CO2p=CO2p_kg*mass; %[kg]
      % use phase
      usemodel='physicsbased';
      switch usemodel
        case 'simple'
          %simple model based on fuel efficiency
          CO2u_km_kg=10.8/100*2.31/1400; %[kg/km/kg]
          CO2u=CO2u_km_kg*model.driveDistTotal*mass; %[kg]
        case 'physicsbased'
          %more advanced model based on drive cycle energy
          crr=0.01*0.85; %[-] %coefficient of effective rolling resistance
          CR=11013; %[m] drive cycle rolling resistance constant
          CA=1227; %[m^2/s^2] drive cycle acceleration constant
          diffEff=0.42; %[-] differential efficiency (petrol)
          Eu_km_kg=(9.81*crr*CR+CA)/CR*1e3/diffEff; %[J/km/kg]
          heatValueFuel=43.5*1e6*0.75; %[J/L]
          CO2Fuel=2.31; %[kg/L]
          fuelUse_km_kg=Eu_km_kg/heatValueFuel; %[L/km/kg]
          CO2u_km_kg=fuelUse_km_kg*CO2Fuel; %[kg/km/kg]
          CO2u=CO2u_km_kg*model.driveDistTotal*mass; %[kg]
      end
      % disposal phase
      CO2d_kg=model.CO2Disp; %[kg/kg]
      CO2d=CO2d_kg*mass; %[kg]
      % end-of-life phase
      CO2e_kg=model.CO2EoL; %[kg/kg]
      CO2e=CO2e_kg*mass; %[kg]
      % total
      LCCO2=[CO2p CO2u CO2d CO2e]*1e-3;
    end
    
    %%
    function LCCost=computeLCCost(model)
      mass=ecoOptimizeFuncs.computeMass(model);
      % production phase
      Costp_kg=model.Cost; %[SEK/kg]
      Costp=Costp_kg*mass; %[SEK]
      % use phase
      usemodel='physicsbased';
      switch usemodel
        case 'simple'
          %simple model based on fuel efficiency
          Costu_km_kg=10.8/100*15/1400; %[SEK/km/kg]
          Costu=Costu_km_kg*model.driveDistTotal*mass; %[J]
        case 'physicsbased'
          %more advanced model based on drive cycle energy
          crr=0.01*0.85; %[-] %coefficient of effective rolling resistance
          CR=11013; %[m] drive cycle rolling resistance constant
          CA=1227; %[m^2/s^2] drive cycle acceleration constant
          diffEff=0.42; %[-] differential efficiency (petrol)
          Eu_km_kg=(9.81*crr*CR+CA)/CR*1e3/diffEff; %[J/km/kg]
          heatValueFuel=43.5*1e6*0.75; %[J/L]
          fuelUse_km_kg=Eu_km_kg/heatValueFuel; %[L/km/kg]       
          priceFuel=15; %[SEK/L]
          Costu_km_kg=fuelUse_km_kg*priceFuel; %[SEK/km/kg]
          Costu=Costu_km_kg*model.driveDistTotal*mass; %[SEK]  
      end
      % disposal phase
      Costd=0; %[SEK]
      % end-of-life phase
      Coste=0; %[SEK]
      % total
      LCCost=[Costp Costu Costd Coste]*1e-3;
    end
    
    %%
    function fval=computeConstraints(model,update)
      method='analytical';
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
          disp('Computional method for constraints not defined')
      end
    end
    
    %%
    function [f0val,fval,df0dx,dfdx] = optFuncs(x,xnam,grads)
      %  This calculates function values and gradients
      global model material materialsData
      fmax=[1e-3]';% -60]';
      scale=[1e3]';% 1]';
      for i=1:numel(x)
        eval(['model.',xnam{i},'=x(i);'])
      end
      material=ecoOptimizeFuncs.blendMaterials(model,materialsData);
      model=ecoOptimizeFuncs.updateMaterialProps(model,material);
      model=ecoOptimizeFuncs.updateDependentVars(model);
      eval(['f0val = sum(ecoOptimizeFuncs.compute',model.objfunc,'(model));'])
      fval  = scale.*[ecoOptimizeFuncs.computeConstraints(model,1)'-fmax];
      if grads
        dx=1e-4;
        df0dx=[];
        dfdx=[];
        for i=1:numel(x)
          eval(['model.',xnam{i},'=x(i)+dx;'])
          material=ecoOptimizeFuncs.blendMaterials(model,materialsData);
          model=ecoOptimizeFuncs.updateMaterialProps(model,material);
          model=ecoOptimizeFuncs.updateDependentVars(model);
          eval(['df0dx = [df0dx; (sum(ecoOptimizeFuncs.compute',model.objfunc,'(model))-f0val)/dx];'])
          dfdx  = [dfdx, scale.*[ecoOptimizeFuncs.computeConstraints(model,1)'-(fval./scale+fmax)]/dx];
          eval(['model.',xnam{i},'=x(i);'])
          material=ecoOptimizeFuncs.blendMaterials(model,materialsData);
          model=ecoOptimizeFuncs.updateMaterialProps(model,material);
          model=ecoOptimizeFuncs.updateDependentVars(model);
        end
      end
    end
    
    %%
    function dispModel(model)
      N=max([numel(model.B) numel(model.H)]);
      B=repmat(model.B,1,N-numel(model.B)+1);
      H=repmat(model.H,1,N-numel(model.H)+1);
      H0=([0 cumsum(H(1:end-1))]-sum(H)/2);
      B0=-B/2;
      C=colormap;
      for i=1:N
        R=rectangle('Position',[B0(i) H0(i) B(i) H(i)]);
        set(R,'Facecolor',C(50*(i-1)+1,:))
      end
      axis([B0(i)-B(i)/10 B0(i)+B(i)+B(i)/10 H0(i)-H(i)/10 H0(i)+H(i)+H(i)/10])
      axis auto, axis equal
      xlabel('b [m]'), ylabel('h [m]')
    end
    
  end
end