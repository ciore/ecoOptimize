%% Beam example of eco design
% This code is written with the intention of developing the authors
% understanding of how the code is structured.
%
%
% Author: Robert Jonsson


%% Clear 
clc
clear all
close all

%% Add material data to path
addpath('.') %Path to material database

%% Add GCMMA algorithm
addpath('/Users/robertkth/Documents/GitHub/GCMMA-MMA-code-1.5')

%% Add constraint solver
addpath('/Users/robertkth/Documents/GitHub/beamEB')

%% Import material data
materialsData = importdata('materialData.mat');

%% Import ecoOptimize functions
import ecoOptimize.*

%% Iterating over materials

for m = 1:9
    
    % Setting up model geometry
    model = initmodel;
    
    % Adding the model material properties
    model = setModelMaterial(model,materialsData(m));
    
    % Extracting the material names
    names(m) = materialsData(m).info;
    
    
    
end


%% Functions used

function model = initmodel
% Create and initialize/restore the model for each material.

    model.fmax = [1e-3];
    model.driveDistTotal = 1e5;
    model.P = -1e4;
    model.L = 1;
    model.B = 0.15*ones(1,3);
    model.H = 0.05*ones(1,3);
    
end


function model = setModelMaterial(model,materialsData)
% Important thing to remember here is that the materialsData input is not
% the full struct. The input is every "column" from the materialsData
% struct.

    model.E=materialsData.youngsModulus{1};
    model.rho=materialsData.density{1};
    model.EProd=materialsData.productionEnergy{1};
    model.EDisp=materialsData.disposalEnergy{1};
    model.EEoL=materialsData.eolEnergy{1};
    model.CO2Prod=materialsData.productionCO2{1};
    model.CO2Disp=materialsData.disposalCO2{1};
    model.CO2EoL=materialsData.eolCO2{1};
    model.Cost=materialsData.price{1};
end

function model = computeOptimalVariable(model)
% Moment of inertia
    I0=model.B.*model.H.^3/12;
    d=([0 cumsum(model.H(1:end-1))]+model.H/2-sum(model.H)/2);
    A=model.B.*model.H;
    I=I0+A.*d.^2;
    
    
    
    
end


