% function to load material parameters into P
% 


function P = LoadMaterialParams(P)

if strcmp(P.beamMat,'diamond') || strcmp(P.beamMat,'dia')
    P.beamMat = 'diamond';
    
    % if statements below for backcompatibility
    
    % material properties
    if ~isfield(P,'E') || (isfield(P,'E') && P.E ~= 1050e9)
        P.E = 1050e9;                           % Young's modulus (Pa)
    end
    if ~isfield(P,'rho') || (isfield(P,'rho') && P.rho ~= 3500)
        P.rho = 3500;                           % density (kg/m3)
    end
    if ~isfield(P,'nu') || (isfield(P,'nu') && P.nu ~= 0.2)
        P.nu = 0.2;                             % Poisson ratio
    end
    
    % Anisotropic elasticity matrix elements (units: Pa)
    c11 = 1076e9;
    c12 = 125e9;
    c44 = 578e9;
    
    % index of refraction in material
    if ~isfield(P,'nbeam') || (isfield(P,'nbeam') && P.nbeam ~= 2.386)
        P.nbeam = 2.386;
    end
    
    % photoelastic constants
    if ~isfield(P,'p11') || (isfield(P,'p11') && P.p11 ~= -0.25)
        P.p11 = -0.25;
    end
    if ~isfield(P,'p12') || (isfield(P,'p12') && P.p12 ~= 0.043)
        P.p12 = 0.043;
    end
    if ~isfield(P,'p44') || (isfield(P,'p44') && P.p44 ~= -0.172)
        P.p44 = -0.172;
    end
    
    % strain susceptibilities - from Sohn & Meesala et al., arXiv (2017)
    if ~isfield(P,'dg') || (isfield(P,'dg') && P.dg ~= 1.3e15)
        P.dg = 1.3e15;    
    end
    if ~isfield(P,'fg') || (isfield(P,'fg') && P.fg ~= -250e12)
        P.fg = -250e12;
    end
    
elseif strcmp(P.beamMat,'silicon') || strcmp(P.beamMat,'Si')
    P.beamMat = 'silicon';
    % material properties
    P.E = 168.5e9;                          % Young's modulus (Pa)
    P.rho = 2330;                           % density (kg/m3)
    P.nu = 0.28;                            % Poisson ratio
    
    % Anisotropic elasticity matrix elements (units: Pa)
    c11 = 166e9;
    c12 = 64e9;
    c44 = 80e9;
    
    % index of refraction in material
    P.nbeam = 3.493;
    
    % photoelastic constants
    P.p11 = -0.094;
    P.p12 = 0.017;
    P.p44 = -0.051;
    
elseif strcmp(P.beamMat,'LN') || strcmp(P.beamMat,'lithium niobate')
    P.beamMat = 'LN';
    % material properties 
    % Anisotropic elasticity matrix elements (units:Pa)
    c11 = 200e9;
    c12 = 55e9;
    c13 = 70e9;
    c14 = 7.9e9;
    c33 = 240e9;
    c44 = 60e9;
    c66 = 72e9;
    % material properties
    P.E = 170e9;                          % Young's modulus (Pa)
    P.rho = 4650;                           % density (kg/m3)
    P.nu = 0.28;                            % Poisson ratio
    % index of refraction in material
    P.nbeam = 3.493;
end

% Anisotropic elasticity matrix (units: Pa) - COMSOL v4+ ordering: 
% 11, 12,22, 13,23,33, 14,24,34,44, 15,25,35,45,55, 16,26,36,46,56,66
% disable P.D if using isotropic material (defined by E, nu, rho)
if isfield(P,'anisoMat') && P.anisoMat
    if strcmp(P.beamMat,'LN')
        P.D = [c11, c12,c11, c13,c13,c33, c14,-c14,0,c44, 0,0,0,0,c44, 0,0,0,0,c14,c66];
    else
        P.D = [c11, c12,c11, c12,c12,c11, 0,0,0,c44, 0,0,0,0,c44, 0,0,0,0,0,c44];
    end
end





end