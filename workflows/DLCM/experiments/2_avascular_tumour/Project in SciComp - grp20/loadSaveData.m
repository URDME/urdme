% Load save data

U = saveData.U;
Usave = saveData.Usave;
tspan = saveData.tspan;
R = saveData.R;
V = saveData.V;
BC1 = saveData.BC1;
BC2 = saveData.BC2;
if exist('saveData.doRadiusBC','var') == 1
    doRadiusBC = saveData.doRadiusBC;
else
    doRadiusBC = false;
end