function model = model_builder(nVoxels,L,D)
%MODEL_BUILDER Neuron model builder.
%   MODEL = MODEL_BUILDER(nVoxels,L,D) creates a straight cylindrical
%   neuron of dimensions length L [um] long and D [um] in
%   diameter. The defaults are the standard Rallpack 1 & 3 cable: L =
%   1 [mm] = 1000 [um], D = 1 [um].
%
%   MODEL = MODEL_BUILDER(nVoxels,TREE) loads a TREES 1.15 file,
%   discretizes it into nVoxels voxels and returns a neuron model
%   structure.
%
%   Field          Description
%   --------------------------------------------------------------------
%   dA             Adjacency matrix, see TREES for an explanation.
%   coord          Coordinates (X,Y,Z).
%   nVoxels        Number of voxels.
%   vol            Volume vector.
%   L              Voxel lengths.
%   surfacearea    Voxel surface areas.
%   D              Voxel diameters.

% S. Engblom 2019-12-02 (Revision, syntax)
% S. Engblom 2018-06-27 (Revision)
% A. Senek 2017-05-31

if nargin == 2 && isstruct(L)
  % tree-syntax
  tree = L;
  if exist('vol_tree','file')
    model = tree;
    model.nVoxels = length(model.X)-1;
    switch sign(nVoxels-model.nVoxels)
     case -1
      disp(['Number of voxels increased to allow for ' ...
            'discretization']);
     case 1
      disp('Reshaping model size to nVoxels');
      L_quotient = nVoxels./length(model.L);
      model = resample_tree(model, model.L(1)/L_quotient);
    end

    % converts to voxels (index 1 = root)
    model.vol = vol_tree(model);
    % first vol will be nonsensical from TREES
    model.vol = model.vol(2:end);
    model.L = len_tree(model);
    % first length will be nonsensical from TREES
    model.L = model.L(2:end);
    model.surfacearea = surf_tree(model);
    % first area will be nonsensical from TREES
    model.surfacearea = model.surfacearea(2:end);
    model.nVoxels= numel(model.L);
  else
    error('Error message: TREES package has to be included'); 
  end
else
  % cable syntax

  % defaults
  if nargin < 3
    D = 1; % [um]
    if nargin < 2
      L = 1000; % [um]
    end
  end
  model = struct;

  % create (linear) tree structure
  model.coord.X = linspace(0,L,nVoxels+1)';
  model.coord.Y = zeros(1+nVoxels,1);
  model.coord.Z = zeros(1+nVoxels,1);
  model.D = D*ones(1+nVoxels,1);

  % directed adjacency matrix, a zero matrix with ones on the -1
  % diagonal
  model.dA = spdiags(ones(nVoxels+1,1),-1,nVoxels+1,nVoxels+1);

  % convert to voxels (index 1 = root)
  model.L = diff(model.coord.X);
  model.surfacearea = pi*model.D(2:end).*model.L;
  model.nVoxels = nVoxels;
  model.vol = pi*(model.D(2:end)/2).^2.*model.L;
end
