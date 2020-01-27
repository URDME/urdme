function pre = findPredecessor(ward,dA)
% FINDPREDECESSOR utility for the NEURON2COMSOL function.

% A. Senek 2017-07-11 (Adaptation)

pre = find(dA(ward,:) == 1);
if ward == 1
  pre = 1;
end
