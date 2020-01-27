function init_u0 = set_u0(init_V,num_ions,prop_relpath)
%SET_U0
%   u0 = SET_U0(init_V) Finds the propensities of the membrane
%   reactions and calculates the expected distribution at init_V mV
%   initial potential.

% A. Senek 2017-05-31

if nargin < 3
    prop_relpath = '';
end
mexmake_uds('channel/lib/HHKNa_prop.c', ...
            'source',{{'channel/lib/HHKfun.c' 'channel/lib/HHNafun.c'}});

t = [0 0.1];
u0(1:13) = 1;
u0 = u0(:);
N = 28;
vol = 1;
ldata = init_V;
gdata = [];
sd = 1;
K = zeros(3,0);
I = zeros(3,0);
S = sparse(0,0);
%   Get the transition probabilities for each reaction
prop_vec = mexrhs(t, u0, N, vol, ldata, gdata, sd, K, I, S);

R1a = prop_vec(1:2:8);
R1b = prop_vec(2:2:8);

R2 = prop_vec(9:end);

%   Create transition matrices
D1 = diag(R1b,-1);
D1 = D1 + diag(R1a,1);

D1 = D1 + diag(-sum(D1,2));

D2 = matrix_builder(R2);
D2 = D2 + diag(-sum(D2,2));

% Find null space to get ratio of channel states
L1 = abs(null(D1'));
L1 = L1./sum(L1(:));
% $$$ L1 = (num_ions(:,1) .* L1')';
L1 = tprod(num_ions(:,1),L1,[2 3],[1 3]);

L2 = abs(null(D2'));
L2 = L2./sum(L2(:));
% $$$ L2 = (num_ions(:,2) .* L2')';
L2 = tprod(num_ions(:,2),L2,[2 3],[1 3]);

init_u0 = round([L1; ...
                 L2]); 
end

function [D2] = matrix_builder(R2)
%   Kinetic scheme of Sodium:
%
%   S1 ⟺ S3 ⟺ S5 ⟺ S7
%   ⟰     ⟰     ⟰     ⟰
%   ⟱     ⟱     ⟱     ⟱
%   S2 ⟺ S4 ⟺ S6 ⟺ S8

N_states = 8; 

D2 = zeros(N_states);

%   S1 -> S2 and S3
D2(1,[2 3]) = R2(1:2);
%   S2 -> S1 and S4 :
D2(2,[1 4]) = R2(3:4);

for i = 1:N_states-4
    %   Even numbers in the rection-system:
    %   S(i = even) -> S(i-2), S(i-1), S(i+2)
    %   Odd: S(i = odd) -> S(i-2), S(i+1), S(i+2)
    sign = mod(i,2);
    
    D2(2 + i, [i, (i+ 2*sign + 1), (4+i)]) = R2((2+3*i):(4+3*i));
end

%   S7 -> S5 and S8 :
D2(7,[5 8]) = R2(17:18);
%   S8 -> S6 and S7
D2(8,[6 7]) = R2(19:20);
end