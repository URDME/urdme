function params = hes1red_params(type,par)
%HES1RED_PARAMS Parameters for reduced Hes1@2 models.
%   PARAMS = HES1RED_PARAMS(TYPE,PAR) returns parameters for reduced
%   type TYPE, assuming the original parameters are PAR, see
%   HES1_PARAMS.
%
%   For the reductions, see the manuscript cited in README.
%
%   Reduction alternatives are as follows:
%     alternative 1: (x,y) = (M,n) [type 1]
%     alternative 2: (x,y) = (P,n) [type 1]
%     alternative 3: (x,y) = (M,N) [type 2]
%     alternative 4: (x,y) = (P,N) [type 2]
%     alternative 5: (x,y) = (M,D) [type 1]
%     alternative 6: (x,y) = (P,D) [type 1]
%     alternative 7: (x,y) = (M,P) [type 3]
% 
%   The output is a struct with the following fields:
%   params.x0 - non-dimensionalising scaling for state x
%               for alternative 1 this is: x = M/M0,
%               for alternative 2 this is: x = P/P0 etc
%               so, for alternative 1: params.x0 = M0
%                   for alternative 2: params.x0 = P0 etc
%   params.y0 - non-dimensionalising scaling for state y
%               for alternative 1 this is: y = n/n0,
%               for alternative 3 this is: y = N/N0 etc
%               so, for alternative 1: scaling.y0 = n0
%                   for alternative 3: scaling.y0 = N0 etc
%   params.tau - scaling of time for each alternative
%   params.a - final parameter 'a'
%   params.b - final parameter 'b'
%   params.v - final parameter 'v'
%
%   Examples:
%     % a single set of parameters
%     par = hes1_params;
%     parred = hes1red_params(1,par);
%
%     % multiple sets
%     NMC = 50;
%     pars = hes1_params(NMC);
%     parsred = cell(1,7);
%     for i = 1:7
%       parsred{i} = hes1red_params(i,pars);
%     end
%
%     % note that a and b are the same throughout:
%     da = 0; db = 0;
%     for i = 1:6
%       da = max(da,norm(parsred{i}.a-parsred{i+1}.a,inf));
%       db = max(db,norm(parsred{i}.b-parsred{i+1}.b,inf));
%     end
%
%  See also HES1_PARAMS, HES1_CONC.

% G. Menz 2024-04-18

% scale given parameters according to given alternative
params = l_hes1_scale(type,par);
params.k = par.k;
params.h = par.h;

end
%-----------------------------------------------------------------------
function scaling = l_hes1_scale(alternative,par)
%L_HES1_SCALE Scalings used for different alternatives of reduced model

% type 1
if alternative == 1
  c1 = (par.alphaD .* par.alphaN .* par.alphaM) ./ (par.muD .* par.muN);
  c2 = (par.alphaP./(par.muP .* par.KM)).^par.k;
  c3 = (par.alphaP./(par.muP .* par.Kn)).^par.h;
  scaling.x0 = (c1 .* par.alphan./(c2 .* par.muM .* par.mun)).^...
               (1./(par.k+1));
  scaling.y0 = par.alphan ./ par.mun;
  scaling.tau = 1 ./ par.muM;
  scaling.a = 1 ./ (c2 .* scaling.x0.^par.k);
  scaling.b = c3 .* scaling.x0.^par.h;
  scaling.v = par.mun ./ par.muM;
% type 1
elseif alternative == 2
  c1 = (par.alphaD .* par.alphaN .* par.alphaM .* par.alphaP) ./ ...
       (par.muD .* par.muN .* par.muM);
  c2 = par.KM.^(-par.k);
  c3 = par.Kn.^(-par.h);
  scaling.x0 = (c1 .* par.alphan./(c2 .* par.muP .* par.mun)).^ ...
               (1./(par.k+1));
  scaling.y0 = par.alphan ./ par.mun;
  scaling.tau = 1 ./ par.muP;
  scaling.a = 1 ./ (c2 .* scaling.x0.^par.k);
  scaling.b = c3 .* scaling.x0.^par.h;
  scaling.v = par.mun ./ par.muP;
% type 2
elseif alternative == 3
  c1 = (par.alphaD .* par.alphaN .* par.alphan) ./ (par.muD .* par.mun);
  c2 = (par.alphaP ./ (par.muP .* par.KM)).^par.k;
  c3 = (par.alphaP ./ (par.muP .* par.Kn)).^par.h;
  scaling.y0 = c1 ./ par.muN;
  scaling.x0 = (par.alphaM .* scaling.y0./(par.muM .* c2)).^ ...
               (1./(par.k+1));
  scaling.tau = 1 ./ par.muM;
  scaling.a = 1 ./ (c2 .* scaling.x0.^par.k);
  scaling.b = c3.*scaling.x0.^par.h;
  scaling.v = par.muN ./ par.muM;
% type 2
elseif alternative == 4
  c1 = par.alphaM .* par.alphaP ./ par.muM;
  c2 = par.KM.^(-par.k);
  c3 = par.alphaD .*  par.alphaN .* par.alphan ./ (par.muD .* par.mun);
  c4 = par.Kn.^(-par.h);
  scaling.y0 = c3 ./ par.muN;
  scaling.x0 = (c1 .* scaling.y0 ./ (c2 .* par.muP)).^ ...
               (1./(par.k+1));
  scaling.tau = 1 ./ par.muP;
  scaling.a = 1 ./ (c2 .* scaling.x0.^par.k);
  scaling.b = c4 .* scaling.x0.^par.h;
  scaling.v = par.muN ./ par.muP;
% type 1
elseif alternative == 5
  c1 = par.alphaN .* par.alphaM ./ par.muN;
  c2 = (par.alphaP ./ (par.muP .* par.KM)).^par.k;
  c3 = par.alphaD .* par.alphan ./ par.mun;
  c4 = (par.alphaP ./ (par.muP .* par.Kn)).^par.h;
  scaling.y0 = c3 ./ par.muD;
  scaling.x0 = (c1 .* scaling.y0 ./ (par.muM .* c2)).^ ...
               (1./(par.k+1));
  scaling.tau = 1 ./ par.muM;
  scaling.a = 1 ./ (c2 .* scaling.x0.^par.k);
  scaling.b = c4 .* scaling.x0.^par.h;
  scaling.v = par.muD ./ par.muM;
% type 1
elseif alternative == 6
  c1 = par.alphaN .* par.alphaM .* par.alphaP ./ (par.muN .* par.muM);
  c2 = par.KM.^(-par.k);
  c3 = par.alphaD .* par.alphan ./ par.mun;
  c4 = par.Kn.^(-par.h);
  scaling.y0 = c3 ./ par.muD;
  scaling.x0 = (c1 .* scaling.y0 ./ (c2 .* par.muP)).^(1./(par.k+1));
  scaling.tau = 1 ./ par.muP;
  scaling.a = 1 ./ (c2 .* scaling.x0.^par.k);
  scaling.b = c4 .* scaling.x0.^par.h;
  scaling.v = par.muD ./ par.muP;
% type 3
elseif alternative == 7
  c1 = par.alphaD .* par.alphaN .* par.alphaM .* par.alphan ./ ...
      (par.muD .* par.muN .* par.mun);
  c2 = par.Kn.^(-par.h);
  c3 = par.KM.^(-par.k);
  scaling.x0 = (c1 ./ (par.muM .* c3 .* (par.alphaP ./ par.muP).^ ...
                par.k)).^(1./(par.k+1));
  scaling.y0 = par.alphaP .* scaling.x0 ./ par.muP;
  scaling.tau = 1 ./ par.muM;
  scaling.a = 1 ./ (c3 .* scaling.y0.^par.k);
  scaling.b = c2 .* scaling.y0.^par.h;
  scaling.v = par.muP ./ par.muM;
else
  error('Choose an alternative between 1 and 7.')
end
end
%-----------------------------------------------------------------------
