function [ISI_vec, t_max_vec] = ISI_calc(V_out, t_vec, thresh)
%ISI_CALC ISI statistics for a firing neuron.
%   [ISI_vec, t_max_vec] = ISI_CALC(V_out, t_vec, thresh) Calculates
%   statistics, such as ISI, for a firing neuron. Where a V_out value
%   > thresh counts as a peak.
%
%   Input V_out can be a (m,n) sized matrix where the individual
%   compartements are column-vectors, and n is the number of
%   compartments.  t_vec has to be a m long vector where t_vec(x)
%   corresponds to the time of V_out(x,:).
%
%   ISI_CALC works by assuming all peaks start below, and then pass,
%   the 0 volt limit. It then removes this "negative" region and
%   calculates the max of each region that still exists.
%
%   ISI_vec(x,1) - contains mean ISI for compartment x
%   ISI_vec(y,2) - contains standard deviation for compartment y

%   Aleksandar Senek 2017-03-13

m = size(V_out,2);
n = size(V_out,1);
ISI_vec = zeros(m,1); 
for k = 1:m
  temp_V = V_out(:,k);
  pos_idx = find(temp_V > thresh);
  
  %   Find all values that belong to separate peaks
  a=diff(pos_idx);
  b=find([a' inf]>1);
  c=diff([0 b]); %    length of the sequences
  d=cumsum(c); %  endpoints of the sequences
  d_ = [0 d];
  t_max_vec = zeros(length(d),1);
  
  for i = 1:length(d)
    if size(d_) == 0
      break;
    end
    peak_s = pos_idx((d_(i) + 1): d_(i+1)); 
    %   Find peaks of potentials and respective times of occurence.
    max_index = find(V_out(:,k) == max(V_out(peak_s,k)));
    t_max_vec(i) = t_vec(max_index);
  end 
  ISI_vec(k,1) = mean(diff(t_max_vec));
  ISI_vec(k,2) = std(diff(t_max_vec));
  ISI_vec(k,3) = length(t_max_vec);
end

%   Might be necessary if thresh is high
%ISI_vec(find(isnan(ISI_vec(:,1))),:) = [];
%ISI_vec(find(ISI_vec(:,2) == 0),:) = [];
