function [curvature] = curvature2D(line_coords_cellarr,desired_coords,hmax,cmin,win,deg)
%CURVATURE2D evaluates the signed curvature of a closed 2D line.
%   CURVATURE2D(line_coords,desired_coords) evaluates the curvature of the
%   line comprised of the cartesian coordinates in line_coords at the
%   coordinates in desired_coords. The line coordinates are run through a
%   Savitsky-Golay filter and the curvature is evaluated using splines and
%   the general analytical expression for parametrised curves in 2D.
%
%   in: line_coords     - cartesian coordinates of closed line in 2D
%                         (for cell array holding the set of closed lines)
%       desired_coords  - cartesian coordinates of output curvature
%   out: curvature      - curvature of line_coords at desired_coords
%
%   Example:
%       theta = -pi:0.1:pi;
%       k = 4;
%       line_coords = [(R + eps.*cos(k.*theta)).*cos(theta); ...
%                      (R + eps.*sin(k.*theta)).*sin(theta)];
%       desired_coords = [(R + eps.*cos(k.*(theta+0.05))).*cos(theta+0.05); ...
%                      (R + eps.*cos(k.*(theta+0.05))).*sin(theta+0.05)];
%       [curvature] = curvature2D(line_coords,desired_coords, 0.01, 0.0);
%       plot(desired_theta, curvature)
%
%       See also CONTOUR, SGOLAYFILT, SPLINE
%

% E. Blom 2023-02-20

% NOTE: sets surface tension to zero to too small contours
curvature_all = []; % arrays for each closed line to be appended
line_coords_all = [];

% find contour sizes
for i = 1:numel(line_coords_cellarr)
    contour_size(i) = length(line_coords_cellarr{i});
end
    
max_size = max(contour_size);
min_size_factor = cmin;  % disregard contours less than this fraction of max_size

% evaluate curvature for each closed line
for i = 1:numel(line_coords_cellarr)
line_coords = line_coords_cellarr{i};

line_len = length(line_coords);
if line_len >= max_size*min_size_factor
    % Savitzky-Golay filter (hereon: sgf) properties
    sgf_wind = min(length(line_coords),win);  % window size, was 21, PDE (51, DLCM)
    sgf_deg = min(sgf_wind-2,deg);    % polynomial degree, was 7, PDE (3, DLCM)

    % sgf on line coordinates (includes overlapping for periodicity)
    lx_sgf = sgolayfilt([line_coords(1,end-sgf_wind+1:end) line_coords(1,:) line_coords(1,1:sgf_wind)], ...
                        sgf_deg,sgf_wind);
    ly_sgf = sgolayfilt([line_coords(2,end-sgf_wind+1:end) line_coords(2,:) line_coords(2,1:sgf_wind)], ...
                        sgf_deg,sgf_wind);

    % find curvature by parametric representation x(t), y(t) 
    % NOTE: highly and suddenly inaccurate for less than ~ 100 points
    x = lx_sgf;
    y = ly_sgf;
    x0 = spline(1:numel(x), x);
    y0 = spline(1:numel(y), y);
    % 1st, 2nd derivative of x(t)
    x1 = mkpp(x0.breaks, [3.*x0.coefs(:,1); 2.*x0.coefs(:,2); x0.coefs(:,3)]);
    x2 = mkpp(x0.breaks, [3.*x0.coefs(:,1); 2.*x0.coefs(:,2)]);
    % 1st, 2nd derivative of y(t)
    y1 = mkpp(y0.breaks, [3.*y0.coefs(:,1); 2.*y0.coefs(:,2); y0.coefs(:,3)]);
    y2 = mkpp(y0.breaks, [3.*y0.coefs(:,1); 2.*y0.coefs(:,2)]);
    curvature = (ppval(x1,x1.breaks).*ppval(y2,y2.breaks) ... 
        - ppval(y1,y1.breaks).*ppval(x2,x2.breaks))./ ...
        (ppval(x1,x1.breaks).^2 + ppval(y1,y1.breaks).^2).^1.5; 

    % filter and remove overlap
    curvature = sgolayfilt(curvature, sgf_deg, sgf_wind);
    curvature = curvature(sgf_wind+1:end-sgf_wind); % errors here can induce rotational bias

    % Attempt to speed up process by limiting max curvature by resolution
    curvature = min(curvature, 1/hmax);
    curvature = max(curvature, -1/hmax);

else
    % disregard small contours (but set to NaN to remember this)
    curvature = NaN.*zeros(1,line_len);
end

% NOTE: could optimize below, knowing final lengths
curvature_all = [curvature_all curvature];
line_coords_all = [line_coords_all line_coords]; 
end


% Find value at coordinate closest to desired one (interpolate by
% averaging if two or more are closest)
for i = 1:numel(desired_coords(1,:))
    idx = find((desired_coords(1,i) - line_coords_all(1,:)).^2 +  ...
                  (desired_coords(2,i) - line_coords_all(2,:)).^2 == ...
              min((desired_coords(1,i) - line_coords_all(1,:)).^2 +  ...
                  (desired_coords(2,i) - line_coords_all(2,:)).^2));
    temp_curv_all = curvature_all(idx); % find the values of the closest pts...
    temp_curv_all = temp_curv_all(~isnan(temp_curv_all)); % ...remove NaN...
    %if isnan(temp_curv_all)
    %    temp_curv_all = 0; % crude fix to avoid 'false' NaN:s
    %end
    desired_curvature(i) = mean(temp_curv_all); % ... take average
end

curvature = desired_curvature;  
end

