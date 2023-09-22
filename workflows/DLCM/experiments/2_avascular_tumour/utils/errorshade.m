function p = errorshade(x,lo,hi,color,alpha)
%ERRORSHADE Shaded error bars along line plot.
%   ERRORSHADE(X,LO,HI) plots(X,[LO HI]) as a shaded area meant to
%   represent an error around the line plot (X,Y) with Y inside
%   [LO,HI].
%
%   ERRORSHADE(...,Color,Alpha) controls the color and transparency of
%   the shade, respectively, via the RGB-vector Color and scalar
%   Alpha.
%
%   In all syntaxes, the patch object itself is returned, i.e., P =
%   ERRORSHADE(...).
%
%   Example:
%     x = -pi:pi/25:pi;
%     y = tan(sin(x))-sin(tan(x));
%     y_ = y.*(1+0.1*(rand(size(y))-0.5))+0.25*randn(size(y));
%     figure, clf,
%     h = plot(x,y_,'r.-'); col = get(h,'Color');
%
%     % standard error formula
%     e = sqrt((0.1*abs(y)).^2/12+0.25^2);
%     errorshade(x,y_-e,y_+e,col);
%
%     % "truth"
%     hold on, plot(x,y,'k--');
%
%     See also PLOT, PATCH, ERRORBAR.      

% S. Engblom 2022-05-20

if nargin < 5
  alpha = 0.2;
  if nargin < 4
    color = [0 0 0];
  end
end

xx = [x fliplr(x)];
yy = [lo fliplr(hi)];
p_ = patch('xdata',xx,'ydata',yy, ...
           'FaceAlpha',alpha,'FaceColor',color, ...
           'LineStyle','none','HandleVisibility','off');
if nargout > 0, p = p_; end