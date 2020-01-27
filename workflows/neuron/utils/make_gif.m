function make_gif(tree, t_vec, V_out, filename)
%MAKE_GIF Creates a .gif-animation of an input tree.
%   MAKE_GIF(TREE,T_VEC,V_OUT,FILENAME) creates a .gif-animation
%   'filename.gif' of the tree TREE with dynamic potential
%   (T_VEC,V_OUT).
%
%   Note: TREES Toolbox needs to be included.
%
%   See also MOVIE2GIF.

% A. Senek 2017-05-31

len_for = length(t_vec) -1;
  
figure(1),
M = struct('cdata',{},'colormap',{});
for j = 1:len_for
  figure(1), clf,
  color = [-70; V_out(j,2:end)'; 40];
  plot_tree(tree, color, [], [], [], '-3l');
  str = (['Time [ms] = ' num2str(t_vec(j),2)]);
  text(-300,300,str);
  drawnow;
  M(j) = getframe;
end
movie2gif(M,{M(1:5:50).cdata},filename,'delaytime',0.05,'loopcount',inf); 
