% Determine if the 'comsol' field in the umod struct is a Comsol 4.x 
% object. 
%

function isit = iscmp4x(comsol)

    isit = isjava(comsol);