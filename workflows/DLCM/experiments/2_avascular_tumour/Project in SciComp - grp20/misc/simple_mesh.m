% Make basic mesh
function [P,E,T] = simple_mesh(deg)
    if deg == 2
    % 2x2 mesh
    P = [0 1 1 0; 
         0 0 1 1]; % Point matrix

    T = [1 1;
         2 3;
         3 4;
         1 1]; % Connectivity matrix

    E = [1 2 3 4;   
         2 3 4 1;
         0 0 0 0;
         1 1 1 1;
         1 2 3 4;
         1 1 1 1; 
         0 0 0 0]; % Edge matrix
    else
    % 3x3 mesh
    P = [0 1 2 0 1 2 0 1 2; 
         0 0 0 1 1 1 2 2 2]; % Point matrix

    T = [1 1 2 2 4 4 5 5;
         2 5 3 6 5 8 6 9;
         5 4 6 5 8 7 9 8;
         1 1 1 1 1 1 1 1]; % Connectivity matrix
    
    E = [1 2 3 6 9 8 7 4;   
         2 3 6 9 8 7 4 1;
         0 0 0 0 0 0 0 0;
         1 1 1 1 1 1 1 1;
         1 2 3 4 5 6 7 8;
         1 1 1 1 1 1 1 1; 
         0 0 0 0 0 0 0 0]; % Edge matrix
     
    end
end