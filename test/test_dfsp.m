% Test the dfsp solver as well as the supported syntaxes for the input 
% arguments
%
% B. Drawert, A. Hellander, 13/12/12. 
% 

clc;

umod.name = 'annihilation';

wd = pwd;
cd ../examples/annihilation
outfile= 'test_dfsp_out.mat';
system(['rm ' outfile])

try
% Execute the solver. This will autimatically generate a propensity file
% from the inline proponsties described in the function annihilation.m
umod = urdme(umod,@annihilation,'outfile',outfile,'report',0);
assert(logical(file_exists(outfile)));
system(['rm ' outfile])

umod = urdme(umod,'annihilation','outfile',outfile,'report',0);
assert(logical(file_exists(outfile)));
system(['rm ' outfile])

umod = urdme(umod,'annihilation',{'Solver','dfsp','verbose',2,'report',0,'outfile',outfile});
assert(logical(file_exists(outfile)));
system(['rm ' outfile])

umod = urdme(umod,'annihilation',{'Solver','dfsp','tau',2.5e-2,'max_jump',10,'verbose',2,'report',0,'outfile',outfile});
assert(logical(file_exists(outfile)));
system(['rm ' outfile])

umod = urdme(umod,'annihilation',{'Solver','dfsp','verbose',1,'DFSP_cache','ann_dfsp_cache_file.mat','report',0,'outfile',outfile});
assert(logical(file_exists(outfile)));
system(['rm ' outfile])
cd(wd)
catch err
    cd(wd);
    rethrow(err)
end

disp('Passed')

