%%URDME_COMPILE. 
%
% Used internally by urdme.m to compile the solvers when 
% running urdme in interactive mode.  
%
% B. Drawert, 2009-2012
%
% urdme_compile(fem,model_name,propensity_file,solver_name,varagin)

function fem = urdme_compile(fem,model_name,propensity_file,solver,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
verbose=0;
if(nargin==5 && varargin{1}~=0)
    verbose=varargin{1};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fprintf('urdme_compile(): model_name=%s propensities=%s solvers=%s verbose=%i \n',model_name,propensity_file,solver,verbose);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% URDME initialization
if isempty(getenv('URDME_ROOT'))
    fprintf('URDME_ROOT environmental variable not set, attempting to set it\n');
    [s,r] = system('urdme_init -r');
    uroot = strtok(r);
    if(isdir(uroot)),setenv('URDME_ROOT',uroot),end
    if(isempty(getenv('URDME_ROOT'))) %Only works from the examples dir
        [s,r] = system('../../bin/urdme_init -r');
        uroot = strtok(r);
        if(isdir(uroot)),setenv('URDME_ROOT',uroot),end
    end
    if(isempty(getenv('URDME_ROOT')))
        error('URDME_ROOT environmental variable not set (are both urdme_init and matlab in your PATH?)');
    end
end
basedir = getenv('URDME_ROOT');
path(path,strcat(basedir,'/msrc'));
path(path,strcat(basedir,'/comsol'));
curdir = pwd;
temp = getenv('PATH');
newpath = [curdir './urdme' ':' curdir  ':' basedir '/bin' ];
if( size(strfind(temp,newpath),1)>0)% have we already done this?
    temp = [temp ':' newpath ];
    setenv('PATH',temp);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% URDME Solver Path
urdme_solver_path = cell(1);
str= getenv('URDME_SOLVER_PATH');
if(isempty(str))
    urdme_solver_path{1} = basedir;
else
    if(verbose>1),fprintf('URDME_SOLVER_PATH=%s\n',str);end
    i=1;
    while true
        [tok,str]=strtok(str,':');
        if(isdir(tok))
            if(verbose>1),fprintf('\t%s found\n',tok);end
            urdme_solver_path{i}=tok;
            i=i+1;
        else
            if(verbose>2),fprintf('\t%s not found\n',tok);end
        end
        if(isempty(str)),break;end
    end
    if(i==1)
        fprintf('WARNING: URDME_SOLVER_PATH environmental variable not set properly\n');
        urdme_solver_path{1} = basedir;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find solver root
%for i=1:size(solvers,2)
    %sol = lower(solvers{i});
sol = lower(solver);
solroot = '';
for j=1:length(urdme_solver_path)
   makefile = strcat(urdme_solver_path{j},'/build/Makefile.',char(sol));
   if(file_exists(makefile))
        solroot=urdme_solver_path{j};
        if(verbose>1),fprintf('using: %s\n',makefile);end
        break
    end
    if(verbose>1),fprintf('not found: %s\n',makefile);end
end
if(isempty(solroot))
    error(sprintf('%s solver not found',sol));
end
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check if ".urdme/" exists, if so, was it made with the same solver, if not delete it
if(exist('.urdme/','dir')>0)
    d=dir('.urdme/Makefile.*');
    if(length(d)>0 && strcmp( d(1).name,['Makefile.',sol])~=1)
        if(verbose),fprintf('found .urdme/%s, wrong solver. Removing .urdme/ directory\n',d(1).name);end
        system('rm -r .urdme/');
        system('mkdir .urdme/');
    end
else
    system('mkdir .urdme/');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if propensity file exists, copy it to ".urdme/model_name.c" (if newer)
if(file_exists(propensity_file))
    dest_file = sprintf('.urdme/%s.c',model_name);
%    if(file_exists(dest_file))
%        dest_fi=dir(dest_file);
%        src_fi =dir(propensity_file);
%        % only copy if the src/dest files have differnting change-date or size
%        if(dest_fi.datenum ~= src_fi.datenum ||  dest_fi.bytes ~= src_fi.bytes)
%            copy_dest_file=1;
%        else
%            copy_dest_file=0;
%        end
%    else
%        copy_dest_file=1;
%    end
    copy_dest_file=1;
    if(copy_dest_file)
        if(verbose>1),fprintf('cp %s %s\n',propensity_file,dest_file);end
        [copy_status,copy_msg,copy_msgId]=copyfile(propensity_file,dest_file);
        if(copy_status==0),error('could not write model propensity file to .urdme directory: %s',copy_msg);end
        if(file_exists(strrep(dest_file,'.c','.o')))
            % Delete the propensity .o file to awoid errors if the time stamp
            % of the newly copied source file is too close to the timestamp
            % of the old object file.
            system(['rm ' strrep(dest_file,'.c','.o')]);
        end
        %fprintf(['if(exist([''urdme_propensity_'',',solver,'])==2);\n']);
        if(exist(['urdme_propensity_',solver])==2)
            if(verbose>1),fprintf(['urdme_propensity_',solver,'(fem,dest_file);']);end
            eval(['urdme_propensity_',solver,'(fem,dest_file);']);
        end
    %else
    %    fprintf('No copy\n');
    end
% else check if inline propensities exist, create ".urdme/model_name.c"
elseif(isfield(fem,'M1')) % inline propensities specified
    dest_file = sprintf('.urdme/%s.c',model_name);
    if(file_exists([solroot,'/msrc/urdme_inline_convert_',solver,'.m']))
        orig_path = path();
        path(orig_path,[solroot,'/msrc/']);
        if(verbose>1)fprintf('executing %s(fem,%s) in %s\n',strcat('urdme_inline_convert_',solver),strcat('.urdme/',model_name),[solroot,'/msrc/']);end
        eval(strcat('urdme_inline_convert_',solver,'(fem,''',dest_file,''')'));
        path(orig_path);
    else
        if(verbose>1)fprintf('executing %s(fem,%s.c) in %s\n','urdme_inline_convert',strcat('.urdme/',model_name),[solroot,'/msrc/']);end
        urdme_inline_convert(fem,dest_file);
    end
% else error
else
    error('Propensity file not found. Specify the file containing the reaction propensities using the ''Propensities'' option.');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% copy the Makefile (if newer)
if(verbose>1),fprintf('cp %s %s\n',makefile,strcat('.urdme/Makefile.',sol));end
copyfile(makefile,strcat('.urdme/Makefile.',sol));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compile Solvers
    %cmd = sprintf('make -f %s/build/Makefile MODEL=%s SOLVER=%s SOLVER_ROOT=%s',basedir,model_name,sol,basedir);
    %cmd = sprintf('make -f %s/build/Makefile.%s URDME_MODEL=%s SOLVER=%s SOLVER_ROOT=%s',solroot,sol,model_name,sol,solroot);
cmd = sprintf('make -f .urdme/Makefile.%s URDME_MODEL=%s SOLVER=%s SOLVER_ROOT=%s',sol,model_name,sol,solroot);
if(verbose>1),fprintf('%s\n',cmd);end
[r,s] = system(cmd);
if(r~=0),error(sprintf('%s\nCMD:%s\n',s,cmd)),end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
