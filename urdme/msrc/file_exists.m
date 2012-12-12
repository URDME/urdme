%  TF = file_exits(filename)
%  returns 1 if the file exists on the filesystem, zero otherwise
%
%
% The Matlab built-in exists() will use the m-file search path to determine if a file exists (an m-file).
% If you want to use it to determine if a file exists in a specific place,
% and not "does this exist anywhere in the search path" it will give you false positives.
%

function TF = file_exists(filename)
    if(isempty(filename))
        TF=0;
        return;
    end
    fid = fopen(filename);
    %fprintf('%i = fopen(%s)\n',fid,filename)
    if(fid==-1)
        TF=0;
    else
        fclose(fid);
        TF=1;
    end
end
