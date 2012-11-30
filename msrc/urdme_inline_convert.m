%function urdme_inline_convert(umod,outputfile)
%   OR
%function urdme_inline_convert(umod,outputfile,modelfile)
%
%INPUT: 
% umod.M1 :: Number of inline propensities.  This must be equal to the
%                 number of reactions
% umod.K  :: (3 x M1) array. Each column corresponds to a reaction. Contains
%                 the rate constants for the reactions. For each colume, only one
%                 of the values can be non-zero.  If it is row 1 then the reaction
%                 is bi-molecular, if it is row 2 then the reaction is mono-molecular,
%                 if it is row 3 the reaction is zero-order.  Note that if the 
%                 "umod.parameter" vector is defined, then the rate constant 
%                 will be replaced by the parameter placeholder.
%       
% umod.I  :: (3 x M1) array.  Each column corresponds to a reaction.  If the 
%                 reaction is bi-molecular, then the value of the first two rows are
%                 the indicies of the reagents of the reaction.  If the reaction is 
%                 mono-molecular, then the value of row 3 is the index of the reagent.
%               
% umod.S  :: (Num-Subdomains x M1) array. Each column corresponds to a reation.
%                 Non-zero rows indicate subdomains where the corresponding reaction
%                 is DISABLED.  If row 1 has a non-zero element, then the reaction
%                 propensity will be zero for all voxel in that subdomain.  Note that
%                 subdomain numbering starts with 1.  Inline propensities are not 
%                 active in a subdomain labeled '0'.
%
%OUTPUT:
% outputfile :: File where generated C propensity functions are written.
%
%OPTIONAL:
% modelfile :: The urdme model '.m' file, typically the calling function. If this
%              option is provided, the output file is written only if the 
%              change date on the '.m' file is newer than the change date on the
%              '.c' file, or if the '.c' file does not exist.
%
function urdme_inline_convert(umod,outputfile,varargin)
    %make sure we have M1,K,S,I defined
    if(~(isfield(umod,'M1')&&isfield(umod,'K')&&isfield(umod,'S')&&isfield(umod,'I')))
        error('urdme_inline_convert(): M1,K,S,I must be defined in umod');
    end
    [Mspecies,Mreactions]=size(umod.N);
    [NumSubdomain,M1]=size(umod.S);
    if(isfield(umod,'parameters') && size(umod.parameters,1)~=Mreactions)
        error('When using both inline propensities and parameters, the first dimension of "umod.parameters" must be equal to "Mreactions"');
    end
    if(Mreactions~=M1)
        error('Number of inline propensities (umod.M1) does not match Mreactions = size(umod.N,1)');
    end
    %%%%%%%%%%%%%%%%
    %find out if  the model file is newer than the propensity file
    if(nargin==3)
        list=dir(outputfile);
        if(length(list)>0)
            %fprintf('%s exists, changedate: %s\n',outputfile,list(1).date);
            listin=dir(varargin{1});
            if(length(listin)>0)
                %fprintf('%s exists, changedate: %s\n',varargin{1},listin(1).date);
                if(datenum(listin(1).date) > datenum(list(1).date))
                    fprintf('model file "%s" newer than propensity file "%s", updating\n',varargin{1},outputfile);
                else
                    fprintf('model file "%s" not newer than propensity file "%s", stopping\n',varargin{1},outputfile);
                    return;  %
                end
            else
                %fprintf('%s does not exist\n',varargin{1});
            end
        else
            %fprintf('%s does not exist\n',outputfile);
        end
    end
    %%%%%%%%%%%%%%%%
    % open outfile
    %fh = 1; %stdout
    fh=fopen(outputfile,'w+');
    if(fh==-1)
        error(sprintf('can not write to file %s',outputfile));
    end
    %%%%%%%%%%%%%%%%
    print_header(fh);
    for rxn_num=1:umod.M1
        print_rxn_fn(fh,rxn_num);
    end
    print_footer(fh);
    if(fh>1)
        fclose(fh);
    end
    %%%%%%%%%%%%%%%%
    function print_header(fh)
        fprintf(fh,'#include <stdlib.h>\n#include <stdio.h>\n#include "propensities.h"\n\n');
    end%print_header()
    function print_rxn_fn(fh,rxn_n)
        fprintf(fh,'double rFun%i(const int *x, double t, const double vol, const double *data, int sd){\n',rxn_n);
        sdv=find(umod.S(:,rxn_n));
        for sdn=1:length(sdv) %disabled subdomains
            fprintf(fh,'if(sd==%i){return 0.0;}\n',sdv(sdn));
        end
        if(isfield(umod,'parameters'))
            fprintf(fh,'return parameters[%i]',rxn_n-1);
        else
            fprintf(fh,'return (%e)',max(umod.K(:,rxn_n)));
        end
        if(umod.K(1,rxn_n)~=0)%bi-molecular
            if(umod.I(1,rxn_n)==umod.I(2,rxn_n)) %homo
            fprintf(fh,'/vol*x[%i]*(x[%i]-1)/2',umod.I(1,rxn_n)-1,umod.I(2,rxn_n)-1);
            else %hetro
            fprintf(fh,'/vol*x[%i]*x[%i]',umod.I(1,rxn_n)-1,umod.I(2,rxn_n)-1);
            end
        elseif(umod.K(2,rxn_n)~=0)%mono-molecular
            fprintf(fh,'*x[%i]',umod.I(3,rxn_n)-1);
        else %zero-order
            fprintf(fh,'*vol');
        end
        fprintf(fh,';\n}\n');
    end%print_rxn_fn
    function print_footer(fh)
        fprintf(fh,'PropensityFun*ALLOC_propensities(void){\n');
        fprintf(fh,'PropensityFun*ptr = (PropensityFun*)malloc(sizeof(PropensityFun)*%i);\n',Mreactions);
        for r=0:Mreactions-1
            fprintf(fh,'ptr[%i]=rFun%i;\n',r,r+1);
        end
        fprintf(fh,'return ptr;\n}\n');
        fprintf(fh,'void FREE_propensities(PropensityFun*ptr){free(ptr);}\n');
    end%print_footer
end%function urdme_inline_convert
