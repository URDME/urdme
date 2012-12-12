#!/usr/bin/perl

use Symbol;


my $URDME_ROOT = $ARGV[0];
if(!defined $URDME_ROOT){
    print "Error: must be call with URDME_ROOT as the first argument";
    exit(0);
}


##filename = which('urdme');
##while(~isempty(filename))
##[pathstr, name, ext] = fileparts(filename);
##rmpath(pathstr);
##filename = which('urdme');
##end

my $cleanup_cmd = "filename = which('urdme');while(~isempty(filename)),[pathstr, name, ext] = fileparts(filename);rmpath(pathstr);filename = which('urdme');end";


my $in = gensym;
my $MATLAB_CMD = 'matlab -nojvm -nosplash';
open($in, '|'.$MATLAB_CMD);
syswrite($in,"$cleanup_cmd\n\n");
syswrite($in,"addpath $URDME_ROOT/msrc\n");
syswrite($in,"savepath\n");
sleep(2);
print "exiting";
syswrite($in,"quit\n");
sleep(2);
close($in);

