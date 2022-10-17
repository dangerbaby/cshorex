% Some housecleaning
clear all
close all

addpath(genpath('mfiles'));

% Name assign for site
name = 'planar';
name = 'barsed';
%name = 'frf_runup';
%name = 'lstf';
addpath(name)

% set the input params or read infile
% read_infile(g)
params

%run simplified cshore 
out=cshore(in);

%plot
%plot_results

 