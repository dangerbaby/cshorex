% Some housecleaning
clear all
close all

addpath(genpath('mfiles'));

% Name assign for site
name = 'planar';
%name = 'barsed';
%name = 'frf_runup';
%name = 'lstf';
%name = 'agate';
%name = 'gee'
%name = 'q3d';
name = 'osu_mangrove';

addpath(name)

% set the input params or read infile
% read_infile(g)
params

%run cshore 
[in,out]=cshorex(in);

%plot
%plot_results

 