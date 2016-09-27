%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Eval_N2_Codes.m
%
% 
%
%-----------------
% 09/27/16 - A.Pickering - andypicke@gmail.com
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

path_to_repo = '/Users/Andy/Standard-Mixing-Routines/'

%addpath(fullfile(path_to_repo,'ThorpeScales'))
addpath(fullfile(path_to_repo,'seawater_ver3_2'))

% load CTD data file
% 'ctd' structure contains downcast CTD profiles from a ~24-hr station on R/V
% Revelle during 2010 IWISE experiment in Luzon Strait.
load(fullfile(path_to_repo,'Data','CTD_N2a.mat'))

% choose cast to use
icast=5

% Method (1) use sw_bfreq from sw ver3.2 library
[n2,q,p_ave]=sw_bfrq(ctd.s(:,icast), ctd.t(:,icast), ctd.p(:), ctd.lat(icast));

figure(1);clf
plot(n2,p_ave)
axis ij
grid on
ylabel('p [db]','fontsize',15)
xlabel('N^2 [rad^{-1}','fontsize',15)

% Method (2) the GSW (TEOS-10) codes



%%