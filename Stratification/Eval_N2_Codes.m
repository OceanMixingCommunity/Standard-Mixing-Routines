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

s=ctd.s(:,icast);
t=ctd.t(:,icast);
p=ctd.p;

% Method (1) use sw_bfreq from sw ver3.2 library
[n2_sw,q,p_sw]=sw_bfrq(ctd.s(:,icast), ctd.t(:,icast), ctd.p(:), ctd.lat(icast));

figure(1);clf
h1=plot(n2_sw,p_sw)
axis ij
grid on
ylabel('p [db]','fontsize',15)
xlabel('N^2 [rad^{-1}','fontsize',15)

%% Method (2) the GSW (TEOS-10) codes (not sure if i'm doing this
% correctly..)
addpath(fullfile(path_to_repo,'gsw_matlab_v2_0'))
addpath(fullfile(path_to_repo,'gsw_matlab_v2_0','library'))
lat=ctd.lat(icast);
lon=ctd.lon(icast);
% 1st need Absolute Salinity
%[SA, in_ocean] = gsw_SA_from_SP(SP,p,long,lat);
[SA, in_ocean] = gsw_SA_from_SP(s,p,lon,lat);
% then need CT
CT = gsw_CT_from_t(SA,t,p);
% then we can calculuate N2
[n2_gsw, p_gsw, in_funnel] = gsw_Nsquared_CT25(SA,CT,p,lat);

%%

figure(1);%clf
hold on
h2=plot(n2_gsw,p_gsw)
axis ij
grid on
ylabel('p [db]','fontsize',15)
xlabel('N^2 [rad^{-1}','fontsize',15)

%% Method (3) - adiabatic leveling from Amy W.


[n2_adlev] = adiabatic_N2(p,t,s,lat);

figure(1);%clf
hold on
h3=plot(n2_adlev,p,'b')
axis ij
grid on
ylabel('p [db]','fontsize',15)
xlabel('N^2 [rad^{-1}','fontsize',15)

legend([h1 h2 h3],'sw32','gsw','adlev')

%% Method (4) - N^2 code from Gunnar Voet

p0=0;
dp=nanmean(diff(p))
[n2_gv,p_gv,dthetadz,dsdz] = g_nsqfcn(s,t,p,p0,dp);

figure(1);%clf
hold on
h4=plot(n2_gv,p_gv)
axis ij
grid on
ylabel('p [db]','fontsize',15)
xlabel('N^2 [rad^2s^{-s}]','fontsize',15)

legend([h1 h2 h3 h4],'sw32','gsw','adlev','gv')


%% save figure

print('N2_comparisons','-dpng')
