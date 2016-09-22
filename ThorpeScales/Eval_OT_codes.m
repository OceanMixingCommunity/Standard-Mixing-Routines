%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Eval_OT_codes.m
%
% Script to run overturns codes in repo on CTD data and compare.
%
%---------------
% 09/22/16 - A.Pickering - andypicke@gmail.com
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

path_to_repo = '/Users/Andy/Standard-Mixing-Routines/'

addpath(fullfile(path_to_repo,'ThorpeScales'))

% load CTD data file
% 'ctd' structure contains CTD profiles from a ~24-hr station on R/V
% Revelle during 2010 IWISE experiment in Luzon Strait.
load(fullfile(path_to_repo,'Data','CTD_N2a.mat'))


% (1) compute_overturns_discrete_AP.m

Params.lat=nanmean(ctd.lat);
Params.plotit=1
icast=1
p=ctd.p(:,icast);
t=ctd.t(:,icast);
s=ctd.s(:,icast);
[Epsout,Lmin,Lot,runlmax,Lttot,p_tmp,n2out,Otnsq_out,OT]=compute_overturns_discrete_AP(p,t,s,Params)

%%

figure(1);clf
semilogx(Epsout,ctd.p,'d-')
axis ij
grid on
ylabel('p [db]')
xlabel('log_{10}[\epsilon]')

%%