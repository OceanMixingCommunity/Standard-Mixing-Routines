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
% 'ctd' structure contains downcast CTD profiles from a ~24-hr station on R/V
% Revelle during 2010 IWISE experiment in Luzon Strait.
load(fullfile(path_to_repo,'Data','CTD_N2a.mat'))

figure(1);clf

ax1=subplot(211);
ezpc(ctd.datenum,ctd.p,ctd.t)
cb=colorbar;
cb.Label.String='^oC';
datetick('x');
ylabel('p [db]')

ax2=subplot(212);
ezpc(ctd.datenum,ctd.p,ctd.s)
cb=colorbar;
cb.Label.String='psu';
datetick('x');
ylabel('p [db]')

linkaxes([ax1 ax2])
%%

% (1) compute_overturns_discrete_AP.m

% run for 1 cast
Params.lat=nanmean(ctd.lat);
Params.plotit=1
icast=1
p=ctd.p(:,icast);
t=ctd.t(:,icast);
s=ctd.s(:,icast);
[Epsout,Lmin,Lot,runlmax,Lttot,p_tmp,n2out,Otnsq_out,OT]=compute_overturns_discrete_AP(p,t,s,Params)

%%

figure(1);clf

subplot(121)
semilogx(t,p,'linewidth',2)
axis ij
grid on
ylabel('p [db]')
xlabel('t [^oC]')

subplot(122)
semilogx(OT.eps,OT.p,'d-')
xlim([1e-11 nanmax(OT.eps)])
axis ij
grid on
ylabel('p [db]')
xlabel('log_{10}[\epsilon]')
%%