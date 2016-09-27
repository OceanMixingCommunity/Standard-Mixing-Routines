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
addpath(fullfile(path_to_repo,'seawater_ver3_2'))

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
icast=8
p=ctd.p;
t=ctd.t(:,icast);
s=ctd.s(:,icast);
OT=compute_overturns_discrete_AP(p,t,s,Params)

%% Plot results for 1 profile

figure(1);clf

ax1=subplot(131);
semilogx(t,p,'linewidth',2)
axis ij
grid on
ylabel('p [db]')
xlabel('t [^oC]')
title('Temperature')

ax2=subplot(132);
plot(OT.Lot,OT.p,'d-')
axis ij
grid on
xlabel('Lot')
title('Thorpe displacement')
set(gca,'YTickLabel',[])

ax3=subplot(133);
semilogx(OT.eps,OT.p,'d')
axis ij
grid on
xlabel('log_{10}[\epsilon]')
title('Dissipation rate')
set(gca,'YTickLabel',[])
set(gca,'XTick',[1e-10 1e-9 1e-8 1e-7 1e-6 1e-5])
xlim([1e-11 nanmax(OT.eps)])

linkaxes([ax1 ax2 ax3],'y')

%%