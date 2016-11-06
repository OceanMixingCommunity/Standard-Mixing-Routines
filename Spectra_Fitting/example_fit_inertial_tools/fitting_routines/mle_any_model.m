function [epsilon, std_error]=mle_any_model(specObs,v,epslim, modelSpec, tplot)
% Uses maximum likelihood for any numbers degrees of freedom  for fitting a model spectrum to spectral observations. Assumes turb obs are gaussian, such that they are chi distributed when the spectra is computed. MLE could be modified for use with other statistical distribution.
%  see Ruddick etal, 2000 which inspired these routines
% Required inputs:
%	SpecObs: structure with spectral observations in same "space" that the model Spectrum function (modelSpec) will return. 
%	The naming of the fields in structure must be:
	%	 k: for wavenumber (or cpm) or even frequency in Hz (Rad/s). this must be consistent with the modelSpec units 
%		P= observed spectra power in same units as modelSpec (typically in (units timeserie)^2/(rad/m)
%	v=degrees of freedom of the spectral estimate  (see Emery & Tompson, Data Analysis book) depends on how you coputed it
% Uses:
%   chi2pdf.m, hence you need matlab's statistical toolbox, but this can be simply "coded" into a "searchable" table in matlab. Bit tedious to type in the values from a stats book.
% 
%__Created by CBluteau, used in getting epsilon from ADVs using inertial subrange model (Bluteau et al, 2011 L&O: methods) and also when fitting Nasmyth to shear spectrum (useful for high epsilon). See Bluteau et al 2016 JTECH for more details

    
if length(v)>1
    v=mode(v);
    warning('You supplied a vector of dof, should be 1 value (same) across the spectra. Using the mode(dof)');
end
    %% search limits

epst=linspace(epslim(1),epslim(2),1000); % very large at first...initial search space
 
%% MLE Computations for first initial search

Pt= modelSpec(epst); %alp*f.^(-5/3)*(epst.*U).^(2/3);% theoretical spectrum in rad/s. Matrix,
[logL ind val]=logLikelihood(epst,specObs.P,Pt,v);



%% Narrow into nearest 2 decade range
epst=logspace(floor(log10(epst(ind)))-1, ceil(log10(epst(ind)))+1,100);
Pt= modelSpec(epst); %alp*f.^(-5/3)*(epst.*U).^(2/3);% theoretical spectrum in rad/s. Matrix,
[logL ind val]=logLikelihood(epst,specObs.P,Pt,v);

 %% NArrow search into nearest 2/10th of a decade (better std_error calc)
dec=floor(log10(epst(ind))); 
emin=(floor(epst(ind)./10.^dec)-0.1).*10^dec;
emax=(ceil(epst(ind)./10.^dec)+.1).*10^dec;
epst=logspace(log10(emin), log10(emax),100);

Pt= modelSpec(epst); %alp*f.^(-5/3)*(epst.*U).^(2/3);% theoretical spectrum in 
[logL, ind, val]=logLikelihood(epst,specObs.P,Pt,v);
      

%% Final estimates and errors

epsilon=epst(ind);
if ind+2>length(logL) || ind-2<1
    std_error=10^10; %when this happens, either the spectra is junk or you've specified too wide search range
else
    h=diff(epst); % step size
    h=h(ind);
    df_logL=(-logL(ind+2)+16*logL(ind+1)-30*logL(ind)+16*logL(ind-1)-logL(ind-2))./(12*h^2); % 2nd order differentiation at peak
    std_error=sqrt((-1/df_logL)); % standard error based on Eq.22 of Ruddnick et al. Must sqrt(var(Kb))
end



%% Optional PLots
if tplot % plots similar to those in Ruddick et al
  
    Pt= modelSpec(epsilon); %alp*f.^(-5/3)*(epst.*U).^(2/3);% theoretical spectrum in 
    figure;
    subplot(2,2,1)
    plot(epst,exp(logL-max(logL)),'bo-','markerface','b','markersize',4); hold on; % plot the maximum likelyhood
    yli=get(gca,'ylim');
    plot([epsilon epsilon],yli,'k--'); % estimated epsilon
    plot(epsilon+[-2*std_error -std_error std_error 2*std_error],[0.5 0.5 0.5 0.5],'r-o');
%	plot(epst,misfit,'--','color',[0 0.75 0.2]);
    xlabel('\epsilon');
    ylabel('Relative likelihood')
    title(['1 Std error: ', num2str(std_error,'%2.1e')]);
    set(gca,'ylim',yli);

    subplot(2,2,2)
    loglog(specObs.k, specObs.P); hold on;
   loglog(specObs.k,Pt,'g');
      yli=get(gca,'ylim');
   % plot(fkolm*[1 1],yli,'k--');
    xlabel('units supplied spectra');
    ylabel('PSD')
    title(['\epsilon =',num2str(epsilon,'%2.1e'),' +/- ',num2str(100*std_error*2/epsilon,'%3.1f'),'%'])
    hl=legend('Obs','MLE fit'); legend('boxoff');set(hl,'fontsize',8)

    subplot(2,2,3)
    [nb bb]=hist((v.*specObs.P./Pt),10); 
    bar(bb,nb/length(Pt));hold on;
    hp=plot(bb,chi2pdf(bb,v),'g'); % only valid for Chi 2 dof v
    hl=legend(hp,'Theoretical pdf'); legend('boxoff');set(hl,'fontsize',8)
    xlabel('{dof} P_{obs}/P_{model}') 
    ylabel('Relative frequency')

    subplot(2,2,4)
    probplot('exponential',v*specObs.P./Pt); box on;
    xlabel('dof P_{obs}./P_{model}')
  
end
