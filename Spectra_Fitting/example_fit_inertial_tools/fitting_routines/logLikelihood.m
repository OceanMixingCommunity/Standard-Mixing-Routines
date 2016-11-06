function [logL ind val]=logLikelihood(epst,Pk,Pt,v)
% Doing the loglikehood calc using the Eq in Ruddick et al 2000 and L&O
% methods paper for fitting a theroretical spectra Pt (matrix) to a observed Pk (vector) with
% degrees of freedom v.
% Needs the chi2pdf function in the stats toolbox..
nPxx=repmat(Pk,[1 length(epst)]); % create matrix to match Pt
z=v*nPxx./Pt;%f(z)
Z=chi2pdf(z,v); %f(z)
Y=log(Z./Pt);
logL=sum(Y,1)+length(Pk)*log(v);
[val ind]=max(exp(logL-max(logL)));