% Minimal working example of fitting spectra with MLE (assuming gaussian observations, then chi distributed in spectral form). 
% You can replace with other models (Nasmyth) and spectral observations
% from a shear probe
load('VerticalSpectra','meanU','f','P','dof')% units of Hz

veldir=3; % suppleid vertical spectra
ind=find(f>3e-1); %lets fit these frequencies 
SpecObs.k=f(ind)./meanU; % in rad/m
SpecObs.P=P(ind)*meanU; % in rad/m to match up with modelSpec. Use vertical direction for now

%%
idmModel=@(epsi)(inertial_model(epsi,SpecObs.k,veldir)); % creating a function handle of the model used during fitting (only epsilon is allowed to vary)

[epsilon, std_error]=mle_any_model(SpecObs,dof,[1e-10 1e-2], idmModel, 1);



