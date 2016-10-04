function [n2,pout,dthetadz,dsdz] = g_nsqfcn(s,t,p,p0,dp)

% G_NSQFCN Calculate buoyancy frequency from a ctd profile
%
%     [N2,POUT,DTHETADZ,DSDZ] = g_nsqfcn(S,T,P,P0,DP) Calculates buoyancy
%     for profiles of salinity, temperature and pressure.
%     The Function: (1) low-pass filters t,s,p over dp,
%                   (2) linearly interpolates  t and s onto pressures, p0,
%                       p0+dp, p0+2dp, ....,
%                   (3) computes upper and lower potential densities at 
%                       p0+dp/2, p0+3dp/2,...,
%                   (4) converts differences in potential density into nsq
%                   (5) returns NaNs if the filtered pressure is not
%                       monotonic.
%     
%     If you want to have the buoyancy frequency in [cyc/s] then
%     calculate sqrt(n2)./(2*pi).
%     For the period in [s] do sqrt(n2).*2.*pi
%
%     INPUT   s  - salinity (concentration units) in a column array,
%             t  - in-situ temperature (deg C) in a column array
%             p  - pressure (dbar) in a column array
%             p0 - the lower bound for pressures (dbar) used for nsq
%             dp - the pressure interval (dbar) over which potential
%                  density is first-differenced on the reference surface
% 
%     OUTPUT  n2       - the square of the buoyancy frequency, in (rad/s)^2
%             pout     - pressure (dbar) at the center of the nsq windows
%             dthetadz - the gradient of potential temperature
%             dsdz     - the salinity gradient
% 
%     Adapted from Gregg? and Alford
%
%     Gunnar Voet
%     gvoet@ucsd.edu
%
%     Last modification: 06/28/2013

G  = 9.80655;
dz = dp;

% Make sure data comes in rows
if isrow(s); s = s'; end
if isrow(t); t = t'; end
if isrow(p); p = p'; end

% Delete negative pressures
i = find(p>=0);
p = p(i);
s = s(i);
t = t(i);

% Exclude nans in t and s
i = find(~isnan(s) & ~isnan(t));
p = p(i);
s = s(i);
t = t(i);

% Put out all nan if no good data left
if isempty(p)
  n2=NaN; pout=NaN; dthetadz=NaN; dsdz=NaN;
else

% Reverse order of upward profiles
if p(length(p))<p(1)
	p = flipud(p);
  t = flipud(t);
  s = flipud(s);
end

% Low pass temp and salinity to match specified dp
dp_data = diff(p);
dp_med  = median(dp_data);
%[b,a]=butter(4,2*dp_med/dp); %causing problems...
a = 1;
b = hanning(2*floor(dp/dp_med));
b = b./sum(b);

tlp = filtfilt(b,a,t); 
slp = filtfilt(b,a,s);
plp = filtfilt(b,a,p);

% Check that p is monotonic
diffp = diff(plp);
i = find(diffp>0);
if length(i)==length(plp)-1;

	pmin = plp(1);
	pmax = plp(length(plp));
  
%   % Sort density if opted for
%   if sort_dens
%     rho = sw_pden(slp,tlp,plp,plp);
%     [rhos, si] = sort(rho,'ascend');
%     tlp = tlp(si);
%     slp = slp(si);
%   end

	% Determine the number of output points
	i    = find(plp<=pmax-dp);
	npts = length(i);

	while p0<=pmin
		p0 = p0+dp;
	end

	% End points of nsq window
	pwin = (p0:dp:pmax)'; 
	t_ep = interp1(plp,tlp,pwin);
	s_ep = interp1(plp,slp,pwin);
	npts = length(t_ep);

	% Compute pressures at end points
	pout=(p0+dp/2:dp:max(pwin))';
  

	% Compute potential density of upper window pts at output pressures
	i_u  = (1:1:length(pwin)-1);
	pd_u = sw_pden(s_ep(i_u),t_ep(i_u),pwin(i_u),pout);

	% Compute potential density of lower window pts at output pressures
	i_l = (2:1:length(pwin));
	pd_l = sw_pden(s_ep(i_l),t_ep(i_l),pwin(i_l),pout);
	
  % Compute buoyancy frequency squared
  n2 = G*(pd_l - pd_u)./(dp.*pd_u);
  
  % Compute gradients
	dthetadz = (sw_ptmp(s_ep(i_u),t_ep(i_u),pwin(i_u),0) - ...
              sw_ptmp(s_ep(i_l),t_ep(i_l),pwin(i_l),0))/dz;
  dsdz     = (s_ep(i_u)-s_ep(i_l))/dz;
else
	disp('  filtered pressure not monotonic')
	n2=NaN; pout=NaN; dthetadz=NaN; dsdz=NaN;
end
end

