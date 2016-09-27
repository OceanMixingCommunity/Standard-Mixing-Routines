function [Epsout,Lmin,Lot,runlmax,Lttot]=compute_overturns(p,t,s,varargin)


%  [Epsout,Lmin,Lot,runlmax,Lttot]= ...
%        compute_overturns2(p,t,s,'lat',lat,'usetemp',usetemp,'minotsize',minotsize,'sigma',sigma,'runlmin',runlmin)

%compute_overturns.m
%
% calculate overturns from a variety of instruments
% p,t,s are vectors of pressure, temperature and salinity typically from downcast only
% usetemp set to 1 to use temperature to compute overturns, otherwise uses density
% minotsize is the minimum overturn size to consider (too small and may be noise)
% sigma is noise level for density.
% runlmin is run length minimum.

% INPUTS:
%   p - pressure vector
%   t - temperature vector
%   s - salinity vector
%   lat - latitude of measurement (needed for BVFREQ calculation)
%   usetemp - variable (either set to 1 if temperature is to be used to compute overturns (default 0))
%   minotsize - minimum overturn size to consider (default 2)
%   sigma - noise level for density (default 5e-4)
%   runlmin - minimum run length (default 0)
%
%
% OUTPUTS:
%   Epsout - epsilon (dissipation rate) calculated by: e = 0.64*Lt^2*N^3
%   Lmin - smallest observable thorpe scale = 2*9.8/N2 *sigma/1027 where N2
%               is the profile averaged N2
%   Lot - maximum possible overturn observed (depth of sampled column of water)
%   runlmax - maximum number of runs
%   Lttot - rms thorpe overturns
%
% same as compute_overturns.m but input is salinity instead of conductivity
% becuase sometimes salinity has been separately despiked already
%
% 
% Default settings
lat = 30;
usetemp =0;
minotsize=2;
sigma=5e-4;
runlmin=0;

% Switch trap parses the varargin inputs
len = length(varargin);
for i = 1:2:len
    switch lower(varargin{i})
        case 'lat'
            lat=varargin{i+1}
        case 'usetemp'
            usetemp=varargin{i+1}
        case 'minotsize'
            minotsize=varargin{i+1}
        case 'sigma'
            sigma=varargin{i+1}
        case 'runlmin'
            runlmin=varargin{i+1}
        otherwise
            % neglect invalid option
    end
end


%% make potential density and temp at reference depths
% Use one depth if total depth range <1200 m (e.g. fast ctd), but use
% several depth reference levels for a broader range (e.g. shipboard ctd)

if (max(p)-min(p))>1500
    dref=1000; refd=(min(p)+dref/2):dref:max(p);
else
    refd=(min(p)+max(p))/2; dref=(max(p)-min(p))+1;
end

%%
Epsout=NaN*p(:);
Lmin=NaN*t; Lot=NaN*t; Lttot=Lot;


for iref=1:length(refd)
    
    %s = sw_salt(c(:)*10/sw_c3515,t(:),p(:));
    pden = sw_pden(s(:),t(:),p(:),refd(iref));
    ptmp = sw_ptmp(s(:),t(:),p(:),refd(iref));

    if usetemp
        V=ptmp;
    else
        V=pden;
    end
    
    % sort the potential density profile
    [xx,isort]=sort(pden);
    
    % smoothed nsq profile
    %if length(t)>600
    %    a=1; b=ones(200,1)/200;
    %    [n2,q,p_ave] = sw_bfrq(nanfilt(b,a,s),nanfilt(b,a,t),p,lat);
    %elseif length(t)>300
    %    a=1; b=ones(100,1)/100;
    %    [n2,q,p_ave] = sw_bfrq(nanfilt(b,a,s),nanfilt(b,a,t),p,lat);
    %else
    
    % calculate the N2/Buoyancy frequency from slightly filtered s and t
    a=1; b=ones(10,1)/10;
    [n2,q,p_ave] = sw_bfrq(nanfilt(b,a,s(isort)),nanfilt(b,a,t(isort)),p,lat);
    %end
    
    % find nan's and find data vecs without nan's
    ig=find(~isnan(pden));
    p0=p(:);  % full depth
    pg=p(ig); 
    ptmp=ptmp(ig); 
    pden=pden(ig); 
    V=V(ig);
    pg=pg(:);
    
    % take into account that T decreases with depth and RHO increases with depth
    sig = sign(nanmedian(diff(V)));
    [tsort,ind]=sort(sig*V);
    tsort=sig*V;
    psort = pg(ind);
    
    % calculate the difference between the unsorted (pg) and sorted depth (psort)
    dz = pg-psort;
    % take the cumulative sum of the depth difference btwn sorted/unsorted
    %    see increasing values where no overturns, and decreases when OTs
    csdz = cumsum(-dz);
    % define a threshhold
    thresh = 0.00000001;
    
    % find start and stop points of the profile (where it exists and has
    %   good data (above the threshhold)
    start = find(csdz(1:end-1)<thresh & csdz(2:end)>=thresh)+1;
    if dz(1)<0
        start = [1;start];
    end;
    stops = find(csdz(1:end-1)>=thresh & csdz(2:end)<thresh)+1;
    Otnsq = NaN*dz; Lt=NaN*dz;
    Lmin0=NaN*pg; Lot0=NaN*pg; 
    runlmax0=NaN*pg;R0tot=NaN*pg;
    for j = 1:length(start);
        % find the indeces of "good" profile as set by threshhold
        ind=clip([(start(j)-1):(stops(j)+1)],1,prod(size(dz)));
        % find indeces of p_ave (pressure of n2) that are between the "good" data
        indp=find(p_ave>min(pg(ind))&p_ave<max(pg(ind)));
        % average N2 of profile
        n2avg=nanmean(n2(indp));
        warning off
        % total depth DeltaZ
        delz=abs(max(pg(ind))-min(pg(ind)));
        % total depth DeltaRHO/delz
        drhodz=(max(pden(ind))-min(pden(ind)))/delz;
        % run length
        stopnow=0; runl=1;
        % find the number of "segments (runl)" separted by overturns 
        ig=find(diff(sign(dz(ind)))==0);
        if ~isempty(ig)
            runl=runl+1;
            ig2=find(diff(ig)==1);
            while stopnow==0
                if isempty(ig2)
                    stopnow=1;
                else
                    ig2=find(diff(ig2)==1);
                    runl=runl+1;
                end
            end
        end
        runlmax0(ind)=runl;
        % calculate smallest possible Thorpe Scale size
        Lmin0(ind)=2*9.8/n2avg*sigma/1027;
        % Find maximum possible overturn (depth of the sampled water column)
        Lot0(ind)=(max(pg(ind))-min(pg(ind)));
        % Find the pot. density difference between the top and bot of water column
        drho=(max(pden(ind))-min(pden(ind)));
        %    if (delz>minotsize)&  (length(ind)>(10*(sigma/drho)^2))  % jody's suggestion
        % additional test from Gargett and Garner 08
        Lpos=length(find((V(ind)-sort(V(ind)))>0)); 
        Lneg=length(find((V(ind)-sort(V(ind)))<0));  % should this be a less than sign???
        R0=min(Lpos/length(ind),Lneg/length(ind));
        
        % Criteria: 
        %   Galbraith & Kelley 1996
        %       - Nyquist says overturns thinner than 2x vertical resolution (del-z) cannot be
        %       sampled
        %       - Lmin from G&K : 5 * (del-z) 
        %       - If density resolution is del-rho, then lower limit form
        %       density is : 2 del-rho / (drhodz) == 2 g/n2avg * del-rho/rho0
        %       - for epsilon, 
        %           based on vertical resolution, minimum detected value is : 25* del-z^2 N^3
        %           based on density resolution,  minimum detected value is : 4* g^2/sqrt(n2avg) * del-rho/ rho0
        
        if (delz>minotsize)&(delz>(2*9.8/n2avg*sigma/1027))&((max(pden(ind))-min(pden(ind)))>(2*sigma))...
                &runl>runlmin&(max(abs(V(ind)-sort(V(ind))))>2*sigma)&R0>0.2
            Otnsq(ind) = 9.8./mean(pden(ind)).*drhodz; 
            temptemp(j)=(max(pg(ind))-min(pg(ind))); 
            Lt(ind)=sqrt(mean(dz(ind).^2));
            R0tot(ind)=R0;
        else
            Otnsq(ind)=NaN; Lmin0(ind)=NaN; Lot0(ind)=NaN; Lt(ind)=NaN; R0tot(ind)=NaN;
        end
    end;
    iz=find(p0>(refd(iref)-dref/2)&p0<=(refd(iref)+dref/2));
    
    [xxx,iun]=unique(pg); Lt=Lt(:);
    Lttot(iz)=interp1(pg(iun),Lt(iun),p0(iz));
    Epsout(iz) = interp1(pg(iun),0.64*Lt(iun).^2.*sqrt(Otnsq(iun)).^3,p0(iz));;
    Lmin(iz)=interp1(pg(iun),Lmin0(iun),p0(iz)); 
    Lot(iz)=interp1(pg(iun),Lot0(iun),p0(iz));
    runlmax(iz)=interp1(pg(iun),runlmax0(iun),p0(iz));
    
    Epsout(isnan(Epsout))=1e-11;
    
    
    
end

