function OT=compute_overturns_discrete_AP(p,t,s,Params)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% function OT=compute_overturns_discrete_AP(p,t,s,Params)
%
% Compute overturns (Thorpe scales) and epsilon for a given density/temperature profile .
%
%---------
% INPUTS:
%---------
% p     : pressure    [db] , M X 1
% t     : temperature [C]  , M X 1
% s     : salinity    [PSU], M X 1
% Params: Structure with (optional) parameters:
%   lat       : latitude [deg], scalar or vector (defaul=30)
%   usetemp   : 1 to use temperature, 0 to use density (default)
%   minotsize : Min OT size (too small and may be noise) (default=2)
%   sigma     : Noise level for density. (default=5e-4)
%   runlmin   : Runlength minimum (default=0)
%   plotit    : Make a summary plot (default=0)
%
%---------
% OUTPUTS:
%---------
% OT : Structure with fields:
%
% eps         :  Epsilon [Wkg^-1], M X 1
% Lttot       :  Thorpe scale [m] , M X 1
% Lmin        :  Min OT [m] resolvable from noise in density and N2 , M X 1
% runlmax     :  max run length
% Lot         :  Vertical size of overturn [m] M X 1
% ptmp        :  Pot. temperature (not actually returned yet?)
% Otnsq_out   :  Mean N2 in overturn region (used to calc. eps),M X 1
% d           :  Displacements when profile is re-ordered (m), M X 1
%
%   refd   : Reference density(s) used to compute potential density
%   Num_OT : # of overturns that passed all criteria/tests
%   'each' fields return only one value for each overturn (for example if an
%   overturn spanned 6 depth values, Lot would return the (same) value of Lot at all
%   6 depth points, while Lot_each would return only the single value. This
%   could be useful for looking at statistics of Lt, eps, etc. since a
%   histogram of the full outputs would be weighted toward larger
%   overturns. Output is M X 1, with all NaNs except for the first
%   Num_OT values, where Num_OT is the # of OT that passed tests.
%   Lot_each
%   Lt_each
%   Otnsq_each
%   eps_each
%
%---------
% TO DO:
% Different ref densities?
%
%
% Dependencies:
% sw_pden.m
% sw_ptmp.m
% sw_bfrq.m
% clip.m
%
%-----------------
% A. Pickering - andypicke@gmail.com
%
% ~ NOTES - AP ~
% Starting with function compute_overturns_discrete.m, which I reeived from
% J. Nash. He originally got code from Jen MacKinnon. 13 Feb 2015 - AP
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

% Set defaults if parameters not specified
if ~isfield(Params,'lat')
    Params.lat=30;
end

if ~isfield(Params,'minotsize')
    Params.minotsize=2;
end

if ~isfield(Params,'runlmin')
    Params.runlmin=0;
end

if ~isfield(Params,'sigma')
    Params.sigma=5e-4;
end

if ~isfield(Params,'usetemp')
    Params.usetemp=0;
end

if ~isfield(Params,'plotit')
    Params.plotit=0;
end

lat      = Params.lat;
minotsize= Params.minotsize;
runlmin  = Params.runlmin;
sigma    = Params.sigma;
usetemp  = Params.usetemp;

% make sure inputs are columnn vectors
t=t(:);
s=s(:);
p=p(:);

%% make potential density and temp at reference depths
% Use one depth if total depth range <1200 m (e.g. fast ctd), but use
% several depth reference levels for a broader range (e.g. shipboard ctd)

% if (max(p)-min(p))>1500
%     dref=1000;
%     refd=(min(p)+dref/2):dref:max(p);
% else
refd=(min(p)+max(p))/2;
dref=(max(p)-min(p))+1;
%end

%%
% Make empty arrays for results
Epsout=NaN*p(:);
Lmin=NaN*t;
Lot=NaN*t;
Lttot=Lot;
Otnsq_out=Lot;
d=Lot;

for iref=1:length(refd)
    
    % compute potential density and temperature
    pden = sw_pden(s(:),t(:),p(:),refd(iref));
    ptmp = sw_ptmp(s(:),t(:),p(:),refd(iref));
    
    % choose temperature or density to use
    if usetemp
        V=ptmp;
    else
        V=pden;
    end
    
    % sort density profile
    [xx,isort]=sort(pden);
    
    % compute N^2 with sorted (stable) density profile
    [n2,q,p_ave] = sw_bfrq(s(isort),t(isort),p,lat);
    
    p0=p(:);  % full depth
    
    % find good (not NaN density values)
    ig=find(~isnan(pden));
    
    if numel(ig)>1 % only do if we have good data
        
        % make vectors with only good values
        pg=p(ig);
        ptmp=ptmp(ig);
        pden=pden(ig);
        V=V(ig);
        
        % make sure presure is a column vector
        pg=pg(:);
        
        % distinguish between up/down casts?
        sig = sign(nanmedian(diff(V)));
        
        % sort profile
        [Vsort,ind]=sort(sig*V);
        % tsort=sig*V; % AP 13 Feb
        psort = pg(ind);
        dz = pg-psort;
        
        csdz = cumsum(-dz);
        thresh = 0.0000001;
        
        % Find start and stop indices for overturn region(s)
        start = find(csdz(1:end-1)<thresh & csdz(2:end)>=thresh)+1;
        if dz(1)<0
            start = [1;start];
        end;
        stops = find(csdz(1:end-1)>=thresh & csdz(2:end)<thresh)+1;
        
        % Return a list of parameters for each overturn region (1 value per overturn)
        % this will be useful for looking at statistics etc. In full profiles,
        % these values are repeated for each depth value in overturn range, so
        % for example histograms of Lt would be weighted towards larger
        % overturns.
        
        clear Lot_each Lt_each Otnsq_each eps_each
        Lot_each   = [] ;
        Lt_each    = [] ;
        Otnsq_each = [] ;
        eps_each   = [] ;
        
        clear start_pass stops_pass
        start_pass=[];
        
        %~
        if Params.plotit==1
            figure(1);clf
            
            subplot(141)
            plot(n2,p_ave)
            xlim([0 nanmax(n2)])
            xlabel('N^2')
            axis ij
            grid on
            
            subplot(142)
            plot(sig*V,pg,Vsort,pg)
            hold on
            plot(V(start),pg(start),'bo')
            plot(V(stops),pg(stops),'rd')
            axis ij
            grid on
            ytloff
            
            subplot(143)
            plot(dz,pg,'o-')
            xlabel('dz')
            axis ij
            grid on
            ytloff
        end
        %~
        
        % make empty arrays to store results for each overturn region
        Otnsq = NaN*dz;
        Lt=NaN*dz;
        Lmin0=NaN*pg;
        Lot0=NaN*pg;
        runlmax0=Lmin0;
        R0tot=Lmin0;
        
        clear j
        for j = 1:length(start);
            clear ind indp n2avg delz drhodz
            ind=clip([(start(j)-1):(stops(j)+1)],1,prod(size(dz)));
            indp=find(p_ave>min(pg(ind))&p_ave<max(pg(ind)));
            n2avg=nanmean(n2(indp));
            warning off
            delz=abs(max(pg(ind))-min(pg(ind)));
            drhodz=(max(pden(ind))-min(pden(ind)))/delz;
            
            % compute run-length
            stopnow=0; runl=1;
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
            
            % Minimum resolvable overturn given N2 and error in density
            Lmin0(ind)=2*9.8/n2avg*sigma/1027;
            
            % Vertical size (pressure) of overturn region
            Lot0(ind)=(max(pg(ind))-min(pg(ind)));
            
            % Max density difference over overturn region
            drho=(max(pden(ind))-min(pden(ind)));
            
            %    if (delz>minotsize)&  (length(ind)>(10*(sigma/drho)^2))  % jody's suggestion
            % additional test from Gargett and Garner 08
            Lpos=length(find((V(ind)-sort(V(ind)))>0));
            Lneg=length(find((V(ind)-sort(V(ind)))<0));
            R0=min(Lpos/length(ind),Lneg/length(ind));
            
            % Check if overturn passes tests/criteria
            if (delz>minotsize) && (delz>(2*9.8/n2avg*sigma/1027))...
                    && ( (max(pden(ind))-min(pden(ind))) > (2*sigma) )...
                    && runl>runlmin ...
                    && ( max(abs(V(ind)-sort(V(ind)))) > 2*sigma  )...
                    && R0>0.2
                
                Otnsq(ind) = 9.8./mean(pden(ind)).*drhodz;
                Lt(ind)=sqrt(mean(dz(ind).^2));
                R0tot(ind)=R0;
                
                Lot_each   = [Lot_each (max(pg(ind))-min(pg(ind))) ];
                Lt_each    = [Lt_each sqrt(mean(dz(ind).^2))];
                Otnsq_each = [Otnsq_each 9.8./mean(pden(ind)).*drhodz ];
                eps_each   = [eps_each 0.64*Lt_each(end).^2.*sqrt(Otnsq_each(end)).^3];
                
                start_pass=[start_pass start(j)];
                
            else % overturn did not pass test
                Otnsq(ind) = NaN;
                Lmin0(ind) = NaN;
                Lot0(ind)  = NaN;
                Lt(ind)    = NaN;
                R0tot(ind) = NaN;
                
            end % if pass tests
            
            
        end; % each overturn
        
        clear iz
        iz=find( p0>(refd(iref)-dref/2) & p0<=(refd(iref)+dref/2)) ;
        
        clear xxx iun
        [xxx,iun]=unique(pg);
        Lt=Lt(:);
                
        % Thorpe scale (rms displacement)
        Lttot(iz)=interp1(pg(iun),Lt(iun),p0(iz));
        
        % Epsilon
        Epsout(iz) = interp1(pg(iun),0.64*Lt(iun).^2.*sqrt(Otnsq(iun)).^3,p0(iz));;
        
        Otnsq_out(iz)=interp1(pg(iun),Otnsq(iun),p0(iz));
        
        Lmin(iz)=interp1(pg(iun),Lmin0(iun),p0(iz));
        
        % Vertical size of reording region
        Lot(iz)=interp1(pg(iun),Lot0(iun),p0(iz));
        
        % Maximum run-length?
        runlmax(iz)=interp1(pg(iun),runlmax0(iun),p0(iz));
        
        % Assign any NaNs in epsilon value of 1e-11
        Epsout(isnan(Epsout))=1e-11;
        
        % return ptmp and pdens (not sorted) also (interp back to original z vector)
        p_tmp=interp1(ptmp,pg,p0);
        
        p_den=interp1(pden,pg,p0);
        
        d=interp1(pg,dz,p0);
        
    else
        % not enough points to interpolate, keep nan
    end
    
end % each ref. density

n2out=interp1(p_ave,n2,p0); % NOTE this maybe not n2 used to compute epsilon in each overturn?

if Params.plotit==1
    subplot(144)
    plot(Lot,p0,'o-')
    xlabel('Lot')
    hold on
    plot(Lmin,p0)
    axis ij
    grid on
    ylim([nanmin(p0) nanmax(p0)])
    ytloff
    linkaxes(get(gcf,'Children'),'y')
end

% return some variables in an output structure
OT=struct();
OT.refd=refd;
OT.Params=Params;
OT.MakeInfo=['Made ' datestr(now) ' w/ ' mfilename ', in ' version];

OT.eps=Epsout;
OT.p=p0;
OT.Lmin=Lmin;
OT.Lot=Lot;
OT.runlmax=runlmax;
OT.Lttot=Lttot;
OT.p_tmp=p_tmp;
OT.n2out=n2out;
OT.Otnsq_out=Otnsq_out;

OT.d=d;
OT.Num_OT=numel(start_pass);
OT.Lot_each=Lot_each;
OT.Lt_each=Lt_each;
OT.Otnsq_each=Otnsq_each;
OT.eps_each=eps_each;
%
return

%%