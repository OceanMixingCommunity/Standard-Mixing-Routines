function [N2] = adiabatic_N2(press,temp,salin,lat)

plev=20;  % choose reference level depth 

n2_bray = repmat(nan, [size(press,1) size(press,2)]);

for ii=1:size(press,2)
    % [j ii]
    disp([num2str(ii) '/' num2str(size(press,2))])
    pgrid = press(:,ii);
    pgrid = pgrid(~isnan(pgrid)); % discard nans
    for jj=1:length(pgrid)
        pmin_lev = max([pgrid(jj)-plev/2 min(press(:,ii))]);
        pmax_lev = min([pgrid(jj)+plev/2 max(press(:,ii))]);
        icyc = find(press(:,ii)>=pmin_lev & press(:,ii)<=pmax_lev);
        if isempty(icyc)
            n2_bray(jj,ii) = nan;
        else
            clear pbar tbar sbar
            pbar = nanmean(press(icyc,ii)); % nominal reference press
            tbar = nanmean(temp(icyc,ii));
            sbar = nanmean(salin(icyc,ii));
            rhobar = sw_pden(sbar,tbar,pbar,pbar);
            theta = sw_ptmp(salin(icyc,ii),temp(icyc,ii),press(icyc,ii),pbar); % potential temperature referenced to pbar
            sv = 1./sw_pden(salin(icyc,ii),theta,pbar,pbar); % specific volume a water parcel would have if moved adiabatically to pbar

            % regress press against de-meaned sv and store coefficients
            order = 1; % use linear fit for now
            [s dum] = polyfit(press(icyc,ii),sv-nanmean(sv),order);
            alpha(:,jj,ii) = s';
            alpha(2,jj,ii) = alpha(2,jj,ii)+nanmean(sv);
            clear s dum
            g = sw_g(lat(ii),sw_dpth(pbar,lat(ii)));
            n2_bray(jj,ii) = -1e-4 * rhobar^2 * g^2 * alpha(1,jj,ii);
        end
    end

end

N2 = n2_bray;
%N2(hrp.N2<0)=NaN;

