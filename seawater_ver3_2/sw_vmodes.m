function [Vert, Hori, Edep, PVel]=sw_vmodes(z,N,clat,nmodes);
%
% SW_VMODES calculate vertical modes in a flat-bottomed ocean.
%
% [Vert,Hori,Edep,PVel]=sw_vmodes(z,N,clat,nmodes);
% Inputs:
% z is in m (higher positive the deeper you go).
% N is in rad/sec is the buoyancy frequency.  Don't pass in imaginary or
% negative buoyancy frequencies - they will be ignored.  Best done using a
% "sorted" buoyancy profile.  
% clat is a central latitude, 
% nmodes is the number of modes to calculate.
%
% Outputs: 
% Vert: Vertical velocity modes, normalized.
% Hori: Horizontal velocity modes, normalized.
% Edep: Equvalent depth of the mode in m (usefull for deformation radii)
% PVel: phase velocity of the mode in m/s.
%
% Vert and Hori are [mxn] matrices with m the length of z, and n the number
% of modes.  Edep and PVel are [1xn] matrices.
%
% CAVEAT: Solves an eigenvalue equation, which means it will create an mxm
% matrix and solve it.  Don't put in a 4000m drop with data every meter;
% smooth and decimate first!
%
% Basically solves 
%   Gzz + ev*nsq*G = 0
% subject to:
%   G(-d)=0;
%   Gz(0)-g*ev*G(0)=0
% Where G is the mode function, ev the eigenvalue, nsq the buoyancy
% frequency squared and g the gravitational constant as a fcn of latitude.
% 
% Originally written by Benno Blumenthal, 23 July 1981, in FORTRAN
% Modified for SUNS by  CC Eriksen, August 1988
% Translated to MatLab 4.2c J. Klymak, March 1997
%

% The algorithm is based around the fact that the ode can be written in a
% tri-diagonal matrix form:
% [(-1/nsq(z))&dzz]G = ev*G 
% Where the upper, lower and diagonal are calculated so as to give a good
% estimate of dzz of G.  
% i.e. -l(i-2)*G(i-2) + d(i-1)*G(i-1) - u(i)*G(i) = ev*(g(i-1))
% Then we can use Matlabs eig function to find the eigenvalue ev.  (the
% original FORTRAN uses an EISPACK routine called RATQR to do the same
% thing)
%
% Results are in good agreement with the results from Blumenthal's original
% code, but they are not precise.  

% save the input z
z_in=z;

% Check whether first point is at z=0;
nsqin=N.*N;

if z(1)>0.01
  z=[0; z];
  N=[N(1); N];
end;
good=find(N>0&~isnan(N)&isreal(N));
N=N(good); z=z(good);
npts=length(N);

% square n profile, compute dz...
nsq=N.*N;
dz=[z(1); diff(z)];


% calculate nbar...
nbar=N(1)*z(2) + N(npts)*(z(npts)-z(npts-1));
diffz_=z(3:npts)-z(1:npts-2);
nbar = nbar+sum(diffz_.*N(2:npts-1));
nbar=nbar./(2*z(npts));
nbarcy=nbar;
nbar=nbarcy/572.9577951;  % conver to cycles per hour...
    
% compute tridiagonal matrix with free bc...
alat=clat*3.141592654/180;
grav=9.78049*(1.0+5.2884e-3*(sin(alat))^2-5.9e-6*(sin(2.0*alat))^2);
grainv=1.0/grav;

d=zeros(1,npts);
l=0*d;u=0*d;
d(1)=grainv/dz(2);
u(2)=d(1);
l(1)=0;
u(1)=0;

deltaz=diff(z);

d(2:npts-1)=2./(nsq(2:npts-1).*dz(2:npts-1).*dz(3:npts));
l(2:npts-1)=2./(nsq(2:npts-1).*dz(2:npts-1).*(dz(3:npts)+dz(2:npts-1)));
u(3:npts)=2./(nsq(2:npts-1).*dz(3:npts).*(dz(3:npts)+dz(2:npts-1)));

u(npts)=0;
l(npts)=0;

i1=npts-1;

% make the indices to the matrix.
j_=[[1:length(d)] [1:length(d)-1] [2:length(d)]];
i_=[[1:length(d)] [2:length(d)] [1:length(d)-1]];

% make the matrix to solve...
M=sparse(i_,j_,[d -l(2:length(d)) -u(2:length(d))],length(d),length(d));

% compute the eigen values
[V,D]=eig(full(M));

% Get the modes out... Throw the lowest eigenvalue out for some reason...
ev=(sort(diag(D)));
ev=ev(2:length(ev));
phase=-1;
dz=dz';

nptsin=length(z_in);
Vert=zeros(nptsin,nmodes);
Hori=zeros(nptsin,nmodes);
PVel=zeros(1,nmodes);
Edep=zeros(1,nmodes);
for imode=1:nmodes
  phase=-phase;
  dz=zeros(npts,1);
  dz(npts)=0;
  i=npts-1;
  dz(i)=1;
  imax=npts-2;
  for i2=1:imax
    dz(i-1)=-((ev(imode)-d(i))*dz(i)+u(i+1)*dz(i+1))/l(i);
    i=i-1;
  end;
  % normalize dz
  clear sum;
  sum_ = nsq(1)*dz(1)*dz(1)*z_in(2);
  difz_=z(3:npts)-z(1:npts-2);
  sum_ = sum_ + sum(nsq(2:npts-1).*dz(2:npts-1).*dz(2:npts-1).*difz_);
  sum_=sum_*0.5 + grav*dz(1)*dz(1);
  a=z(npts)./(ev(imode)*sum_);
  a=sqrt(abs(a));
  a=phase*a*sign(ev(imode));
  dz=dz.*a;

  % interpolate/extrapolate onto original z grid... 
  dz=interp1(z,dz,z_in);
  edepth=grainv/ev(imode);
  phasev=sqrt(abs(1/ev(imode)));
  const=nbar/phasev;
  phasev=phasev*sign(ev(imode));
  
  emhor=NaN*z_in;
  emver=NaN*z_in;
  % get the vertical and horizontal velocity modes...
  delta1=z_in(2:nptsin-1)-z_in(1:nptsin-2);
  delta2=z_in(3:nptsin)-z_in(2:nptsin-1);
  d1sq=delta1.*delta1;d2sq=delta2.*delta2;
  emver(2:nptsin-1)=const.*dz(2:nptsin-1);
  emhor(2:nptsin-1)=...
      -(d1sq.*dz(3:nptsin)-d2sq.*dz(1:nptsin-2)+...
      (d2sq-d1sq).*dz(2:nptsin-1))./...
      (delta1.*delta2.*(delta1+delta2));
  emver(1)=const*dz(1);
  emhor(1)=-(dz(2)-dz(1))./(z_in(2)-z_in(1));
  emver(nptsin)=const*dz(nptsin);
  emhor(nptsin)=-dz(nptsin-1)./(z_in(nptsin-1)-z_in(nptsin));
  
  % put everybody in their matrices...
  Vert(1:length(emver),imode)=emver;
  Hori(1:length(emver),imode)=emhor;
  PVel(imode)=phasev;
  Edep(imode)=edepth;
  
  Vert(:,imode)=Vert(:,imode)/max(abs(Vert(:,imode)));
  Hori(:,imode)=Hori(:,imode)/max(abs(Hori(:,imode)));

end;



