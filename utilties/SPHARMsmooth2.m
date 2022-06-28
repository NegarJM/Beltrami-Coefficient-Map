 

function [surf_smooth, fourier]=SPHARMsmooth2(surf,sphere, L,sigma)
%---------------------------------------------------------------------------------------
%[surf_smooth, fourier]=SPHARMsmooth2(surf,sphere, L,sigma)
%
% WARNING: This version is different from SPHARMsmooth.m in
% http://www.stat.wisc.edu/~mchung/softwares/weighted-SPHARM/SPHARMsmooth.m
% SPHARMsmooth.m is based on a fixed connectivity across all surface. If
% every subjects have different mesh topology, use SPHARM.smooth2.m.
%
% Parameters: You need to read reference [1] to understand prameters.
%
%    L               : The maximal degree of the weighted-SPHARM representation.
%                       Read the paper below to find it optimally.
%
%  sigma          : bandwith of weighted-SPHARM representation
%                       It is the bandwidth of heat kernel on unit sphere.
%                       When sigma=0, it is the traditional SPHARM representation.
%                       range beween 0.0001 and 0.01 will be sufficient for cortical surfaces.
%
%  sphere         : Spherical mesh obtained from LAPLACEcontour.m and REGUARLIZEarea.m
%
% surf_smooth  : The weighted-SPHARM result.
%
% fourier   : The estimated SPHARM coefficients (Fourier coefficients) given as a structured array
%               containg coeff.x, coeff.y, coeff.z
%               coeff.x is the SPHARM coefficients of x-cooridinate given as (L+1) by (2*L+1) matrix
%               coeff.x(3,:) is the 2nd degree SPHARM coefficients of all order.
%
%
% (C) Moo K. Chung, 2006-2008
%     Department of Biostatistics and Medical Informatics
%     Waisman Laboratory for Brain Imaging
%     University of Wisconsin-Maison
%
% email://mkchung@wisc.edu
% http://www.stat.wisc.edu/softwares/weighted-SPHARM/weighted-SPHARM.html
%
% If you use this code, please reference [1].
% You need to read the paper to understand the notations and the algorithm.
%
% Reference:
% [1] Chung, M.K., Dalton, K.M., Shen, L., L., Evans, A.C., Davidson, R.J. 2007.
% Weighted Fourier series representation and its application to quantifying
% the amount of gray matter. IEEE Transactions on Medical Imaging, 26:566-581.
%
%
% Update history Sept 19 2006; July 5, 2007; January 22, 2008.
%----------------------------------------------------------------------------------------

coord=surf.vertices;
n_vertex = size(coord,1);   % the number of vertices in a mesh.

[theta varphi]=EULERangles(sphere);

x=coord(:,1);
y=coord(:,2);
z=coord(:,3);

%INITIALIZATION

% xestiamte is the weighted-SPHARM of x-coordinate.
xestimate=zeros(n_vertex,1);
yestimate=zeros(n_vertex,1);
zestimate=zeros(n_vertex,1);

% betax is the Fourier coefficients of x-coordinate.
betax=zeros(L+1,2*L+1);
betay=zeros(L+1,2*L+1);
betaz=zeros(L+1,2*L+1);


%0-TH DEGREE. 
%Step 2 in the iterative resiual fitting (IRF) algorithm. See reference [1].

Y=Y_l(0,theta',varphi')';
Ycommon=inv(Y'*Y)*Y';

betal=Ycommon*x;
betax(1,1)=betal';
xsmooth=Y*betal;
xestimate=xsmooth;

betal=Ycommon*y;
betay(1,1)=betal';
ysmooth=Y*betal;
yestimate=ysmooth;

betal=Ycommon*z;
betaz(1,1)=betal';
zsmooth=Y*betal;
zestimate=zsmooth;


%l-TH DEGREE ITERATION

for l=1:L
 
    % Step 4: residual. See the paper for detail
    x_j = x-xestimate;
    y_j = y-yestimate;
    z_j = z-zestimate;

    Y=Y_l(l,theta',varphi')';
    % real(Y_l) gives negative harmonics imag(Y_l) gives positive harmonics
    
    Y=[real(Y)   imag(Y(:,2:(l+1)))];
    Ycommon=inv(Y'*Y)*Y';
    if(~isempty(find(isnan(Ycommon))))

        break;
    end

    % Step 5: refitting the residual. See the paper for detail
    betal=Ycommon*x_j;
    betax(l+1,1:2*l+1)=betal';
    xsmooth=Y*betal;
    xestimate=xestimate + exp(-l*(l+1)*sigma)*xsmooth;

    betal=Ycommon*y_j;
    betay(l+1,1:2*l+1)=betal';
    ysmooth=Y*betal;
    yestimate=yestimate + exp(-l*(l+1)*sigma)*ysmooth;

    betal=Ycommon*z_j;
    betaz(l+1,1:2*l+1)=betal';
    zsmooth=Y*betal;
    zestimate=zestimate + exp(-l*(l+1)*sigma)*zsmooth;
end;


%output the results in a proper shape
temp=[xestimate; yestimate; zestimate];
surf_smooth.vertices=squeeze(reshape(temp,n_vertex,3));
surf_smooth.faces=surf.faces;

fourier.x=betax;
fourier.y=betay;
fourier.z=betaz;
end

%---------------------------------------------------------------
function [theta,varphi]=EULERangles(surf);

n_vertex=size(surf.vertices,1);
c=mean(surf.vertices);  %mass center
surf.vertices=surf.vertices-kron(ones(n_vertex,1),c);  % translation

[theta,varphi,r] = cart2sph(surf.vertices(:,1),surf.vertices(:,2),surf.vertices(:,3));

% MATLAB coordinate systems are different from the convention used in the
% TMI paper.
temp = theta;
theta = pi/2 - varphi;
varphi = pi + temp;

%figure_wire(surf,'yellow')
end
%-----------------------------------------------------------------
function Y_l=Y_l(l,theta,varphi)
% computes spherical harmonics of degree l.
sz=length(theta);

m=0:l;
CLM=[];
exp_i_m=[];
sign_m=[];
SIGNM=[];
Pn=[];

for k = 0:(2*l)
    fact(k+1) = factorial(k);
end
clm = sqrt(((2*l+1)/(2*pi))*(fact(l-abs(m)+1)./fact(l+abs(m)+1)));
CLM=kron(ones(1,sz),clm');

for k = 0:l
    exp_i_m(k+1,:)= exp(i*k*varphi);
    sign_m(k+1) = (-1)^k;
end
exp_i_m(1,:)=exp_i_m(1,:)/sqrt(2);

SIGNM=kron(ones(1,sz),sign_m');
Pn=legendre(l,cos(theta));
Y_l=CLM.*SIGNM.*Pn.*exp_i_m;
end
