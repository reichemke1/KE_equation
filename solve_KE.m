function [ stream_all] = solve_KE(lat,level,t_mean,pot_mean,v_mean,uv_eddy,vt_eddy,tdt_moist_mean,tdt_rad_mean,red_lat)

% if you use this script in your research please cite one/more of our studies who used
% this script as well:

%Chemke R. (2022) Npj climate and atmospheric science
%Chemke R. (2021) Geophysical Research Letters
%Chemke R. & Polvani L. M. (2021) Geophysical Research Letters
%Chemke R. & Polvani L. M. (2019) Nature Geoscience
%Chemke R., Polvani L. M. & Deser C. (2019) Geophysical Research Letters
%Chemke R. & Polvani L. M. (2018) Geophysical Research Letters

% This function calculates the solution for the KE equation. To get the
% relative contribution of each term in the eqaution to a change in the streamfunction one has to design
% other scripts.

%it signifcantly reduces the calculation time if you put single fields and not double

%this fucntion includes local calls to merdional and vertcial gradients on a sphere: mrdnl_gradient_yz, vrtcl_gradient_yz.

%input parameters:

%lat - latitude
%level - pressure
%t_mean - zonal mean (lat,level)
%pot_mean - zonal mean potenital temperature (lat,level)
%v_mean - zonal mean meridional wind (lat,level)
%uv_eddy - zonal mean eddy momentum flux (lat,level); u'v'
%vt_eddy - zonal mean eddy heat flux (lat,level); v'T'
%tdt_moist_mean - zonal mean latent heating (lat,level) in K/s
%tdt_rad_mean - zonal mean radiative heating (lat,level) in K/s
%red_lat - number of latitudes to remove from the pole.


omega=7.292e-5;
radius=6371000;
Rgas=287.04;
g=9.81;
% cp=1004;

    
cosd_yz_v=repmat(cosd(lat),[1 length(level)]);
rho=repmat(level',[length(lat) 1])./(Rgas.*t_mean);
S2=-(1./(rho.*pot_mean)).*vrtcl_gradient_yz(pot_mean,level);
vt_t=(1./cosd_yz_v).*mrdnl_gradient_yz(vt_eddy.*cosd_yz_v,lat,radius);
emfd=(1./cosd_yz_v.^2).*mrdnl_gradient_yz(uv_eddy.*cosd_yz_v.^2,lat,radius);
eddy_heat=-(Rgas./(repmat(level',[length(lat) 1]))).*mrdnl_gradient_yz(vt_t,lat,radius);
eddy_mom=2.*omega.*repmat(sind(lat),[1 length(level)]).*vrtcl_gradient_yz(emfd,level);
tdt_rad_mean_y=(Rgas./(repmat(level',[length(lat) 1]))).*mrdnl_gradient_yz(tdt_rad_mean,lat,radius);
tdt_moist_mean_y=(Rgas./(repmat(level',[length(lat) 1]))).*mrdnl_gradient_yz(tdt_moist_mean,lat,radius);
%
cori_term=2.*omega.*repmat(sind(lat),[1 length(level)]).*v_mean;
fric1=emfd-cori_term;
fric=-2.*omega.*repmat(sind(lat),[1 length(level)]).*vrtcl_gradient_yz(fric1,level);
%


f=2*omega*sind(lat);
y=(pi./180).*lat.*radius;
dz1=(level(2:end)+level(1:end-1))./2;
dz2=(dz1(2:end)-dz1(1:end-1));
dz=(level(2:end)-level(1:end-1));
dy1=(y(2:end)+y(1:end-1))./2;
dy2=(dy1(2:end)-dy1(1:end-1));
dy=(y(2:end)-y(1:end-1)).*(cosd((lat(2:end)+lat(1:end-1))./2));

ver=repmat(g.*f.^2./(2*pi*radius.*cosd(lat)),[1 length(level)]);
hor=repmat(g./(2*pi*radius),[length(lat) length(level)]).*S2;


z1=(level(2)-level(1));
z2=(level(end)-level(end-1));
y1=(y(2)-y(1)).*(cosd((lat(2)+lat(1))./2));
y2=(y(end)-y(end-1)).*(cosd((lat(end)+lat(end-1))./2));


h1=1;
h2=length(level);

l1=red_lat;
l2=length(lat)-red_lat+1;

%on the main diagonal you put the operator for each level, lat after lat, on
%its sides you put the derivative in level. on the sides blocks (main diagonal +- length(level)) you put the
%derivate in lat (these will be just main diagonal). for the vector on the
%RHS just reshape all levels lat after lat. then solve.
k=0;
for l=l1:l2
for n=h1:h2
    k=k+1;   
    %start with the main diagonal and vert derv.
    if n==h1
        if l==l1
        A(k,k)=(1./z1).*(-1./dz(n)-1./z1).*ver(l,n)+...
               (1./y1).*(-1./dy(l)-1./y1).*hor(l,n);
         
        elseif l==l2
        A(k,k)=(1./z1).*(-1./dz(n)-1./z1).*ver(l,n)+...
               (1./y2).*(-1./dy(l-1)-1./y2).*hor(l,n);
        elseif l>l1 && l<l2
        A(k,k)=(1./z1).*(-1./dz(n)-1./z1).*ver(l,n)+...
               (1./dy2(l-1)).*(-1./dy(l-1)-1./dy(l)).*hor(l,n);
        end
        
        A(k,k+1)=(1./z1).*(1./dz(n)).*ver(l,n);  
        
    elseif n==h2
        if l==l1
         A(k,k)=(1./z2).*(-1./dz(n-1)-1./z2).*ver(l,n)+...
               (1./y1).*(-1./dy(l)-1./y1).*hor(l,n);
        elseif l==l2
        A(k,k)=(1./z2).*(-1./dz(n-1)-1./z2).*ver(l,n)+...
               (1./y2).*(-1./dy(l-1)-1./y2).*hor(l,n);
        elseif l>l1 && l<l2
        A(k,k)=(1./z2).*(-1./dz(n-1)-1./z2).*ver(l,n)+...
               (1./dy2(l-1)).*(-1./dy(l-1)-1./dy(l)).*hor(l,n);
        end
        
        A(k,k-1)=(1./z2).*(1./dz(n-1)).*ver(l,n);
        
    elseif n>h1 && n<h2
    if l==l1
    A(k,k)=(1./dz2(n-1)).*(-1./dz(n-1)-1./dz(n)).*ver(l,n)+...
           (1./y1).*(-1./dy(l)-1./y1).*hor(l,n);    
    elseif l==l2
    A(k,k)=(1./dz2(n-1)).*(-1./dz(n-1)-1./dz(n)).*ver(l,n)+...
           (1./y2).*(-1./dy(l-1)-1./y2).*hor(l,n); 
    elseif l>l1 && l<l2
    A(k,k)=(1./dz2(n-1)).*(-1./dz(n-1)-1./dz(n)).*ver(l,n)+...
           (1./dy2(l-1)).*(-1./dy(l-1)-1./dy(l)).*hor(l,n);
    end
    A(k,k-1)=(1./dz2(n-1)).*(1./dz(n-1)).*ver(l,n);
    A(k,k+1)=(1./dz2(n-1)).*(1./dz(n)).*ver(l,n);        
    
    end
    
    
    if l==l1
    A(k,k+length(h1:h2))=(1./y1).*(1./dy(l)).*hor(l,n);
    elseif l==l2
    A(k,k-length(h1:h2))=(1./y2).*(1./dy(l-1)).*hor(l,n);
    elseif l>l1 && l<l2
    A(k,k-length(h1:h2))=(1./dy2(l-1)).*(1./dy(l-1)).*hor(l,n);
    A(k,k+length(h1:h2))=(1./dy2(l-1)).*(1./dy(l)).*hor(l,n);
    end
end    
end


D=(eddy_heat+eddy_mom+tdt_rad_mean_y+tdt_moist_mean_y+fric);
D=reshape(D(l1:l2,h1:h2).',length(lat(l1:l2))*length(level(h1:h2)),1);
stream=A\D;
[ x, err, iter, flag ] = sor( A, stream,D,.5, 1e3, 1e-4 );% 1 is guass siedel
stream_all1=(reshape(x,length(h1:h2),length(lat(l1:l2))))';
stream_all=nan(length(lat),length(level));
stream_all(l1:l2,h1:h2)=stream_all1;



end

