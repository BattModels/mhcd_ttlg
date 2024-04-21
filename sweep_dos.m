%% Load Data
% Load angles
%load('q_dict.mat')
% Load k-rate data
clear all
cd '/Users/mbabar/Desktop/PhD/Analysis/MHCD_tTLG'
load('Eo_var/_0.07/k_data_1.0_?_0.82_?_0.0.mat') % Rates

% Load dos data
dir0 = 'sweep/';
end_str = '_kcut_6.9282_qtype_1_nq_484_zip.mat';

fd_list = dir(append(dir0));

% Given angles, load dos data
dos_max=zeros(1,length(q12_list));
for j=1:length(q12_list)
    theta12 = q12_list(j);
    theta23 = q23_list(j);
    f_name = append('dos_q12_',string(theta12),'_q23_',string(theta23),'_');
    for i=1:length(fd_list)
        str = fd_list(i).name;
        if contains(str, f_name)
            load(append(dir0,str))
            break
        end
    end
    ind_range = find(E_list>-0.5 & E_list<0.5);
    %dos_max(j) = max(abs(dos_tot));
    dos_max(j) = max(abs(dos_tot(ind_range)));
end

%load(append(dir0,'dos_q12_',string(theta12),'_q23_',string(theta23),end_str))


%% Interpolate and plot

x = q12_list;
y = q23_list;
%v = kox_list;
v = kred_list;
%v = dos_max;

k_aba = 2.6323040539732973e-5; % Untwisted ABA Bernal stacked eqb. ko 

[xq,yq] = meshgrid(1:.01:5, 1:.01:5);

% Duplicate values
[theta,uid]  = unique([x(:), y(:)], 'rows');
% Non-unique elements
nuid = 1:1:length(x);
nuid(uid) = [];
% Check non-unique elements
% q12_list(nuid)
% q23_list(nuid)

% Avoid duplicate values
ux = theta(:,1);
uy = theta(:,2);
uv = v(uid);
vq = griddata(ux,uy,uv,xq,yq, 'cubic');

%vq = griddata(x,y,v,xq,yq,'cubic'); % To include duplicate points averaged

surf(xq,yq,(vq+transpose(vq))/(2*k_aba))
shading flat
xlim([1 5])
ylim([1 5])
colorbar
set(gca,'FontName','Times New Roman','FontSize',20,'LineWidth',3,'ytick',1:5,'xtick',1:5);
xlabel('$\theta_{12}$($^\circ$)','interpreter','latex')
ylabel('$\theta_{23}$($^\circ$)','interpreter','latex')
%title('Rate k (cm s$^{-1}$) Ruhex, Vapp=0.2 V','interpreter','latex')
title('$k_{ox} (\eta = 0.0)$','interpreter','latex')
colormap(jet)
view(2)
saveas(gcf,'kred_0.0.png')
savefig('kred_0.0.fig')