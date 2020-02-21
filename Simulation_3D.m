clear 
close all

rng(13);
%%%%%%%%%%%%%%%%%%
% load the resting shape vertices location S, triangulation ff and parameters P 
load('resting_shape_3D.mat')
[S,P.index,P.index2] = fix_points(S,ff,P);
%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%

P.d_max_int = 1.000e-04; %micro m, equivalent to d_{tol}
P.B_0 =20; %maximum initial number of barbed ends in an initial foci

% Actin dynamics
P.gamma =1; %1/s, corresponding to gamma_cap 
P.tau_1 = 30; %s, corresponding to 1/gamma_uncap
P.tau_2 =1; %s, corresponding to 1/gamma_sever
% P.k_2 = 20;%1/s
P.j =5;%5;%1/s
P.delta_y = 0.0022;%micro M - length of actin monomer


P.a = 3.8000;
% P.k_m_1 = 9.9420;
% P.k_p_1 = 0.1160;
% P.phi = 0.0350;

P.k_bT = 4.1*0.001;%pN micro m
P.k_on = 11.6;%1/microM s barbed-end monomer assembly rate constant

%%% actin force
P.alpha =3.8;%pN
P.sigma = 0.3;
P.cte = 10;%30;%micro/m^2, equivalent to phi
P.tau_prob = 1/40;%micro m, equivalent to lambda 
P.tau_n =1/10;% 1/s, equivalent to gamma_n
P.ini_foci =2;% number of initial foci
% P.ini_foci_0 = P.ini_foci;
P.zeta = 0.004;%micro m^2/s PN, strrength of force update

% Time & spatial mesh

P.d_mem = 0.1; %distance from the membrame of Arp2/3 micro m 
P.d_arp23 = 0.1;%micro m distance of Arp2/3 complex from PSD

P.K = length(S);%number of vertices

P.t_initial = 0;%%s
P.t_end =20*60;% s 
aux_t = P.t_initial:P.delta_t:P.t_end;

%%%%%%%%%%%%%%%%%%%%%%
% For saving data
aux_S = cell(length(aux_t),1); %shape vertices
aux_ff = cell(length(aux_t),1);%triangulation
save_a_points = cell(length(aux_t),1);%location of the nucleation points
save_ini_fil = cell(length(aux_t),1);%corresponding vertices in the membrane
save_line = cell(length(aux_t),1);%V^{i,k}
save_aux_fil = cell(length(aux_t),1);%total number of barbed ends in the spine
save_volume = zeros(length(aux_t),1);%save volume

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
% initiation

% select initial foci
[P.ini_fil,P.a_points,P.line] = initial_B_3D(S,P);

% determine the number of barbed ends at each initial foci
B = randi(P.B_0,length(P.ini_fil),1);
filaments = [];
for l = 1:length(P.ini_fil)
    filaments = [filaments;P.ini_fil(l)*ones(B(l),1) zeros(B(l),1) ones(B(l),1)];
end
filaments = cap_ends(filaments,P);
aux_fil = store_data(filaments,P);


aux_S{1} = S;
aux_ff{1} = ff;
save_a_points{1} = P.a_points;
save_ini_fil{1} = P.ini_fil;
save_line{1} = P.line;
save_aux_fil{1} = aux_fil;
ff_length = [P.a_points aux_t(1)*ones(size(P.ini_fil))];
save_volume(1) = volume_sphere(ff,S);



j=2;

while j <= length(aux_t)
    F_mem = membrane_force_3D(S,ff,P);
    P.w = sqrt(F_mem(:,1).^2 + F_mem(:,2).^2 + F_mem(:,3).^2).*P.delta_y/P.k_bT;
% nucleate new foci
    [ff_length,P.ini_fil,P.a_points,P.line,filaments,B] = nucleation_3D(ff_length,aux_t(j),filaments,B,S,P);
% cap,uncap,remove filaments    
    [aux_new] = create_filaments(filaments,P.a,aux_fil,P);
    [filaments] = remove_old(filaments,P); 
    [filaments] = cap_ends(filaments,P);
    filaments = [filaments; aux_new];
% actualize if a focus dies out    
    aux_n = unique(filaments(:,1));
    if ~isempty(aux_n) 
        aux_nn = [];
        for ll = 1:length(aux_n)
            aux_nn = [aux_nn;find(P.ini_fil == aux_n(ll))];
        end
        P.a_points = P.a_points(aux_nn,:);
        P.ini_fil = P.ini_fil(aux_nn);
        P.line = P.line(aux_nn,:);
        B = B(aux_nn);
    else
        P.a_points = [];
        P.ini_fil = [];
        P.line = [];
        B = [];
    end
    
%    calculate displacement of poins 
    S =  solve_system_threshold_3D_rk(S,ff,aux_fil,1,P);
% remessh 
    [ff,S] = remeshing(ff, S, int32(P.index2), P.delta_S, int32(3));
      
    [S,P.index,P.index2] = fix_points(S,ff,P);

%     save data
    aux_S{j} = S;
    aux_ff{j} = ff;
    save_a_points{j} = P.a_points;
    save_ini_fil{j} = P.ini_fil;
    save_line{j} = P.line;
    save_aux_fil{j} = aux_fil;
    save_volume(j) = volume_sphere(ff,S);

    P.K = length(S);
% update corresponding after remeshing
    [P.ini_fil,filaments] = move_index_3D(S,filaments,P);
    aux_fil = store_data(filaments,P);
    

    j = j + 1;
end

figure(1)
for jj = 1:j
    figure(1)
    clf
    lll = aux_S{jj}; 
    ff_ll = aux_ff{jj};
    trimesh( ff_ll, lll(:,1), lll(:,2), lll(:,3),'edgecolor','k');
    title(['time = ' num2str(aux_t(jj),'%.1f') ' s  '])
    axis([-1.5 1.5 -1.5 1.5 -1.5 1.5])
    set(gca,'fontsize',20);
    xlabel('x (\mu m)')
    ylabel('y (\mu m)')
    zlabel('z (\mu m)')   
end
figure;
plot(aux_t(1:j-1)/60,save_volume(1:j-1),'linewidth',2)
set(gca,'fontsize',20)
xlabel('t (min)')
ylabel('Volume')


aux_ini_fil = zeros(j-1,1);
save_mean_B = zeros(j-1,1);
for ll = 1:j-1
    aux_ini_fil(ll) = length(save_ini_fil{ll});
    save_mean_B(ll) = sum(save_aux_fil{ll})./aux_ini_fil(ll); 
end
save_mean_B(aux_ini_fil==0) = 0;
figure;
plot(aux_t(1:j-1)/60,aux_ini_fil,'linewidth',0.5);
set(gca,'fontsize',20)
xlabel('t (min)')
ylabel('Number of active actin foci')

figure;
plot(aux_t(1:j-1)/60,save_mean_B,'linewidth',0.5);
set(gca,'fontsize',20)
xlabel('t (min)')
ylabel('Mean Barbed ends')

