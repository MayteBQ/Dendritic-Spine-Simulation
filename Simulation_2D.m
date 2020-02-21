
load('resting_shape_2D.mat')
% load S matrix containing resting shape and P an structure with the
% parameters (some change name, see paper and Simulation 3D)

P.t_end =15*60;% s 
aux_t = P.t_initial:P.delta_t:P.t_end;

% For saving data 
aux_S = cell(length(aux_t),1);%spine shape
aux_area = zeros(length(aux_t),1);%spine area
aux_foci = zeros(length(aux_t),1);%number of foci
aux_B = zeros(length(aux_t),1);%number of barbed ends
filaments = [];



% choose the initial nucleation points 
[P.ini_fil,P.a_points,P.line] = initial_B_2D(S,P);

% assignate a number of barbed ends
P.K = length(S);
B = randi(P.B_0,length(P.ini_fil),1);
for l = 1:length(P.ini_fil)
    filaments = [filaments;P.ini_fil(l)*ones(B(l),1) zeros(B(l),1) ones(B(l),1)];
end
filaments = cap_ends(filaments,P);


ff_length = [P.a_points aux_t(1)*ones(size(P.ini_fil))];%calculate foci lifetime
aux_foci(1) = length(P.ini_fil);
aux_area(1) =area_S(S);
aux_S{1} = S;
aux_fil = store_data(filaments,P);
aux_B(1) = sum(aux_fil);
 

j=2;
while j<=length(aux_t)
    F_mem = force_membrane_2D(S,P);
    P.w = sqrt(F_mem(:,1).^2 + F_mem(:,2).^2).*P.delta_y/P.k_bT;
% Actin related events
    [ff_length,P.ini_fil,P.a_points,P.line,filaments] = nucleation_2D(aux_t(j),ff_length,filaments,S,P);

    [aux_new] = create_filaments(filaments,P.a,aux_fil,P);
    [filaments] = remove_old(filaments,P); 
    [filaments] = cap_ends(filaments,P);
    filaments = [filaments; aux_new];
%     check if a foci died out and update the corresponding index of the
%     corresponding vertex in the membrane and the lifetime of that foci
    aux_n = unique(filaments(:,1));
    aux_dif = 1:length(P.ini_fil);
    if ~isempty(aux_n) 
        aux_nn = [];
        for ll = 1:length(aux_n)
            aux_nn = [aux_nn;find(P.ini_fil == aux_n(ll))];
        end
        C = setdiff(aux_dif,aux_nn)';
        ff_length = [ff_length;P.a_points(C,:) aux_t(j)*ones(size(C))];
        P.a_points = P.a_points(aux_nn,:);
        P.ini_fil = P.ini_fil(aux_nn);
        P.line = P.line(aux_nn,:);

    else
        ff_length = [ff_length;P.a_points aux_t(j)*ones(length(P.ini_fil),1)];
        P.a_points = [];
        P.ini_fil = [];
        P.line = [];
    end

    S = solve_system_threshold_RK(S,aux_fil,1,P);
    [S_new,P] = remesh_2D_iter(S,P);

    P.K = length(S_new);
    [P.ini_fil,filaments] = move_index_2D(S_new,filaments,P);
    S = S_new;
    aux_fil = store_data(filaments,P);

    aux_area(j) = area_S(S);
    aux_foci(j) = length(P.ini_fil);
    aux_B(j) = sum(aux_fil);
    aux_S{j} = S;

    j = j+1;
end

for l=1:(10/P.delta_t):j-1
    figure(100);clf
    set(gca,'fontsize',20)
    hold on 
    lll = aux_S{l};  
    plot([lll(:,1); lll(1,1)],[lll(:,2); lll(1,2)],'k')
    title(['t = ' num2str(floor(aux_t(l)/60),'%.0f') ' min  ' num2str(rem(aux_t(l),60),'%.0f') ' s'])
    axis([-1.5 1.5 -1.5 1.5]) 
    set(gcf,'color','w');
end


figure
plot(aux_t/60,aux_area)
ylabel('area [{\mu}m^2]')
xlabel('time [min]')

figure
plot(aux_t/60,aux_foci)
ylabel('Number of Foci')
xlabel('time [min]')

figure
mean_B = aux_B./aux_foci;
mean_B(aux_foci == 0) = 0;
plot(aux_t/60,mean_B)
ylabel('mean Barbed ends')
xlabel('time [min]')



aux = ff_length;
time = [];
while ~isempty(aux)
    A = (aux(:,1) == aux(1,1));
    B = (aux(:,2) == aux(1,2));
    C = A.*B;
    D = aux(C==1,:);
    if sum(C) > 1
        time = [time; D(1,:) D(2,3)];
        aux = aux(C==0,:);
    elseif sum(C) == 1
        time = [time; D(1,:) aux_t(end)];
        aux = aux(C==0,:);
    end
end
auxx = (time(:,4)-time(:,3));
figure
histogram(auxx)
xlabel('lifetime [s]')
ylabel('foci')
set(gca,'fontsize',20)