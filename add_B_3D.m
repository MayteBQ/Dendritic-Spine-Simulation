function [ff_length,ini_fil,a_points,line,aux_fil,aux_B] = add_B_3D(ff_length,t,S_0,P)
    tri = delaunayn(S_0);
%     generate ramdom points
    aux2 = max(max(S_0(:,1)),abs(min(S_0(:,1))));
    aux2 = -aux2 + aux2*2*rand(1000,3);
    aux_in_all = aux2(~isnan(tsearchn(S_0,tri, aux2)),:);
%     calculate the distance
    d = abs(sqrt((aux_in_all(:,1)-0).^2 + (aux_in_all(:,2)-0).^2 + (aux_in_all(:,3)-P.h_PSD).^2) - P.r_PSD);
    aux_in_all = aux_in_all(d>=P.d_arp23,:);
    d = d(d>=P.d_arp23);
%     calculate the probability
    aux_d = table(d,aux_in_all);
    aux_d = unique(aux_d);
    aux_d = aux_d{:,:};
    aux_exp = exp(-aux_d(:,1)/P.tau_prob);
    aux_exp2 = ones(size(aux_exp))'*aux_exp;
    aux_d = [aux_exp/aux_exp2 aux_d];

    aux_ind2 = [];
    while isempty(aux_ind2)
        %     choose a nucleation point
        aux_in = rand;
        aux_cumulative = cumsum(aux_d(:,1));
        aux_in = find(aux_cumulative >= aux_in,1);
        if isempty(aux_in)
            aux_in = length(aux_cumulative);
        end
        aux_ind = setdiff(P.index,P.ini_fil);
        aux_d_S = sqrt((S_0(P.index,1)-aux_d(aux_in,3)).^2 + (S_0(P.index,2)-aux_d(aux_in,4)).^2 + (S_0(P.index,3)-aux_d(aux_in,5)).^2 );
        aux_ind2 = find(aux_d_S <= P.d_mem);
    end
  %     choose direction to the membrane and the corresponding vertex
    aux_ini_fil = randi(length(aux_ind2));
    aux_ini_fil = aux_ind(aux_ind2(aux_ini_fil));


    ini_fil = [P.ini_fil;aux_ini_fil];
    aux_l = aux_d(aux_in,3:5);
    a_points = [P.a_points;aux_l];
    m = S_0(aux_ini_fil,:);
    if S_0(aux_ini_fil,1)>=0
        if S_0(aux_ini_fil,2)>=0
            if S_0(aux_ini_fil,3)>=0
                cuad = 1;
            else
                cuad = 2;
            end
        else 
            if S_0(aux_ini_fil,3)>=0
                cuad = 3;
            else
                cuad = 4;
            end
        end
    else 
        if S_0(aux_ini_fil,2)>=0
            if S_0(aux_ini_fil,3)>=0
                cuad = 5;
            else
                cuad = 6;
            end
        else 
            if S_0(aux_ini_fil,3)>=0
                cuad = 7;
            else
                cuad = 8;
            end
        end
    end
    line = [P.line;m cuad];
    aux_B = 1;
    aux_fil = [aux_ini_fil 0 1];
    ff_length = [ff_length; aux_l t];
end
    