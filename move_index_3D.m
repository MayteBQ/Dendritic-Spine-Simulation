function [aux_ini_fil,filaments] = move_index_3D(S_new,filaments,P)
% after remeshing update the closest mesh vertex to the intersection
% with the line from the nucleation location


    aux_ini_fil = zeros(size(P.ini_fil));

    
    for l = 1:length(P.ini_fil)
        if P.line(l,4) == 1
            ind = find(S_new(P.index,1) >= 0 & S_new(P.index,2) >= 0  & S_new(P.index,3) >= 0);
        elseif P.line(l,4) == 2
             ind = find(S_new(P.index,1) >= 0 & S_new(P.index,2) >= 0  & S_new(P.index,3) < 0);
        elseif P.line(l,4) == 3
             ind = find(S_new(P.index,1) >= 0 & S_new(P.index,2) < 0  & S_new(P.index,3) >= 0);
        elseif P.line(l,4) == 4
             ind = find(S_new(P.index,1) >= 0 & S_new(P.index,2) < 0  & S_new(P.index,3) < 0);
        elseif P.line(l,4) == 5
             ind = find(S_new(P.index,1) < 0 & S_new(P.index,2) >= 0  & S_new(P.index,3) >= 0);
        elseif P.line(l,4) == 6
             ind = find(S_new(P.index,1) < 0 & S_new(P.index,2) >= 0  & S_new(P.index,3) < 0);
        elseif P.line(l,4) == 7
             ind = find(S_new(P.index,1) < 0 & S_new(P.index,2) < 0  & S_new(P.index,3) >= 0);
        else
            ind = find(S_new(P.index,1) < 0 & S_new(P.index,2) < 0  & S_new(P.index,3) < 0);
        end
        nn=size(ind);
        aux1 = cross(repmat(P.line(l,1:3)-P.a_points(l,:),nn(1),1),P.a_points(l,:)-S_new(P.index(ind),:));
        aux1 = sqrt(diag(aux1*aux1'));
        aux2 = P.line(l,1:3)-P.a_points(l,:);
        aux2 = sqrt(aux2*aux2');
        aux_d = aux1/aux2;
        aux_ini_fil(l) = P.index(ind(aux_d == min(aux_d)));
        filaments(filaments(:,1)==P.ini_fil(l),1) = aux_ini_fil(l);
    end


end