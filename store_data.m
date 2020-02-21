function f = store_data(filaments,P)
% store the number of barbed ends in each foci, labelled by the
% corresponding vertex
    f = zeros(P.K,1);
    for j=1:P.K
        aux = find(filaments(:,1)==j);
        if ~isempty(aux)
            f(j) = sum(filaments(aux,2)==0);
        end
    end 
end