function [S_new,P] = remesh_2D_iter(S,P)

    S_new = S;
    for h=1:P.iter
        [S_new,P.index,P.index2,P.index3,P.index_psd,P.index_n] = remesh_2D(S_new,P);
    end
end