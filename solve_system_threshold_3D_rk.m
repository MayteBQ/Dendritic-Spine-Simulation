function S_int = solve_system_threshold_3D_rk(S,ff,aux_fil,actin_dyn,P)

    delta_t = P.delta_t;
    count = 0;

    S_int = solve_system_RK_3D(S,ff,aux_fil,actin_dyn,delta_t,P);
    S_int(P.index2,:) = S(P.index2,:);
    d_max = max(sqrt((S_int(:,1)-S(:,1)).^2 + (S_int(:,2)-S(:,2)).^2 + (S_int(:,3)-S(:,3)).^2));
%     check if the diplacement of the vertices is less than d_tol,
%     otherwise half the time interval
    while d_max >  P.d_max_int
        delta_t = delta_t/2;
        S_int = solve_system_RK_3D(S,ff,aux_fil,actin_dyn,delta_t,P);
        S_int(P.index2,:) = S(P.index2,:);
        d_max = max(sqrt((S_int(:,1)-S(:,1)).^2 + (S_int(:,2)-S(:,2)).^2 + (S_int(:,3)-S(:,3)).^2));
        count = count+1;
    end
 
    if count > 0
        for j=1:(2^count-1)
            S_int_new = solve_system_RK_3D(S_int,ff,aux_fil,actin_dyn,delta_t,P);
            S_int_new(P.index2,:) = S_int(P.index2,:);
            S_int = S_int_new;
            
        end
    end


end


    
    
        