function f = f_fil_3D(S,B,P)
% calculating the barbed ends force
    f = zeros(P.K,3);
    if ~isempty(P.a_points)
        x_diff = S(:,1)-S(:,1)';
        y_diff = S(:,2)-S(:,2)';
        z_diff = S(:,3)-S(:,3)';
        d = sqrt(x_diff.^2 + y_diff.^2 + z_diff.^2);
        W = (P.alpha/(P.sigma*sqrt(2*pi)))*exp(-(d.^2)./(2*P.sigma^2));
        B =W*B;
        
        for l=1:size(P.a_points,1)
            aux_dist = sqrt((S(:,1)-P.a_points(l,1)).^2 + (S(:,2)-P.a_points(l,2)).^2 + (S(:,3)-P.a_points(l,3)).^2);
            aux_dir = [S(:,1)-P.a_points(l,1) S(:,2)-P.a_points(l,2) S(:,3)-P.a_points(l,3)];
            f = f +  B.*(aux_dir./aux_dist);
        end
    
    end





end