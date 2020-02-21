function f = f_fil_2D(S,B,P)
% Calculate filament force
    f = zeros(P.K,2);
    if ~isempty(P.a_points) 
        [x1,x2] = meshgrid(S(:,1),S(:,2));
        x_diff = x1' - x1;
        y_diff = x2 - x2';
        d = sqrt(x_diff.^2 + y_diff.^2);
        W = (P.alpha/(P.sigma*sqrt(2*pi)))*exp(-(d.^2)./(2*P.sigma^2));
        B =W*B;
        
        for l=1:size(P.a_points)
            aux_dist = sqrt((S(:,1)-P.a_points(l,1)).^2 + (S(:,2)-P.a_points(l,2)).^2);
            aux_dir = [S(:,1)-P.a_points(l,1) S(:,2)-P.a_points(l,2)];
            f(:,1) = f(:,1) +  B.*(aux_dir(:,1)./aux_dist);
            f(:,2) = f(:,2) +  B.*(aux_dir(:,2)./aux_dist);
        end
    
    end
end