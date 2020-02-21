function [ff_length,ini_fil,a_points,line,filaments,B] = nucleation_3D(ff_length,t,filaments,B,S,P)
    ini_fil = P.ini_fil;
    a_points = P.a_points;
    line = P.line;

    if rand < P.delta_t*P.tau_n

        [ff_length,ini_fil,a_points,line,aux_fil,aux_B] = add_B_3D(ff_length,t,S,P);
        filaments = [filaments;aux_fil];
        B = [B;aux_B];
        
    end
end