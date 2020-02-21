function [S_new,index,index2,index3,index_psd,index_n] = remesh_2D(S,P)
  

    
    new_S = S(1:P.index2(1),:);
    aux = size(new_S,1);
    i = 1;
    while  i <= aux - 1
        dd = sqrt((new_S(i+1,1)-new_S(i,1))^2 + (new_S(i+1,2)-new_S(i,2))^2);
        if dd < P.low
            new_S(i+1,:) = [];
            if i+1 < size(new_S,1)
                dd = sqrt((new_S(i+1,1)-new_S(i,1))^2 + (new_S(i+1,2)-new_S(i,2))^2); 
                while dd > P.high
                    new_S = [new_S(1:i,:);
                        [(new_S(i+1,1)+new_S(i,1))/2,(new_S(i+1,2)+new_S(i,2))/2];
                        new_S(i+1:end,:)];
                    dd = sqrt((new_S(i+1,1)-new_S(i,1))^2 + (new_S(i+1,2)-new_S(i,2))^2);  
                end
            end
        elseif dd > P.high
            while dd > P.high
                new_S = [new_S(1:i,:);
                    [(new_S(i+1,1)+new_S(i,1))/2,(new_S(i+1,2)+new_S(i,2))/2];
                    new_S(i+1:end,:)];
                dd = sqrt((new_S(i+1,1)-new_S(i,1))^2 + (new_S(i+1,2)-new_S(i,2))^2);  
            end
        end         
        aux = size(new_S,1);
        i = i+1;
    end
    if i==aux
        dd = sqrt((S(P.index2(1),1)-new_S(i,1))^2 + (S(P.index2(1),2)-new_S(i,2))^2);
        if dd < P.low
            new_S(i,:) = [];
            dd = sqrt((S(P.index2(1),1)-new_S(end,1))^2 + (S(P.index2(1),2)-new_S(end,2))^2); 
        end
        if dd > P.high 
            while dd > P.high
                new_S = [new_S;
                [(S(P.index2(1),1)+new_S(end,1))/2,(S(P.index2(1),2)+new_S(end,2))/2]];
                dd = sqrt((S(P.index2(1),1)-new_S(end,1))^2 + (S(P.index2(1),2)-new_S(end,2))^2);  
            end
        end
    end
        
   
    S_new = [new_S;S(P.index2(1),:);S(P.index_psd,:);S(P.index2(2),:)];
    index = 1:size(new_S,1); 
    index2 = size(new_S,1)+1;
    index2 = [index2; (index2+length(P.index_psd)+1)];
    index_psd =(index2(1)+1):(index2(2)-1);

    
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  
    new_S = S((P.index2(2)+1):P.index3(1),:);
    dd = sqrt((new_S(1,1)-S(P.index2(2),1))^2 + (new_S(1,2)-S(P.index2(2),2))^2);
    if dd < P.low
        new_S(1,:) = [];
        dd = sqrt((new_S(1,1)-S(P.index2(2),1))^2 + (new_S(1,2)-S(P.index2(2),2))^2); 
        while dd > P.high
            new_S = [[(new_S(1,1)+S(P.index2(2),1))/2,(new_S(1,2)+S(P.index2(2),2))/2];
                new_S(1:end,:)];
            dd = sqrt((new_S(1,1)-S(P.index2(2),1))^2 + (new_S(1,2)-S(P.index2(2),2))^2); 
        end
    elseif dd > P.high
        while dd > P.high
            new_S = [[(new_S(1,1)+S(P.index2(2),1))/2,(new_S(1,2)+S(P.index2(2),2))/2];
                new_S(1:end,:)];
            dd = sqrt((new_S(1,1)-S(P.index2(2),1))^2 + (new_S(1,2)-S(P.index2(2),2))^2);  
        end
    end         
    aux = size(new_S,1);
    
    i = 1;
     while  i <= aux - 1
        dd = sqrt((new_S(i+1,1)-new_S(i,1))^2 + (new_S(i+1,2)-new_S(i,2))^2);
        if dd < P.low
            new_S(i+1,:) = [];
            if i+1 < size(new_S,1)
                dd = sqrt((new_S(i+1,1)-new_S(i,1))^2 + (new_S(i+1,2)-new_S(i,2))^2); 
                while dd > P.high
                    new_S = [new_S(1:i,:);
                        [(new_S(i+1,1)+new_S(i,1))/2,(new_S(i+1,2)+new_S(i,2))/2];
                        new_S(i+1:end,:)];
                    dd = sqrt((new_S(i+1,1)-new_S(i,1))^2 + (new_S(i+1,2)-new_S(i,2))^2);  
                end
            end
        elseif dd > P.high
            while dd > P.high
                new_S = [new_S(1:i,:);
                    [(new_S(i+1,1)+new_S(i,1))/2,(new_S(i+1,2)+new_S(i,2))/2];
                    new_S(i+1:end,:)];
                dd = sqrt((new_S(i+1,1)-new_S(i,1))^2 + (new_S(i+1,2)-new_S(i,2))^2);  
            end
        end         
        aux = size(new_S,1);
        i = i+1;
    end
    if i==aux
        dd = sqrt((S(P.index3(1),1)-new_S(i,1))^2 + (S(P.index3(1),2)-new_S(i,2))^2);
        if dd < P.low
            new_S(i,:) = [];
            dd = sqrt((S(P.index3(1),1)-new_S(end,1))^2 + (S(P.index3(1),2)-new_S(end,2))^2); 
        end
        if dd > P.high 
            while dd > P.high
                new_S = [new_S;
                [(S(P.index3(1),1)+new_S(end,1))/2,(S(P.index3(1),2)+new_S(end,2))/2]];
                dd = sqrt((S(P.index3(1),1)-new_S(end,1))^2 + (S(P.index3(1),2)-new_S(end,2))^2);  
            end
        end
    end
      
    S_new = [S_new;new_S; S(P.index3(1),:); S(P.index_n,:); S(P.index3(2),:)];
    index = [index (index2(end)+1):(index2(end)+size(new_S,1))];
    index3 = index2(end)+size(new_S,1) + [1 length(P.index_n)+2];
    index_n =  (index2(end)+size(new_S,1)+2):(index2(end)+length(P.index_n)+size(new_S,1)+1);
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    new_S = S((P.index3(2)+1):end,:);
    dd = sqrt((new_S(1,1)-S(P.index3(2),1))^2 + (new_S(1,2)-S(P.index3(2),2))^2);
    if dd < P.low
        new_S(1,:) = [];
        dd = sqrt((new_S(1,1)-S(P.index3(2),1))^2 + (new_S(1,2)-S(P.index3(2),2))^2); 
        while dd > P.high
            new_S = [[(new_S(1,1)+S(P.index3(2),1))/2,(new_S(1,2)+S(P.index3(2),2))/2];
                new_S(1:end,:)];
            dd = sqrt((new_S(1,1)-S(P.index3(2),1))^2 + (new_S(1,2)-S(P.index3(2),2))^2); 
        end
    elseif dd > P.high
        while dd > P.high
            new_S = [[(new_S(1,1)+S(P.index3(2),1))/2,(new_S(1,2)+S(P.index3(2),2))/2];
                new_S(1:end,:)];
            dd = sqrt((new_S(1,1)-S(P.index3(2),1))^2 + (new_S(1,2)-S(P.index3(2),2))^2);  
        end
    end         
    
    aux = size(new_S,1);
    i = 1;
    
     while  i <= aux - 1
        dd = sqrt((new_S(i+1,1)-new_S(i,1))^2 + (new_S(i+1,2)-new_S(i,2))^2);
        if dd < P.low
            new_S(i+1,:) = [];
            if i+1 < size(new_S,1)
                dd = sqrt((new_S(i+1,1)-new_S(i,1))^2 + (new_S(i+1,2)-new_S(i,2))^2); 
                while dd > P.high
                    new_S = [new_S(1:i,:);
                        [(new_S(i+1,1)+new_S(i,1))/2,(new_S(i+1,2)+new_S(i,2))/2];
                        new_S(i+1:end,:)];
                    dd = sqrt((new_S(i+1,1)-new_S(i,1))^2 + (new_S(i+1,2)-new_S(i,2))^2);  
                end
            end
        elseif dd > P.high
            while dd > P.high
                new_S = [new_S(1:i,:);
                    [(new_S(i+1,1)+new_S(i,1))/2,(new_S(i+1,2)+new_S(i,2))/2];
                    new_S(i+1:end,:)];
                dd = sqrt((new_S(i+1,1)-new_S(i,1))^2 + (new_S(i+1,2)-new_S(i,2))^2);  
            end
        end         
        aux = size(new_S,1);
        i = i+1;
    end
    if i==aux
        dd = sqrt((S(1,1)-new_S(i,1))^2 + (S(1,2)-new_S(i,2))^2);
        if dd < P.low
            new_S(i,:) = [];
            dd = sqrt((S(1,1)-new_S(end,1))^2 + (S(1,2)-new_S(end,2))^2); 
        end
        if dd > P.high 
            while dd > P.high
                new_S = [new_S;
                [(S(1,1)+new_S(end,1))/2,(S(1,2)+new_S(end,2))/2]];
                dd = sqrt((S(1,1)-new_S(end,1))^2 + (S(1,2)-new_S(end,2))^2);  
            end
        end
    end
    
    S_new = [S_new;new_S];
    index = [index (index3(end)+1):(index3(end)+size(new_S,1))];

end