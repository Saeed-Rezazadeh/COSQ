function D = COSQ_Distortion (f , numLevel , codebook , delta , Pr , T , b) 
summation = 0; 
parfor i = 1 : numLevel
    for j = 1 : numLevel
        u_index = find (T (: , 2) == i) ; 
        u = T(u_index , 1) ; 
    
        summation = summation + Pr (b(i) , j) * delta * sum (f(u_index) .* (u - codebook (j)) .^ 2) ;   
    end 
end 
D = summation ; 
end 