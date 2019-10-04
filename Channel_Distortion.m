function D_c = Channel_Distortion (f , T , numLevel , delta , Pr , codebook , b)
summation = 0 ; 

for i = 1 : numLevel
    for j = 1 : numLevel
        u_index = find (T (: , 2) == i) ;
        P_i = delta * sum (f(u_index)) ;
        summation = summation + Pr(j , b(i)) * P_i * (codebook (b(i)) - codebook(j)) ^ 2;
    end
end
D_c = summation ;
end