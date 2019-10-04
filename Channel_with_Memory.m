function Pr = Channel_with_Memory(numLevel , epsilon , delta)
Pr = zeros(numLevel , numLevel) ; 
n = log2(numLevel) ; 
Pr_z = [(1 - epsilon + delta) / (1 + delta)  , epsilon / (1 + delta) ;
    (1 - epsilon) / (1 + delta)  , (epsilon + delta) / (1 + delta)] ;
for x = 1 : numLevel
    for y = 1 : numLevel
        binary_x = de2bi(x - 1 , n , 'left-msb') ;
        binary_y = de2bi(y - 1 , n , 'left-msb') ;
        binary_z = xor(binary_x , binary_y) ;
        
        if (binary_z (1) == 1) 
            Probability = epsilon ; 
        else
            Probability = 1 - epsilon ; 
        end 
        
        for i = 2 : n
            Probability = Probability * Pr_z(binary_z(i - 1) + 1 , binary_z(i) + 1) ; 
        end 
        Pr(x , y) = Probability ; 
    end
end
end