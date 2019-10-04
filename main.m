clc ;
clear all;
close all

%% The Results.txt
% This file contains the ultimate SDR values for different channel parameters delta and epsilon.
% Also, since the ACOSQ is an iterative algorithm, the distirtion
% value for every iteration is provided in this file for given a epsilon and delta.

FileID = fopen ('Results.txt' , 'a') ;
%% Choose the number of quantozation levels
QuantizationLevels = [2 4 8 16] ;

mu = 0 ; % the source's mean value
sigma = 1 ; % the source's standard devision.

alpha = 500000  ; % The size of the training set to choose the initial codebook.


delta = 5 ; % determine the channel noise correlation. Choose between 0 5 and 10. 

%% Channel's cross-over probability epsilon
epsilon = unique ([10 ^ -6 10^-5 : 2 * 10^-5 : 10^-4 , 10 ^ -4 : 10^-4  : 10^-3 ,10 ^ -3 ,  0.005 0.01 0.05 0.1]);

% Since the design converges to a locally optimal solution, to avoid
% bad local optimums, we use a increase-decrease method.
SIZE = length(epsilon) ;
noise =  [1 : SIZE , SIZE : -1 : 1  , 1 : SIZE , SIZE : -1 : 1 ] ;

% The variable resolution determines the accuracy of the Riemann summation
resolution = 2 ^ 11 ;

% Initialize the parameters
SDR = zeros (SIZE , 1) ;
Total_SDR = zeros (length (epsilon) , 1) ;

[indexed_T , delta_u] = initialization (sigma , mu , alpha , resolution) ;
% Compute the source pdf. We herein consider a zero-mean
% unit-variance Gaussian source distribution.
u = indexed_T(: , 1) ;
f =  1 ./ (sqrt (2 .* pi)) .* exp (-u .^ 2 ./ 2) ;
f = f./(sum(f) * delta_u) ;
for h = 1 : length(QuantizationLevels)
    
    numLevel = QuantizationLevels (h);
    
    % Use splitting algorithm to choose the initial cdebook. 
    codebook = init_codebook(numLevel , f , delta_u , indexed_T , alpha) ;
    
    Pr = Channel_with_Memory(numLevel ,  0 , delta);
    [~, D_s , indexed_T , codebook] = ...
        COSQ (f , Pr , numLevel , indexed_T , codebook , delta_u , [1 : numLevel]) ;
    
    Pr = Channel_with_Memory(numLevel ,  epsilon(1) , delta);
    %% Simulated Annealing process
    mapping = randperm (numLevel);
    b = mapping ;
    T_0 = 10 ;
    T_f = 0.00025 ;
    a = 0.97 ;
    N_cut = 200 ;
    
    
    D_c = Channel_Distortion (f , indexed_T , numLevel , delta_u , Pr , codebook , b ) ;
    iterations = 0 ;
    perturbations = 0 ;
    sucessful = 0 ;
    delta_D_c = 1 ;
    Threshold = 0.0001 ;
    while ( T_0 > T_f && perturbations < 50000)
        rand_index = randperm (numLevel , 2) ;
        temp_b = b ;
        temp_b([rand_index(1) rand_index(2)]) = temp_b([rand_index(2) rand_index(1)]) ;
        b_prime = temp_b ;
        
        D_c_prime = Channel_Distortion (f , indexed_T , numLevel , delta_u , Pr , codebook , b_prime ) ;
        
        delta_D_c = D_c_prime - D_c ;
        acceptance_Pr = exp (-delta_D_c / T_0) ;
        if (delta_D_c <= 0)
            b = b_prime ;
            D_c = D_c_prime ;
            
            sucessful = sucessful + 1 ;
            perturbations = 0 ;
        elseif (rand <= acceptance_Pr)
            b = b_prime ;
            D_c = D_c_prime ;
            
            perturbations = perturbations + 1 ;
        else
            perturbations = perturbations + 1 ;
        end
        if (iterations >= N_cut || sucessful >= 5 )
            T_0 = a * T_0 ;
            iterations = 0 ;
            sucessful = 0 ;
        end
        iterations = iterations + 1 ;
    end
    
    fprintf (FileID , 'End of SA algorithm\n') ;
    codebook = codebook (b) ;
    fprintf (FileID , 'b = %d\n' , b) ;
    for  k = 1 : length(noise)
        
        i = noise(k) ;
        
        Pr = Channel_with_Memory(numLevel ,  epsilon(i) , delta);
        [SDR(k), D , indexed_T, codebook] = COSQ (f , Pr , numLevel , indexed_T , codebook , delta_u , b ) ;
        sdr (h , k) = SDR (k) ;
        save('codebook' , 'codebook') ;
        Data = ['Data\Data_' num2str(k) '_numLevel_' num2str(numLevel) '_delta_' num2str(delta)] ;
        save(Data , 'indexed_T' , 'codebook' , 'delta' , 'b') ;
    end
    fprintf (FileID , 'numLevel = %d\n' , numLevel) ;
end
clear D ;
for h = 1 : length(QuantizationLevels)
    for i = 1 : SIZE
        index = find (noise == i) ;
        hold_var  = sdr (h , index) ;
        hold_var = hold_var (:) ;
        [final_SDR(h , i)  , index] = max(hold_var) ;
        fprintf (FileID , 'i = %d\n' , i) ;
        D(h , i) = 10 ^(- final_SDR(h , i) / 10) ;
        fprintf (FileID , '\nD = %f' , D(h , i)) ;
        fprintf (FileID , '\n%f\n' , final_SDR(h , i)) ;
        fprintf (FileID , '\nindex = %d\n' , index) ;
    end
    fprintf (FileID , 'numLevel = %d\n' , QuantizationLevels(h)) ;
end
Data = ['NACOSQ_delta_' num2str(delta)] ;
save(Data , 'final_SDR' , 'sdr' , 'epsilon')
