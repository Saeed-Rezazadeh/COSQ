function [indexed_T , width] = initialization (sigma , u , alpha , resolution)

MAX = 4 ; 
MIN = -4 ;
width = (MAX - MIN) / resolution ;
indexed_T = zeros (resolution + 1, 2) ;
indexed_T(: , 1) = MIN : width : MAX ;
end
