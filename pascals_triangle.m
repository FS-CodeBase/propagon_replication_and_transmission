function pt =  pascals_triangle(n)
% % FUNCTION: pascals_triangle
% % Author: Fabian Santiago 
% % Description: Creates a square matrix of the values for Pascal and
% %              returns an anonymous function for the values in each row
% %              of the triangle up to the nth row (n). 

% Pre-allocate a square matrix for the values
tmp_pt = zeros(n+1);

% Compute the values of Pascal's triangle
for i = 0:n
    for j = 0:n
        if j <= i
        tmp_pt(j+1,i+1) = nchoosek(i,j);
        end
    end
end

% Function that returns the non-zero values in a row of pascals triangle.
pt =@(i) tmp_pt(1:(i+1),i+1); 