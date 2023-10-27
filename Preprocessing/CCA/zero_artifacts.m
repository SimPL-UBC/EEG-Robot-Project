function [result] =  zero_artifacts(original, to_delete)
% ZERO_ARTIFACTS(original, to_delete)   Sets the rows of input matrix
% corresponding to desired noisy channels to remove to zero
%       original: BSS-decomposed array
%       to_delete: channels/components to remove (array)

[numrows, numcols] = size(original);
for i = 1:length(to_delete)
   original(to_delete(i), :) = zeros(1, numcols);
end
result = original;

end