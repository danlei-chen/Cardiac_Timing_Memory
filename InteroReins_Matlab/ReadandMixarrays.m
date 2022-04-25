function out = ReadandMixarrays(subStim, numImagesSelected_subarrays, num_sameConsecutiveCategory)

%This function takes in multiple horizontal arrays (subarrays), selects a number of imagesfrom each array (numSelected_insubarrays), combines them,
%and then randomizes the order of the elements in the new array such that no more than a number images of same category are shown in a row (num_sameCategory).

%Concatonates the two horizontal arrays into one long horizontal array
all_arrays = subStim;

%Measures the size of the new array
num_elements=length(C);

%Randomly generates new indecies to generate the randomize the elements
indexvector=randperm(num_elements);

%Initialize the output vector
out=zeros(1,num_elements);

%Writes the ourput vector
for i=1:num_elements
    out(i)=C(indexvector(i));    
end
