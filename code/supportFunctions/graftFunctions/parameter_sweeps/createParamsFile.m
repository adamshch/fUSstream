clear 

ranges{1} = [200 500 1000 2000]; % hiddenSize
ranges{2} = [3 5 10 15 20];      % nEigs
ranges{3} = [1e-2 1e-3 1e-4];    % lambda1
ranges{4} = [1e-2 1e-3 1e-4];    % lambda2    
ranges{5} = [1e-2 1e-3 1e-4];    % lambda3


paramTable = fAllComb(ranges{:});
tableSize = size(paramTable);


f = fopen('jobs.txt', 'w');
rowString = '';
rowString = strcat(rowString, '%d \t');
rowString = strcat(rowString, '%d \t');
for j = 3:tableSize(2)
    rowString = strcat(rowString, '%6.2e \t');
end
rowString = strcat(rowString, '\r\n');

for j = 1:size( paramTable,1)
    fprintf(f, rowString, paramTable(j,:));
end
fclose(f);
