clear; 
% Sample data
s = input('Enter any even number of lattice points\n'); 
x = [s-1:-2:1,2:2:s]; 
xx=-(s/2)+1:1:s/2; 
n = numel(xx); % Number of elements in x
y = zeros(1, n); % Create a vector of zeros for y
%labels = cell(1, n);

for i = 1:n
    id = find(x == i); % Find index of element i in x
    if ~isempty(id)
        latpt(i) = x(id);
        if id+1 <= length(x)
            right(i) = x(id+1);
        else
            right(i) = i+2;
        end
        if id-1 >= 1
            left(i) = x(id-1);
        else
            left(i) = i+2;
        end
    end
end
latpt;
left;
right;   

table(latpt', left', right', 'VariableNames', {'lattice points', 'left', 'right'})



% Plotting the graph with only dots
plot(xx, y, '.','MarkerSize', 10);

for i = 1:n
    label= sprintf('%d', x(i)); 
    labels{i} = label; 
    text(xx(i), y(i), labels{i}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');  
end

% Adding labels and title
xlabel('X-axis');
ylabel('Y-axis');
title('Graph with Labeled Dots');
grid on