data = importdata('../hw08/code/f21.dat');
xdata = data.data(:,1);
ydata = data.data(:,2);

x = 475:775;
n = size(xdata, 1);
L = ones(n, size(x, 2));

for i = 1:n
    for j = 1:n
        if (i ~= j)
            L(i, :) = L(i, :) .* (x - xdata(j)) / (xdata(i) - xdata(j));
        end
    end
end
y = 0;
for i = 1:n
    y = y + ydata(i) * L(i, :);
end

figure
plot(x, y)