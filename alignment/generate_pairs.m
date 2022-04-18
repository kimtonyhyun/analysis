function pairs = generate_pairs(num_sessions)

pairs = zeros(num_sessions*(num_sessions-1)/2, 2);
idx = 0;
for i = 1:num_sessions-1
    for j = i+1:num_sessions
        idx = idx + 1;
        pairs(idx,:) = [i j];
    end
end