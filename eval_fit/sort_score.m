function sorted_score = sort_score(score)

[sorted_vals_n, sort_idx] = sort(score.vals_n, score.sort_order);
sorted_score = [sort_idx sorted_vals_n];