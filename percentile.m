
function val = percentile(vec,p)

    idx = floor((1-p)*length(vec));
    vec_sort = sort(vec);
    val = vec_sort(idx);
    
end