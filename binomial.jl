
function bottomup(n, k)
    if k < 0 || n < k
        return 0
    elseif k == 0 || k == n
        return 1
    end

    if 2k > n
        k = n - k
    end
    table = [1, 1]

    # Prev row begin
    prb = 2

    for i in 2:n
        push!(table, 1)
        crb = length(table)
        for j in 1:fld(i - 1, 2)
            push!(table, table[prb+j-1] + table[prb+j])
        end
        if iseven(i)
            push!(table, 2 * table[crb-1])
        end
        prb = crb
    end

    table[prb+k]
end
