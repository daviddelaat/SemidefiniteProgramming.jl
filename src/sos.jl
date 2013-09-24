function addsoscon!(sdp::SparseSDP,
                    id,
                    fs,
                    gs,
                    nvars=1,
                    constant=0.0) 
    
    setrhs!(sdp, (id, tuple(ones(Int, nvars)...)), -constant)

    for (bi, i, j, pol) in fs
        @forcartesian monomial size(pol) begin
            v = pol[monomial...]
            if abs(v) != 0.0
                setcon!(sdp, (id, tuple(monomial...)), bi, i, j, pol[monomial...])
            end
        end
    end
    
    for i = 1:length(gs)
        @forcartesian mg size(gs[i][1]) begin
            @forcartesian mx gs[i][2] begin
                @forcartesian my gs[i][2] begin
                    v = -gs[i][1][mg...]
                    monomial = mg-1 + mx-1 + my-1 + 1
                    if abs(v) != 0.0
                       setcon!(sdp, (id, tuple(monomial...)), (id, i), tuple(mx...), tuple(my...), v)
                    end
                end
            end
        end
    end
end
