using SparseArrays

function jacpattern(l)
    np = length(l); njace = 5*l[end]-7*np + np*(np-1);
    I, J = [ones(njace) for ii = 1:2]

    pos = 1; lprev = 1;

    for kk = 1:np
        for ii = lprev:l[kk]
            if ii <= lprev + 1
                range = (lprev-ii):2
            elseif ii == lprev + 2
                range = -1:2
            elseif ii >= l[kk] - 1
                range = -2:(l[kk]-ii)
            else
                range = -2:2
            end
                
            for jj = range
                I[pos] = ii; J[pos] = ii + jj;
                pos += 1;
            end
        end
        lprev = l[kk] + 1
    end

    for ii = 2:np
        for jj = 1:ii-1
            I[pos] = l[ii]; J[pos] = l[jj];
            I[pos+1] = l[jj]; J[pos+1] = l[ii];
            pos += 2
        end
    end

    nf = pos - njace - 1

    if nf != 0
        print("ERRO"); print("\n");
        print(nf); print("\n");
    end            

    jc = Float64.(sparse(I, J, ones(njace)));
    return jc
end