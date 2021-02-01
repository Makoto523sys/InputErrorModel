using LinearAlgebra
using Printf;
using Plots;
#using StatsBase;

struct point
    x::Real
    y::Real
end

function gen_cmatrix(data, mx, my)
    M = Matrix{Float64}(undef, 2, 2)
    n = length(data);
    M[1,1] = sum([(data[i].x - mx)^2 for i = 1:n]) / n;
    M[1,2] = M[2, 1] = sum([(data[i].x - mx) * (data[i].y - my) for i = 1:n]) / n;
    M[2,2] = sum([(data[i].y - my)^2 for i = 1:n]) / n;
    return M;
end

f(A, B, C, x) = -A/B*x - C/B;

function main()
    print("データ数->")
    n = parse(Int64, readline());
    data = Vector{point}(undef, n);
    for i = 1:n
        print("x座標:");
        x = parse(Int64, readline());
        print("y座標:");
        y = parse(Int64, readline());
        data[i] = point(x, y);
    end
    #重心を求める
    mx = sum([v.x for v in data]) / n;
    my = sum([v.y for v in data]) / n;
    M = gen_cmatrix(data, mx, my);
    la, v = eigen(M);
    A ,B = v[:, argmin(la)][1], v[:, argmin(la)][2]
    C = -mx*A - my*B;
    @printf("(A, B, C) = (%.4f, %.4f, %.4f)\n", A, B, C);
    plt = plot([data[i].x for i = 1:n], [data[i].y for i = 1:n], seriestype=:scatter, legend=false);
    plot!(1:200, f.(A, B, C, 1:200), legend=false)
    savefig(plt, "test.png")
end

main()
