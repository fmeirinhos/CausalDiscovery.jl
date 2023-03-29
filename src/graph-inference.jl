"""
    infer_graph(x::AbstractMatrix, test::IndependenceTest)

Infers the causal graph underlying the set of variables along the columns of `x`.
"""
function infer_graph(x::AbstractMatrix, test::IndependenceTest)
    # Find fixed point
    fix(f, x) =
        let
            x′ = f(copy(x))
            x′ == x ? x : fix(f, x′)
        end

    # Infer skeleton
    g, S = skeleton(complete_graph(size(x, 2)), independence_test(x, test))

    # Orient skeleton
    dg = fix(rule3! ∘ rule2! ∘ rule1!, orient_colliders(g, S))

    dg
end

# Infers the skeleton graph
function skeleton(g::SimpleGraph{T}, independence) where {T}
    S = Dict{Edge{T},Vector{T}}()

    for ℓ ∈ 1:nv(g)
        for e ∈ collect(edges(g))
            U = setdiff(union(neighbors(g, src(e)), neighbors(g, dst(e))), src(e), dst(e))
            for s ∈ combinations(U, ℓ - 1)
                if independence(src(e), dst(e), s) # Remove conditionally-independent edge
                    rem_edge!(g, e)
                    S[e] = s
                    S[reverse(e)] = s
                    @info "Removing $(e) given $(s)"
                    break
                end
            end
        end
    end
    g, S
end

# For non-adjacent a and b, a -- c -- b, if c is not in their separating set, c must be a collider, as conditioning on it introduces dependence
function orient_colliders(g::SimpleGraph, S)
    dg = DiGraph(g)

    for c ∈ vertices(g)
        for (a, b) ∈ filter(x -> !has_edge(g, x...), combinations(neighbors(g, c), 2))
            if c ∉ S[Edge(a, b)]
                rem_edge!(dg, Edge(c, a))
                rem_edge!(dg, Edge(c, b))
                @info "Removing Edge $a <= $c => $b (Rule 0)"
            end
        end
    end
    dg
end

# Orients b -- c into b -> c whenever a -> b such that a and c are nonadjacent
function rule1!(dg::SimpleDiGraph)
    for b ∈ vertices(dg)
        for a ∈ filter(a -> !has_edge(dg, b, a), inneighbors(dg, b)) |> collect
            for c ∈ filter(c -> are_nonadjecent(dg, a, c), outneighbors(dg, b)) |> collect
                if rem_edge!(dg, c, b)
                    @info "Removing $(Edge(c, b)) (Rule I)"
                end
            end
        end
    end
    dg
end

# Orients a -- b into a -> b whenever there is chain a -> c -> b
function rule2!(dg::SimpleDiGraph)
    for c ∈ vertices(dg)
        for a ∈ filter(a -> !has_edge(dg, c, a), inneighbors(dg, c)) |> collect
            for b ∈ filter(b -> !has_edge(dg, b, c), outneighbors(dg, c)) |> collect
                if has_edge(dg, a, b) && rem_edge!(dg, b, a)
                    @info "Removing $(Edge(b, a)) (Rule II)"
                end
            end
        end
    end
    dg
end

# Orients a -- b into a -> b whenever there are two chains a -- c -> b and a -- d -> b such that c and d are nonadjacent
function rule3!(dg::SimpleDiGraph)
    for b ∈ vertices(dg)
        for (c, d) ∈ filter(x -> are_nonadjecent(dg, x...), combinations(directed_inneighbors(dg, b), 2)) |> collect
            for a ∈ filter(a -> is_bidirected(dg, a, d) && is_bidirected(dg, a, c), setdiff(common_neighbors(dg, c, d), b)) |> collect
                if has_edge(dg, a, b) && rem_edge!(dg, b, a)
                    @info "Removing $(Edge(b, a)) (Rule III)"
                end
            end
        end
    end
    dg
end

are_nonadjecent(dg, a, b) = !has_edge(dg, a, b) && !has_edge(dg, b, a)

is_bidirected(dg, a, b) = has_edge(dg, a, b) && has_edge(dg, b, a)

directed_inneighbors(dg, a) = setdiff(inneighbors(dg, a), outneighbors(dg, a))
