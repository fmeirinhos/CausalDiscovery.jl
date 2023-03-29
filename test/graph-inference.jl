using CausalDiscovery.Graphs

make_graph(graph) = foldl((g, e) -> (add_edge!(g, e...); g), graph; init = DiGraph(maximum(maximum, graph)))

@test begin
    (; x, graph) = let
        N = 100_000
        x1 = randn(N)
        x2 = x1 + randn(N) * 0.25
        x3 = x1 + randn(N) * 0.25
        x4 = x2 + x3 + randn(N) * 0.25
        x5 = x4 + randn(N) * 0.25

        x = [x1 x2 x3 x4 x5]
        graph = [1 => 2, 1 => 3, 2 => 4, 3 => 4, 4 => 5]

        (; x, graph)
    end

    issubset(make_graph(graph), infer_graph(x, GaussianTest(α = 0.01)))
end

@test begin
    (; x, graph) = let
        N = 100_000
        x1 = randn(N)
        x2 = x1 + randn(N) * 0.25
        x3 = x2 + randn(N) * 0.25
        x4 = x1 + x3 + randn(N) * 0.25

        x = [x1 x2 x3 x4]
        graph = [1 => 2, 2 => 3, 1 => 4, 3 => 4]

        (; x, graph)
    end

    issubset(make_graph(graph), infer_graph(x, GaussianTest(α = 0.01)))
end

@test begin
    (; x, graph) = let
        N = 100_000
        x1 = randn(N)
        x2 = x1 + randn(N) * 0.25
        x3 = (x1 + x2) / 2 + randn(N) * 0.25
        x4 = x2 + randn(N) * 0.25
        x5 = x4 + randn(N) * 0.25
        x6 = (x1 + x3 + x5) / 3 + randn(N) * 0.25
        x7 = x6 + randn(N) * 0.25

        x = [x1 x2 x3 x4 x5 x6 x7]
        graph = [1 => 2, 1 => 3, 1 => 6, 2 => 3, 2 => 4, 3 => 6, 4 => 5, 5 => 6, 6 => 7]

        (; x, graph)
    end

    issubset(make_graph(graph), infer_graph(x, GaussianTest(α = 0.01)))
end

@test begin
    (; x, graph) = let
        N = 100_000
        x1 = randn(N)
        x2 = x1 + randn(N) * 0.25
        x3 = (x1 + x2) / 2 + randn(N) * 0.25

        x5 = randn(N)
        x6 = x5 + randn(N) * 0.25

        x4 = (x3 + x6) / 2 + randn(N) * 0.25
        x7 = (x3 + x4 + x6) / 3 + randn(N) * 0.25
        x8 = (x6 + x7) / 2 + randn(N) * 0.25

        x = [x1 x2 x3 x4 x5 x6 x7 x8]
        graph = [1 => 2, 1 => 3, 2 => 3, 3 => 4, 3 => 7, 5 => 6, 6 => 4, 6 => 8, 6 => 7, 4 => 7, 7 => 8]
        
        (; x, graph)
    end

    issubset(make_graph(graph), infer_graph(x, GaussianTest(α = 0.01)))
end
