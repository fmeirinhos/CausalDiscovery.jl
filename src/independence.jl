abstract type IndependenceTest end

"""
    GaussianTest(α)

Gaussian conditional independence test

# Parameters
    - `α::Real`: Significance level of H0
"""
Base.@kwdef struct GaussianTest <: IndependenceTest
    α
end

"""
    independence_test(x::AbstractMatrix, g::GaussianTest)

Conditional independence test assuming all variables are multivariate Gaussian.

Returns an unnamed `function(ix, iy, iz)` which calculates X ⊥ Y | Z, given their
respective indices, `ix`, `iy` and `iz`, along the columns of `x`.

# Notes
    https://en.wikipedia.org/wiki/Partial_correlation
"""
function independence_test(x::AbstractMatrix, g::GaussianTest)
    # Correlation matrix
    C = cor(x)

    # Number of samples
    N = size(x, 1)

    function(ix, iy, iz)
        # ρXY,Z
        ρ = let
            P = pinv(C[[ix, iy, iz...], [ix, iy, iz...]], sqrt(eps(one(eltype(C)))))
            -P[1, 2] / sqrt(P[1, 1] * P[2, 2])
        end

        # The H0 null hypothesis ρXY,Z = 0 can be accepted if
        sqrt(N - length(iz) - 3) * abs(atanh(clamp(ρ, -1, 1))) < quantile(Normal(), 1 - g.α / 2)
    end
end
