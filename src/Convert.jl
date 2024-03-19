module Convert
export cart2coes, cart2coes_validated, coes2cart, coes2cart_validated

using LinearAlgebra

# TODO(melmer): Add release branch

# TODO(melmer): Use StaticArray

"""
    cart2coes(cartesian_state_vector, μ, ε=1e-12)

Based on Algorithm 9 of "Fundamentals of Astrodynamics and Applications".

For each of Ω and ω, if it's undefined, set it to zero.

Units in: [km, km, km, km/s, km/s, km/s]

Units out: [km, n.d., rad, rad, rad, rad]
"""
function cart2coes(cart::Vector{<:AbstractFloat}, μ::AbstractFloat, ε::AbstractFloat=1e-12)
    r_vec = cart[1:3]
    v_vec = cart[4:6]

    r = norm(r_vec)
    v = norm(v_vec)

    h_vec = r_vec × v_vec
    h = norm(h_vec)  # QUESTION(melmer): Warn if h near zero?

    k_hat = [0, 0, 1]
    n_vec = k_hat × h_vec
    n = norm(n_vec)
    equatorial = abs(n) < ε

    e_vec = (r_vec * (v^2 - μ/r) - v_vec * (r_vec ⋅ v_vec)) / μ
    e = norm(e_vec)
    parabolic = abs(e - 1) < ε
    circular = abs(e) < ε

    a = undef
    if parabolic
        a = Inf
    else
        a = 1 / (2/r - v^2/μ)
    end

    i = acos(h_vec[3]/h)

    Ω = undef
    if equatorial
        Ω = 0.0
    else
        Ω = acos(n_vec[1] / n)
        if n_vec[2] < 0.0
            Ω = 2π - Ω
        end
    end

    ω = undef
    if circular
        ω = 0.0
    elseif equatorial
        ω = acos(e_vec[1] / e)
        if e_vec[2] < 0.0
            ω = 2π - ω
        end
    else
        ω = acos(n_vec ⋅ e_vec / (n * e))
        if e_vec[3] < 0.0
            ω = 2π - ω
        end
    end

    ν = undef
    if circular && equatorial
        ν = acos(r_vec[1] / r)
        if r_vec[2] < 0.0
            ν = 2π - ν
        end
    elseif circular
        ν = acos(n_vec ⋅ r_vec / (n * r))
        if r_vec[3] < 0.0
            ν = 2π - ν
        end
    else
        ν = acos(e_vec ⋅ r_vec / (e * r))
        if r_vec ⋅ v_vec < 0.0
            ν = 2π - ν
        end
    end

    return [a, e, i, ω, Ω, ν]
end

# Validate input before calling cart2coes.
function cart2coes_validated(cart::Vector{<:AbstractFloat}, μ::AbstractFloat, ε::AbstractFloat=1e-12)
    if !all(isfinite, cart) || !isfinite(μ)
        @warn "Invalid input passed to cart2coes_validated."
        # TODO(melmer): Throw exception instead?
        return [NaN, NaN, NaN, NaN, NaN, NaN]
    end
    return cart2coes(cart, μ, ε)
end

"""
    coes2cart(coes, μ, ε=1e-12)

Based on Algorithm [TODO(melmer)] of "Fundamentals of Astrodynamics and
Applications".

Units in: [km, n.d., rad, rad, rad, rad]

Units out: [km, km, km, km/s, km/s, km/s]
"""
function coes2cart(coes::Vector{<:AbstractFloat}, μ::AbstractFloat, ε::AbstractFloat=1e-12)
    a, e, i, ω, Ω, ν = coes

    circular = abs(e) < ε
    if circular
        ω = 0.0
    end

    equatorial = i < ε || π - i < ε
    if equatorial
        Ω = 0.0
    end

    p = a * (1 - e^2)

    r_vec_pqw = [
        p * cos(ν) / (1 + e * cos(ν))
        p * sin(ν) / (1 + e * cos(ν))
        0.0
    ]
    v_vec_pqw = [
        -sqrt(μ / p) * sin(ν)
         sqrt(μ / p) * (e + cos(ν))
        0.0
    ]

    _11 = cos(Ω) * cos(ω) - sin(Ω) * sin(ω) * cos(i)
    _12 = -cos(Ω) * sin(ω) - sin(Ω) * cos(ω) * cos(i)
    _13 = sin(Ω) * sin(i)
    _21 = sin(Ω) * cos(ω) + cos(Ω) * sin(ω) * cos(i)
    _22 = -sin(Ω) * sin(ω) + cos(Ω) * cos(ω) * cos(i)
    _23 = -cos(Ω) * sin(i)
    _31 = sin(ω) * sin(i)
    _32 = cos(ω) * sin(i)
    _33 = cos(i)
    ijk_pqw = [
        _11 _12 _13
        _21 _22 _23
        _31 _32 _33
    ]

    r_vec = ijk_pqw * r_vec_pqw
    v_vec = ijk_pqw * v_vec_pqw

    return [r_vec; v_vec]
end

# Validate input before calling coes2cart
function coes2cart_validated(coes::Vector{<:AbstractFloat}, μ::AbstractFloat, ε::AbstractFloat=1e-12)
    if !all(isfinite, coes) || !isfinite(μ)
        @warn "Invalid input passed to coes2cart_validated"
        # TODO(melmer): Throw exception instead?
        return [NaN, NaN, NaN, NaN, NaN, NaN]
    end
    a, e, i, ω, Ω, ν = coes
    parabolic = abs(e - 1) < ε  # Would check if a isn't finite but already did
    if parabolic
        throw(ArgumentError(e,
            "coes2cart cannot convert parabolic trajectories. Consider using " *
            "a parameterization that isn't singular at zero eccentricity."
        ))
    end
    if abs(a) < ε
        @warn "Invalid a: $a"
    end
    if e < 0
        @warn "Invalid e: $e"
    end
    if i < 0 || i > π
        @warn "Invalid i: $i"
    end
    if ω < 0.0 || ω > 2π
        @warn "Invalid ω: $ω"
    end
    if Ω < 0.0 || Ω > 2π
        @warn "Invalid Ω: $Ω"
    end
    if ν < 0.0 || ν > 2π
        @warn "Invalid ν: $ν"
    end
    return coes2cart(coes, μ, ε)
end

end
