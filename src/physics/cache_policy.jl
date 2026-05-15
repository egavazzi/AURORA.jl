struct CachePolicy
    force_recompute::Bool
    save_cache::Bool
    cache_root::Union{Nothing, String}
end

function CachePolicy(;
                     force_recompute::Bool = false,
                     save_cache::Bool = true,
                     cache_root::Union{Nothing, String} = nothing)
    return CachePolicy(force_recompute, save_cache, cache_root)
end

cache_version_string() = string(pkgversion(AURORA))
