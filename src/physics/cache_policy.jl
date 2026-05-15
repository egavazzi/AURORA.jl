default_cache_root() = pkgdir(AURORA, "internal_data")

"""
    CachePolicy(; force_recompute=false, save_cache=true,
                cache_root=default_cache_root())

Control how internal scattering and cascading caches are reused.

- `force_recompute`: ignore compatible cache files already present on disk
- `save_cache`: skip writing newly computed caches when set to `false`
- `cache_root`: parent directory that contains the `e_cascading/` and `e_scattering/`
  cache subdirectories. By default this is `pkgdir(AURORA, "internal_data")`.
"""
struct CachePolicy
    force_recompute::Bool
    save_cache::Bool
    cache_root::String
end

function CachePolicy(;
                     force_recompute::Bool = false,
                     save_cache::Bool = true,
                     cache_root::String = default_cache_root())
    return CachePolicy(force_recompute, save_cache, cache_root)
end

cache_version_string() = string(pkgversion(AURORA))
