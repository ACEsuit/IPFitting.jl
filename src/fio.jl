module FIO

using HDF5, JSON
import FileIO

struct FHDF5 end
struct FJSON end
struct FJLD2 end

# WARNING: this seems dumb. maybe we can use empty array instead?
const NOTHING_H5 = UInt8(0)


function _format(fname)
   if length(fname) >= 3 && fname[end-2:end] == ".h5"
      return FHDF5()
   end
   if length(fname) >= 5 && fname[end-4:end] == ".json"
      return FJSON()
   end
   if length(fname) >= 5 && fname[end-4:end] == ".jld2"
      return FJLD2()
   end
   error("unsupported file format `$fname`")
end

save(fname, args...) = save(fname, Dict(args...))
save(fname, D::Dict) = _save(_format(fname), fname, D)
load(fname, args...; kwargs...) = _load(_format(fname), fname, args...; kwargs...)

function _load(format, fname, args...; kwargs...)
   D = _load(format, fname; kwargs...)
   if length(args) == 1
      return D[args[1]]
   end
   return tuple([D[key] for key in args]...)
end

# ------------- JSON LOAD AND SAVE -------------

function _load(::FJSON, fname)
   return JSON.parsefile(fname)
end

function _save(::FJSON, fname, D::Dict)
   f = open(fname, "w")
   JSON.print(f, D)
   close(f)
   return nothing
end

# ------------- JLD2 LOAD AND SAVE -------------

function _load(::FJLD2, fname, args...)
   @error "Usage of .jld2 is no longer permitted"
end

function _save(::FJLD2, fname, args...)
   @error "Usage of .jld2 is no longer permitted"
end

# ----------- HDF5 LOAD AND SAVE

function _cleanup!(X::Vector)
   for (i, x) in enumerate(X)
      X[i] = _cleanup!(x)
   end
   return X
end

function _cleanup!(D::Dict)
   for (key, val) in D
      D[key] = _cleanup!(val)
   end
   return D
end

_cleanup!(val::Any) = val
_cleanup!(val::UInt8) = (val === NOTHING_H5 ? nothing : val)



_readall(fid, mmap) =
      Dict( n => _readfrom(fid, mmap, n) for n in names(fid) )


function _load(::FHDF5, fname, args...; mmap=false)
   fid = h5open(fname, "r")
   ret = try
      if length(args) == 0
         ret = _readall(fid, mmap)
      else
         ret = [_readfrom(fid, mmap, a) for a in args]
      end
      close(fid)
      ret
   catch e
      close(fid)
      rethrow(e)
   end
   _cleanup!(ret)
   if ret isa Dict
      return ret
   elseif length(ret) == 1
      return ret[1]
   else
      return tuple(ret...)
   end
end

_readfrom(g, mmap, a::Tuple) = _readfrom(g, mmap, a...)
_readfrom(g, mmap, a::String, I::Tuple) = g[a][I...]
_readfrom(g, mmap, a::String) = mmap ? _readmmap(g[a]) : read(g, a)

_readmmap(g::HDF5Group) =
   Dict( n => _readmmap(g[n]) for n in names(g) )
_readmmap(g::HDF5Dataset) = readmmap(g)


function _save(::FHDF5, fname, D::Dict)
   h5open(fname, "w") do fid
      _hdf5_write_group!(fid, D)
   end
end

function _hdf5_write_group!(g, D)
   for (key, val) in D
      @assert key isa String
      if val isa Dict
         sub_g = g_create(g, key)
         _hdf5_write_group!(sub_g, val)
      elseif val == nothing
         g[key] = NOTHING_H5
      else
         g[key] = val
      end
   end
end


end
