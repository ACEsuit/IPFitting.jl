module FIO

using HDF5

struct FHDF5 end

const NOTHING = UInt8(0)

function _format(fname)
   if length(fname) >= 3
      if fname[end-2:end] == ".h5"
         return FHDF5()
      end
   end
   error("unsupported file format `$fname`")
end

save(fname, args...) = save(fname, Dict(args...))
save(fname, D::Dict) = _save(_format(fname), fname, D)
load(fname, args...) = _load(_format(fname), fname, args...)


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
_cleanup!(val::UInt8) = (val === NOTHING ? nothing : val)

function _load(::FHDF5, fname, args...)
   fid = h5open(fname, "r")
   ret = try
      if length(args) == 0
         ret = read(fid)
      else
         ret = [_readfrom(fid, a) for a in args]
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

_readfrom(g, a::Tuple) = _readfrom(g, a...)
_readfrom(g, a::String, I::Tuple) = g[a][I...]
_readfrom(g, a::String) = read(g, a)



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
         g[key] = NOTHING
      else
         g[key] = val
      end
   end
end


end
