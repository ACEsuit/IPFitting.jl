module FIO

using HDF5

struct FHDF5 end

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

function _load(::FHDF5, fname, args...)
   fid = h5open(fname, "r")
   if length(args) == 0
      D = read(fid)
      close(fid)
      return D
   end
   ret = [_readfrom(fid, a) for a in args]
   close(fid)
   return length(ret) == 1 ? ret[1] : tuple(ret...)
end

_readfrom(fid, a::String) = read(fid, a)


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
      else
         g[key] = val
      end
   end
end


end
