function __init__()
    @static if Sys.iswindows()
      paths = [joinpath(ENV["ProgramFiles"], "Point Grey Research", "Spinnaker", "bin", "vs2015")]
      libspinnaker = "SpinnakerC_v140.dll"
      libspinvideo = ""
    elseif Sys.islinux()
      paths = ["/usr/lib" "/opt/spinnaker/lib"]
      libspinnaker = "libSpinnaker_C.so"
      libspinvideo = "libSpinVideo_C.so"
    elseif Sys.isapple()
      paths = ["/usr/local/lib"]
      libspinnaker = "libSpinnaker_C.dylib"
      libspinvideo = "libSpinVideo_C.dylib"
    else
      @error "Spinnaker SDK only supported on Linux, Windows and MacOS platforms"
      return
    end
    libSpinnaker_C_path = ""
    libSpinVideo_C_path = ""
    for path in paths
      libSpinnaker_C_path = joinpath(path, libspinnaker)
      libSpinVideo_C_path = joinpath(path, libspinvideo)
      if isfile(libSpinnaker_C_path) && isfile(libSpinVideo_C_path)
        libSpinnaker_C[] = libSpinnaker_C_path
        libSpinVideo_C[] = libSpinVideo_C_path
      end
    end

    if libSpinnaker_C[] == "" || libSpinVideo_C[] == ""
        @error "Spinnaker SDK cannot be found. This package can be loaded, but will not be functional."
        return
      end
      try
        libSpinnaker_C_handle = dlopen(libSpinnaker_C[])
        !Sys.iswindows() && (libSpinVideo_C_handle = dlopen(libSpinVideo_C[]))
      catch ex
        bt = catch_backtrace()
        @error "Spinnaker SDK cannot be dlopen-ed"
        showerror(stderr, ex, bt)
      end
      try
        global spinsys = System()
      catch ex
        bt = catch_backtrace()
        @error "Spinnaker SDK loaded but Spinnaker.jl failed to initialize"
        showerror(stderr, ex, bt)
      end
    end
end