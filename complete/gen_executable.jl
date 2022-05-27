using PackageCompiler
create_sysimage(["GLMakie","Spinnaker","Images","Gtk"], sysimage_path="sys_plots.so", precompile_execution_file="cameracheck.jl")