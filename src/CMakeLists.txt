add_executable(
    riemann.app 
    bc_mod.f90
    data_mod.f90
    roe_mod.f90
    fluxr_mod.f90
    init_mod.f90
    io_mod.f90
    misc_mod.f90
    params_mod.f90
    riemann.f90
)

set(project_mod_dir ${PROJECT_SOURCE_DIR}/include) 
set_target_properties(
    riemann.app   
    PROPERTIES Fortran_MODULE_DIRECTORY 
    ${project_mod_dir}  
)
