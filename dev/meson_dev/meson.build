project('hello', 'fortran')

#
# Compiler flags
#
dialect = ''
bounds =  ''
if meson.get_compiler('fortran').get_id() == 'gcc'
    dialect = [ '-ffree-form','-std=f2008','-fimplicit-none']
    bounds = '-fbounds-check'
endif
if meson.get_compiler('fortran').get_id() == 'intel'
    dialect = ['-stand f08','-free','-implicitnone']
    bounds = '-check bounds'
endif
if meson.get_compiler('fortran').get_id() == 'pgi'
    dialect = [ '-Mfreeform','-Mdclchk',
                '-Mstandard','-Mallocatable=03']
    bounds = '-C'
endif

add_global_arguments(dialect, language : 'fortran')
add_global_arguments(bounds, language : 'fortran')

#
# Done.
#
sources = ['src/hello.f90']
executable('prog', sources)
