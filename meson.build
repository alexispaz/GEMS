project('GeMS', 'fortran',
     meson_version : '>=1.1',
)

add_project_arguments('-DPACKAGE_VERSION="' + meson.project_version() + '"', language: 'fortran')

if get_option('sprng')
#   sprng_proj = mod.add_project('sprng/configure', verbose: true)
#   sprng_dep = sprng_proj.dependency('libsprng.a')
   sprng_src = '../src/Random_sprng.F90'
else
#   sprng_dep = dependency('', required: false)
  sprng_src = '../src/Random.F90'
endif 

# Dependencias externas via wrap
fpt_dep = dependency('fpt', fallback: ['fpt', 'fpt_dep'])
fparser_dep = dependency('fparser', fallback: ['fparser', 'fparser_dep'])

subdir('src')
