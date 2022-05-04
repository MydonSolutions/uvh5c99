project('uvh5c99', 'c',
  version: '0.5',
)

cc = meson.get_compiler('c')
m_dep = cc.find_library('m', required : true)

src_lst = []
inc_lst = [
	include_directories('include'),
]
dep_lst = [
  m_dep,
  dependency('hdf5'),
]

subdir('external_src')
subdir('include')
subdir('src')

lib_uvh5 = library(
  'uvh5',
  src_lst,
  include_directories: inc_lst,
  dependencies: dep_lst,
  install: true,
)

lib_uvh5_dep = declare_dependency(
  include_directories: inc_lst,
  dependencies: dep_lst,
  link_with: lib_uvh5,
)

subdir('tests')
