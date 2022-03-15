project('uvh5c99', 'c',
  version: '0.4',
)

cc = meson.get_compiler('c')
m_dep = cc.find_library('m', required : true)

radiointerferometryc99_inc = include_directories('/home/sonata/dev/radiointerferometryc99/include/')
radiointerferometryc99_dep = cc.find_library(
  'radiointerferometryc99',
  dirs : ['/home/sonata/dev/radiointerferometryc99/builddir/'],
)

src_lst = []
inc_lst = [
	include_directories('include'),
  radiointerferometryc99_inc
]
dep_lst = [
  m_dep,
  dependency('hdf5'),
  dependency('erfa'),
  radiointerferometryc99_dep,
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
