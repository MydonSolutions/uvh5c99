project('uvh5c99', 'c',
  version: '0.8.1',
  default_options: [
    'c_std=c99',
    'buildtype=release',
    'b_ndebug=if-release'
  ]
)

cc = meson.get_compiler('c')
m_dep = cc.find_library('m', required : true)

src_lst = []
inc_lst = [
	include_directories('include'),
]
dep_lst = [
  m_dep,
  subproject('h5dsc99').get_variable('lib_h5dsc99_dep'),
  subproject('radiointerferometryc99').get_variable('lib_radiointerferometry_dep'),
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
