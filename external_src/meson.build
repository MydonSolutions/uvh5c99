radiointerferometryc99_dep = cc.find_library(
  'radiointerferometryc99',
  dirs : [meson.current_source_dir()+'/radiointerferometryc99/builddir/'],
)

dep_lst += [
  radiointerferometryc99_dep
]

inc_lst += [
  include_directories('tomlc99'),
  include_directories('radiointerferometryc99/include'),
]

src_lst += files([
  'tomlc99/toml.c',
])