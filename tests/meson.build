test('write', executable(
  'write', ['write.c'],
	dependencies: lib_uvh5c99_dep,
	),
	args : files(['telinfo_ata.toml', 'obsinfo.toml']),
	is_parallel: false
)
test('uvw_phased', executable(
  'uvw_phased', ['uvw_phased.c', 'uvw_phased.h'],
	dependencies: lib_uvh5c99_dep,
	),
	args : files(['telinfo_ata.toml', 'obsinfo.toml']),
	is_parallel: false
)
test('visdata_from_xgpu', executable(
  'visdata_from_xgpu', ['visdata_from_xgpu.c', 'matrix_h.h', 'visdata.h'],
	dependencies: lib_uvh5c99_dep,
	),
	args : files(['telinfo_ata.toml', 'obsinfo.toml']),
	is_parallel: false
)
