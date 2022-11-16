# UVH(df)5 in C

Based on V1.1 [this memo](https://github.com/RadioAstronomySoftwareGroup/pyuvdata/blob/69b9fa15669a7c449f65c50cdc4f5f01901085c0/docs/references/uvh5_memo.pdf) (see Section 5.4, Version History: Version 1.1).

## Compilation and installation

`$ git submodule update --init`

`# apt install libhdf5-dev cmake liberfa-dev`

`$ pip install meson ninja`

`$ meson build && cd build && ninja test && ninja install`
