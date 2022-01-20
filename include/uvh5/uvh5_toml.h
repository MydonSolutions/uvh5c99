#ifndef UVH5_TOML_H
#define UVH5_TOML_H

#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include "toml.h"
#include "uvh5.h"
#include "uvh5_printer.h"

void UVH5toml_parse_telescope_info(char* file_path, UVH5_header_t* UVH5header);
void UVH5toml_parse_observation_info(char* file_path, UVH5_header_t* UVH5header);

#endif