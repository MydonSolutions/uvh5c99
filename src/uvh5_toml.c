#include "uvh5/uvh5_toml.h"

// Parse a '%d+:%d+:%d+' sexgesimal string to a double
// Return zero if success, else non-zero
double sexagesimal_to_double(const char *sexagesiamal_string) {
    double f = 0.0;
    char *pend;
    char *pstart = (char *) sexagesiamal_string;
    if (*sexagesiamal_string == '-') {
        pstart += 1;
    }
    do {
        f += strtof(pstart, &pend);
        f *= 60;
        pstart = pend+1;
    } while(*(pend) == ':');
    f /= 60*60*60;
    if (*sexagesiamal_string == '-') {
        f *= -1.0;
    }
    return f;
}

// Parse a sexgesimal value (string or double) at `location`, populating the `sexagesimal` double
// Return zero if success, else non-zero
int _UVH5toml_sexagesimal_in(toml_table_t* parent, const char* location, double* sexagesimal) {
    toml_datum_t toml_datum = toml_string_in(parent, location);
    if (!toml_datum.ok) {
        toml_datum = toml_double_in(parent, location);
        if (!toml_datum.ok) {
            UVH5print_error(__FUNCTION__, "cannot read location '%s'", location);
            return 1;
        }
        *sexagesimal = toml_datum.u.d;
        UVH5print_verbose(__FUNCTION__, "%s: %f", location, *sexagesimal);
    }
    else {
        // process sexagesimal string to double
        *sexagesimal = sexagesimal_to_double(toml_datum.u.s);
        free(toml_datum.u.s);
        UVH5print_verbose(__FUNCTION__, "%s: %f", location, *sexagesimal);
    }
    return 0;
}

// Parse the string at the array's index `idx`, populating `string_out`
// Return zero if success, else non-zero
int _UVH5toml_string_at(toml_array_t* parent, const int idx, char** string_out) {
    toml_datum_t toml_datum = toml_string_at(parent, idx);
    if (!toml_datum.ok) {
        UVH5print_error(__FUNCTION__, "cannot read index %d", idx);
        return 1;
    }
    else {
        *string_out = malloc(strlen(toml_datum.u.s)+1);
        memcpy(*string_out, toml_datum.u.s, strlen(toml_datum.u.s)+1);
        free(toml_datum.u.s);
        UVH5print_verbose(__FUNCTION__, "[%d]: %s", idx, *string_out);
    }
    return 0;
}

// Parse the string at the array's index `idx`, populating `string_out`, with limited `length`
// Return zero if success, else non-zero
int _UVH5toml_nstring_at(toml_array_t* parent, const int idx, char* string_out, size_t length) {
    toml_datum_t toml_datum = toml_string_at(parent, idx);
    if (!toml_datum.ok) {
        UVH5print_error(__FUNCTION__, "cannot read idx %s", idx);
        return 1;
    }
    else {
        length = length <= strlen(toml_datum.u.s) ? length : strlen(toml_datum.u.s);
        memcpy(string_out, toml_datum.u.s, length);
        free(toml_datum.u.s);
        UVH5print_verbose(__FUNCTION__, "[%d]: %s", idx, string_out);
    }
    return 0;
}

// Parse the string at `location`, populating `string_out`
// Return zero if success, else non-zero
int _UVH5toml_string_in(toml_table_t* parent, const char* location, char** string_out) {
    toml_datum_t toml_datum = toml_string_in(parent, location);
    if (!toml_datum.ok) {
        UVH5print_error(__FUNCTION__, "cannot read location", location);
        return 1;
    }
    else {
        *string_out = malloc(strlen(toml_datum.u.s)+1);
        memcpy(*string_out, toml_datum.u.s, strlen(toml_datum.u.s)+1);
        free(toml_datum.u.s);
        UVH5print_verbose(__FUNCTION__, "%s: %s", location, *string_out);
    }
    return 0;
}

// Parse the string at `location`, populating `string_out`, with limited `length`
// Return zero if success, else non-zero
int _UVH5toml_nstring_in(toml_table_t* parent, const char* location, char* string_out, size_t length) {
    toml_datum_t toml_datum = toml_string_in(parent, location);
    if (!toml_datum.ok) {
        UVH5print_error(__FUNCTION__, "cannot read location '%s'", location);
        return 1;
    }
    else {
        memset(string_out, '\0', length);
        length = length <= strlen(toml_datum.u.s)+1 ? length : strlen(toml_datum.u.s)+1;
        memcpy(string_out, toml_datum.u.s, length);
        free(toml_datum.u.s);
        UVH5print_verbose(__FUNCTION__, "%s: %s", location, string_out);
    }
    return 0;
}

// Parse the double at array's index `idx`, populating `double_out`
// Return zero if success, else non-zero
int _UVH5toml_double_at(toml_array_t* parent, const int idx, double* double_out) {
    toml_datum_t toml_datum = toml_double_at(parent, idx);
    if (!toml_datum.ok) {
        UVH5print_error(__FUNCTION__, "cannot read index %d", idx);
        return 1;
    }
    else {
        *double_out = toml_datum.u.d;
        UVH5print_verbose(__FUNCTION__, "[%d]: %f", idx, *double_out);
    }
    return 0;
}
int _UVH5toml_float_at(toml_array_t* parent, const int idx, float* float_out) {
    double intermediate;
    if(_UVH5toml_double_at(parent, idx, &intermediate)){
        UVH5print_error(__FUNCTION__, "cannot read index %d", idx);
    }
    *float_out = (float) intermediate;
    return 0;
}

// Parse the double at `location`, populating `double_out`
// Return zero if success, else non-zero
int _UVH5toml_double_in(toml_table_t* parent, const char* location, double* double_out) {
    toml_datum_t toml_datum = toml_double_in(parent, location);
    if (!toml_datum.ok) {
        UVH5print_error(__FUNCTION__, "cannot read location '%s'", location);
        return 1;
    }
    else {
        *double_out = toml_datum.u.d;
        UVH5print_verbose(__FUNCTION__, "%s: %f", location, *double_out);
    }
    return 0;
}
int _UVH5toml_float_in(toml_table_t* parent, const char* location, float* float_out){
    double intermediate;
    if(_UVH5toml_double_in(parent, location, &intermediate)) {
        UVH5print_error(__FUNCTION__, "cannot read location '%s'", location);
    }
    *float_out = (float) intermediate;
    return 0;
}

// Parse the int at array's index `idx`, populating `int_out`
// Return zero if success, else non-zero
int _UVH5toml_int_at(toml_array_t* parent, const int idx, int* int_out) {
    toml_datum_t toml_datum = toml_int_at(parent, idx);
    if (!toml_datum.ok) {
        UVH5print_error(__FUNCTION__, "cannot read index %d", idx);
        return 1;
    }
    else {
        *int_out = toml_datum.u.i;
        UVH5print_verbose(__FUNCTION__, "[%d]: %d", idx, *int_out);
    }
    return 0;
}

// Parse the int at `location`, populating `int_out`
// Return zero if success, else non-zero
int _UVH5toml_int_in(toml_table_t* parent, const char* location, int* int_out) {
    toml_datum_t toml_datum = toml_int_in(parent, location);
    if (!toml_datum.ok) {
        UVH5print_error(__FUNCTION__, "cannot read location", location);
        return 1;
    }
    else {
        *int_out = toml_datum.u.i;
        UVH5print_verbose(__FUNCTION__, "%s: %d", location, *int_out);
    }
    return 0;
}

// Process the toml_table_t `parent`, using entries to populate arguments:
//  "diameter" -> `ant_diameter`
//  "number" -> `ant_id`
//  "name" -> `ant_name`
//  "position[0, 1, 2]" -> `ant_pos[0, 1, 2]` (x, y, z)
int _UVH5toml_antenna_table_in(
    toml_table_t* parent,
    UVH5_antinfo_t* ant
){
    if(_UVH5toml_double_in(parent, "diameter", &ant->diameter)) {
        UVH5print_warn(__FUNCTION__, "cannot read diameter");
    }
    if(_UVH5toml_string_in(parent, "name", &ant->name)) {
        UVH5print_error(__FUNCTION__, "cannot read name");
    }
    if(_UVH5toml_int_in(parent, "number", &ant->number)) {
        UVH5print_error(__FUNCTION__, "cannot read number");
    }
    toml_array_t* toml_ant_position = toml_array_in(parent, "position");
    if(toml_ant_position == NULL) {
        UVH5print_error(__FUNCTION__, "cannot read position");
    }
    
    int ant_pos_count = toml_array_nelem(toml_ant_position);
    if(ant_pos_count != 3) {
        UVH5print_error(__FUNCTION__, "unexpected position length of %d", ant_pos_count);
    }
    for (size_t i = 0; i < 3; i++)
    {
        if(_UVH5toml_double_at(toml_ant_position, i, ant->position+i)) {
            UVH5print_error(__FUNCTION__, "cannot read position float at index %d", i);
            ant->position[i] = 0.0f;
        }
    }
    return 0;
}

// Parse the telescope_info toml file, populating the relevant fields in `header`.
// The 'tests/telinfo_ata.toml' is exemplary.
//
//  "telescope_name" 			-> header->telescope_name
//  "latitude" 					-> header->latitude
//  "longitude" 				-> header->longitude
//  "altitude" 					-> header->altitude
//  "reference_antenna_name"	-> header->extra_keywords->reference_antenna_name (optional, inconsequential)
//  "antenna_diameter"			-> header->antenna_diameters[...]
//  len("[[antennas]]")			-> header->Nants_telescope
//  "[[antennas]]"[i]			-> header->antenna_numbers[i]
//									 ...-> header->antenna_names[i]
//									 ...-> header->antenna_positions[i]
//		See `_UVH5toml_antenna_table_in`
//
// All of the associated pointer-arrays in `header` are allocated before being
// populated.
//
// The antenna_positions are converted to an "xyz" frame based on
// the frame indicated by "antenna_position_frame".
void UVH5toml_parse_telescope_info(const char* file_path, UVH5_header_t* header) {
    FILE* fp;
    char errbuf[200];

    fp = fopen(file_path, "r");
    if (!fp) {
        UVH5print_error(__FUNCTION__, "cannot open '%s' (%s)", file_path, strerror(errno));
        return;
    }

    toml_table_t* conf = toml_parse_file(fp, errbuf, sizeof(errbuf));
    fclose(fp);
    if (!conf) {
        UVH5print_error(__FUNCTION__, "cannot parse %s", errbuf);
        return;
    }
    
    char* name;
    double latitude;					/* The latitude of the telescope site, in degrees. */
    double longitude;					/* The longitude of the telescope site, in degrees. */
    double altitude;					/* The altitude of the telescope site, in meters. */
    char* antenna_position_frame;       /* The antenna position frame: {"XYZ", "ENU", "ECEF"}. */
    double default_antenna_diameter;    /* The fallback antenna diameter, in meters. */
    int nantenna;                       /* The number of antenna (the length of the `antenna` field). */
    UVH5_antinfo_t* antenna;            /* The antenna information with positions. */

    _UVH5toml_string_in(conf, "telescope_name", &name);
    _UVH5toml_sexagesimal_in(conf, "latitude", &latitude);
    _UVH5toml_sexagesimal_in(conf, "longitude", &longitude);
    _UVH5toml_double_in(conf, "altitude", &altitude);

    _UVH5toml_double_in(conf, "antenna_diameter", &default_antenna_diameter);
    
    toml_array_t* toml_antennas_array = toml_array_in(conf, "antennas");
    if (!toml_antennas_array) {
        UVH5print_error(__FUNCTION__, "missing [[antennas]]");
        return;
    }
    
    nantenna = toml_array_nelem(toml_antennas_array);
    antenna = malloc(sizeof(UVH5_antinfo_t)*nantenna);

    for (size_t i = 0; i < nantenna; i++)
    {
        toml_table_t* toml_antenna_info = toml_table_at(toml_antennas_array, i);
        if(toml_antenna_info) {
            UVH5print_verbose(__FUNCTION__, "Antenna: %ld", i);
            _UVH5toml_antenna_table_in(
                toml_antenna_info,
                antenna+i
            );
        }
        else {
            UVH5print_error(__FUNCTION__, "cannot access antenna #%ld info", i);
        }
    }

    if(_UVH5toml_string_in(conf, "antenna_position_frame", &antenna_position_frame)) {
        // Not specified
        antenna_position_frame = NULL;
    }

    toml_free(conf);
    UVH5set_telescope_info(
        name,
        latitude,
        longitude,
        altitude,
        antenna_position_frame,
        default_antenna_diameter,
        nantenna,
        antenna,
        header
    );
    free(name);
    free(antenna_position_frame);
    for (int i = 0; i < nantenna; i++) {
        free(antenna[i].name);
    }
    free(antenna);
}

/*
* Parse the observation_info toml file, populating the relevant fields in `header`.
* The 'tests/obsinfo.toml' is exemplary.
*
* "input_map" is accessed and passed to `UVH5toml_parse_input_map`,
* after some preliminary calculation of the number of polarisations present.
* From the length of the "input_map" and inferred number of pol's, the 
* `header->Nants_data` adn `header->Nbls` are populated, and associated pointer-arrays
* allocated.
* This also populates `header->polarization_array`.
* Finally, `UVH5toml_parse_input_map` is called
*/
void UVH5toml_parse_observation_info(const char* file_path, UVH5_header_t* header) {
    FILE* fp;
    char errbuf[200];

    fp = fopen(file_path, "r");
    if (!fp) {
        UVH5print_error(__FUNCTION__, "cannot open '%s' (%s)", file_path, strerror(errno));
        return;
    }

    toml_table_t* conf = toml_parse_file(fp, errbuf, sizeof(errbuf));
    fclose(fp);
    if (!conf) {
        UVH5print_error(__FUNCTION__, "cannot parse %s", errbuf);
        return;
    }
    
    toml_array_t* toml_input_mapping = toml_array_in(conf, "input_map");
    if(toml_input_mapping) {
        int Npol_ant = 0; // Initial assumption
        UVH5print_verbose(__FUNCTION__, "input_map length: %d", toml_array_nelem(toml_input_mapping));

        char *ant_name0 = NULL, *ant_name1 = NULL;
        char pols_ant[3] = {'\0'}; // at most 2 unique pols
        // iterate through the input_map entries which are string tuples of (antenna_name, polarisation)
        // assume that entries are grouped by antenna
        // iterate until the antenna name changes, at which point all of the polarisations should be collected in `pols_ant`.
        do {
            if (ant_name1 != NULL) {
                free(ant_name1);
            }
            toml_array_t* toml_input_ant_name_pol = toml_array_at(toml_input_mapping, Npol_ant);
            _UVH5toml_string_at(toml_input_ant_name_pol, 0, &ant_name1);
            _UVH5toml_nstring_at(toml_input_ant_name_pol, 1, pols_ant + Npol_ant, 1);
            if ( ant_name0 == NULL) {
                ant_name0 = malloc(strlen(ant_name1));
                memcpy(ant_name0, ant_name1, strlen(ant_name1)+1);
            }
        } while(strcmp(ant_name0, ant_name1) == 0 && ++Npol_ant <= 2);
        free(ant_name0);
        free(ant_name1);
        if (Npol_ant > 2) {
            UVH5print_error(__FUNCTION__, "The first antenna is repeated more than twice indicating more than 2 polarisaitons.");
            return;
        }
        pols_ant[Npol_ant] = '\0';

        UVH5print_info(__FUNCTION__, "Assuming every antenna has %d polarisations ('%s'), as the first does.", Npol_ant, pols_ant);
        int nantenna = toml_array_nelem(toml_input_mapping)/Npol_ant;
        UVH5print_info(__FUNCTION__, "\tLeads to Nants_data: %d.", nantenna);

        char** antenna_names = malloc(nantenna*sizeof(char*));
        toml_array_t* toml_input_ant_name_pol;
        for(int a = 0; a < nantenna; a++) {
            toml_input_ant_name_pol = toml_array_at(toml_input_mapping, Npol_ant*a);
            _UVH5toml_string_at(toml_input_ant_name_pol, 0, antenna_names+a);
        }
        
        UVH5set_observation_info(
            nantenna,
            antenna_names,
            pols_ant,
            header
        );

        for(int a = 0; a < nantenna; a++) {
            free(antenna_names[a]);
        }
        free(antenna_names);
    }
    else {
        UVH5print_error(__FUNCTION__, "cannot read location 'input_map'");
    }

    toml_free(conf);
}
