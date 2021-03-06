#include "definitions.h"

//extern int sdatio_debug;

extern "C" void sdatio_init(struct sdatio_file * sfile, char * fname);
extern "C" void sdatio_free(struct sdatio_file * sfile);
/* Open a new datafile for writing. fname is the name of the file 
 * The stuct sfile is used to store the state information
 * of the file.*/
extern "C" void sdatio_create_file(struct sdatio_file * sfile);

/* Create a new dimension in the file sfile. Dimension names must
 * be a single letter. */
extern "C" void sdatio_add_dimension(struct sdatio_file * sfile, 
													 char * dimension_name, 
													 int size,
													 char * description,
													 char * units);

/* Print out a nice list of all the dimensions defined so far*/
extern "C" void sdatio_print_dimensions(struct sdatio_file * sfile);


/* Close the file and free all memory associated with sfile*/
extern "C" void sdatio_close(struct sdatio_file * sfile);

/* Ensure all variables are written to disk in case of crashes*/
extern "C" void sdatio_sync(struct sdatio_file * sfile);

/* Define a variable in the given file. Dimension list 
 * is a character string listing (in order) the dimension names
 * (which are all single characters) e.g. "xyx".*/
extern "C" void sdatio_create_variable(struct sdatio_file * sfile,
														int variable_type,
														char * variable_name,
														char * dimension_list,
														char * description,
														char * units);

/* Write to the given variable. address should be the address of the start of the array */
extern "C" void sdatio_write_variable(struct sdatio_file * sfile, char * variable_name, void * address);
/* REad the given variable. address should be the address of the start of the array */
extern "C" void sdatio_read_variable(struct sdatio_file * sfile, char * variable_name, void * address);

/* Write to the given variable. address should be the address of the start of the array. Indexes should be an array the same size as the number of dimensions of the variable. Using the second form is quicker as the first form requires a search for the variable at every write*/
extern "C" void sdatio_write_variable_at_index(struct sdatio_file * sfile, char * variable_name, int * indexes, void * address);
extern "C" void sdatio_write_variable_at_index_fast(struct sdatio_file * sfile, struct sdatio_variable * svar, int * indexes, void * address);

/* Print out a nice list of all the variables defined so far*/
extern "C" void sdatio_print_variables(struct sdatio_file * sfile);

/* Return a pointer the struct containing all the metadata of the given variable */
extern "C" struct sdatio_variable * sdatio_find_variable(struct sdatio_file * sfile, char * variable_name);

/* Return a pointer the struct containing all the metadata of the given dimension */
extern "C" struct sdatio_dimension * sdatio_find_dimension(struct sdatio_file * sfile, char * dimension_name);

/* Increment the start of the specified infinite dimension */
extern "C" void sdatio_increment_start(struct sdatio_file * sfile, char * dimension_name);

extern "C" void sdatio_set_dimension_start(struct sdatio_file * sfile, char * dimension_name, int start);

extern "C" void sdatio_set_count(struct sdatio_file * sfile, char * variable_name, char * dimension_name, int * count);

/* Returns 1 if the given variable exists, 0 otherwise */
extern "C" int sdatio_variable_exists(struct sdatio_file * sfile, char * variable_name);

/* Open an existing file, reading all dimension and variable metadata */
extern "C" void sdatio_open_file(struct sdatio_file * sfile );
