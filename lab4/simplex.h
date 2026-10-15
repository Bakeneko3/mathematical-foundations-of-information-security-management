#ifndef SIMPLEX_H
#define SIMPLEX_H

#include <stddef.h>

#ifndef SIMPLEX_EPS
#define SIMPLEX_EPS 1e-9
#endif

#ifndef SIMPLEX_ITER_LIMIT
#define SIMPLEX_ITER_LIMIT 10000
#endif

struct simplex_table;

struct simplex_table *simplex_table_create(int rows, int cols,
										   const double *values);

void simplex_table_destroy(struct simplex_table *table);

size_t simplex_table_rows(const struct simplex_table *table);

size_t simplex_table_cols(const struct simplex_table *table);

int simplex_table_get(const struct simplex_table *table,
					  size_t row, size_t col, double *out);

int simplex_table_set(struct simplex_table *table,
					  size_t row, size_t col, double value);

struct simplex_table *simplex_table_copy(const struct simplex_table *src);

#endif /* SIMPLEX_H */