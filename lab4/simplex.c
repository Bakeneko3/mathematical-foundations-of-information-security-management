#include <stdlib.h>
#include <string.h>

#include "simplex.h"

static inline int simplex_is_zero(double x)
{
	return x > -SIMPLEX_EPS && x < SIMPLEX_EPS;
}

static inline int simplex_is_negative(double x)
{
	return x < -SIMPLEX_EPS;
}

static inline int simplex_is_positive(double x)
{
	return x > SIMPLEX_EPS;
}

struct simplex_table {
	size_t rows;
	size_t cols;

	double *data;
};

struct simplex_table *simplex_table_create(int rows, int cols,
										   const double *values)
{
	struct simplex_table *table;
	size_t size;

	if (rows < 2 || cols < 2 || !values) {
		return NULL;
	}

	table = malloc(sizeof(*table));
	if (!table) {
		return NULL;
	}

	table->rows = (size_t)rows;
	table->cols = (size_t)cols;

	size = table->rows * table->cols * sizeof(*table->data);

	table->data = malloc(size);
	if (!table->data) {
		free(table);
		return NULL;
	}

	memcpy(table->data, values, size);

	return table;
}

void simplex_table_destroy(struct simplex_table *table)
{
	if (!table) {
		return;
	}

	free(table->data);
	free(table);
}

size_t simplex_table_rows(const struct simplex_table *table)
{
	if (!table) {
		return 0;
	}

	return table->rows;
}

size_t simplex_table_cols(const struct simplex_table *table)
{
	if (!table) {
		return 0;
	}

	return table->cols;
}

int simplex_table_get(const struct simplex_table *table,
					  size_t row, size_t col, double *out)
{
	if (!table || !out) {
		return -1;
	}

	if (row >= table->rows || col >= table->cols) {
		return -1;
	}

	*out = table->data[row * table->cols + col];

	return 0;
}

int simplex_table_set(struct simplex_table *table,
					  size_t row, size_t col, double value)
{
	if (!table) {
		return -1;
	}

	if (row >= table->rows || col >= table->cols) {
		return -1;
	}

	table->data[row * table->cols + col] = value;

	return 0;
}

struct simplex_table *simplex_table_copy(const struct simplex_table *src)
{
	struct simplex_table *dest;

	if (!src) {
		return NULL;
	}

	dest = simplex_table_create((int)src->rows, (int)src->cols, src->data);
	if (!dest) {
		return NULL;
	}

	return dest;
}

static int simplex_is_consistent(const struct simplex_table *table)
{
	size_t j, i;
	double b, a;

	if (!table) {
		return -1;
	}

	for (j = 1; j < table->rows; j++) {
		simplex_table_get(table, j, table->cols - 1, &b);

		if (!simplex_is_negative(b)) {
			continue;
		}

		for (i = 0; i < table->cols - 1; i++) {
			simplex_table_get(table, j, i, &a);

			if (simplex_is_negative(a)) {
				break;
			}
		}

		if (i == table->cols - 1) {
			return 0;
		}
	}

	return 1;
}

static int simplex_is_bounded(const struct simplex_table *table)
{
	size_t i, j;
	double a;

	if (!table) {
		return -1;
	}

	for (i = 0; i < table->cols - 1; i++) {
		simplex_table_get(table, 0, i, &a);

		if (!simplex_is_negative(a)) {
			continue;
		}

		for (j = 1; j < table->rows; j++) {
			simplex_table_get(table, j, i, &a);

			if (simplex_is_positive(a)) {
				break;
			}
		}

		if (j == table->rows) {
			return 0;
		}
	}

	return 1;
}

static int simplex_is_feasible(const struct simplex_table *table)
{
	size_t j;
	double b;

	if (!table) {
		return -1;
	}

	for (j = 1; j < table->rows; j++) {
		simplex_table_get(table, j, table->cols - 1, &b);

		if (simplex_is_negative(b)) {
			return 0;
		}
	}

	return 1;
}

static int simplex_is_optimal(const struct simplex_table *table)
{
	size_t i;
	double a;

	if (!table) {
		return -1;
	}

	for (i = 0; i < table->cols - 1; i++) {
		simplex_table_get(table, 0, i, &a);

		if (simplex_is_negative(a)) {
			return 0;
		}
	}

	return 1;
}

static int simplex_find_pivot_step3(const struct simplex_table *table,
				    size_t *k, size_t *l)
{
	size_t j, i;
	size_t pool;
	double b, a;
	double min_b, min_a;

	if (!table || !k || !l) {
		return -1;
	}

	min_b = 0;
	pool = 0;

	for (j = 1; j < table->rows; j++) {
		simplex_table_get(table, j, table->cols - 1, &b);

		if (j == 1 || simplex_is_negative(b - min_b)) {
			min_b = b;
			*k = j;
			pool = 1;
		} else if (simplex_is_zero(b - min_b)) {
			pool++;
			if (rand() % pool == 0) {
				*k = j;
			}
		}
	}

	if (!simplex_is_negative(min_b)) {
		return -1;
	}

	min_a = 0;
	pool = 0;

	for (i = 0; i < table->cols - 1; i++) {
		simplex_table_get(table, *k, i, &a);

		if (i == 0 || simplex_is_negative(a - min_a)) {
			min_a = a;
			*l = i;
			pool = 1;
		} else if (simplex_is_zero(a - min_a)) {
			pool++;
			if (rand() % pool == 0) {
				*l = i;
			}
		}
	}

	return 0;
}


