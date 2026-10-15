/* SPDX-License-Identifier: GPL-2.0 */
/*
 * simplex.c - реализация библиотеки симплекс-метода.
 */

#include <stdlib.h>
#include <string.h>

#include "simplex.h"

struct simplex_table {
	size_t rows;
	size_t cols;

	double *data;	/* rows * cols, построчно */
};

struct simplex_table *simplex_table_create(size_t rows, size_t cols)
{
	struct simplex_table *table;

	if (rows == 0 || cols == 0)
		return NULL;

	table = malloc(sizeof(*table));
	if (!table)
		return NULL;

	table->data = malloc(rows * cols * sizeof(*table->data));
	if (!table->data) {
		free(table);
		return NULL;
	}

	table->rows = rows;
	table->cols = cols;

	return table;
}

void simplex_table_destroy(struct simplex_table **table)
{
	if (!table || !*table)
		return;

	free((*table)->data);
	free(*table);
	*table = NULL;
}

int simplex_table_init(struct simplex_table *table, const double *values)
{
	if (!table || !values)
		return -1;

	memcpy(table->data, values,
	       table->rows * table->cols * sizeof(*table->data));

	return 0;
}

int simplex_table_print(const struct simplex_table *table, FILE *stream)
{
	size_t i, j;

	if (!table || !stream)
		return -1;

	for (i = 0; i < table->rows; i++) {
		for (j = 0; j < table->cols; j++)
			fprintf(stream, "%8.2f ", table->data[i * table->cols + j]);
		fputc('\n', stream);
	}

	return 0;
}

size_t simplex_table_rows(const struct simplex_table *table)
{
	if (!table)
		return 0;

	return table->rows;
}

size_t simplex_table_cols(const struct simplex_table *table)
{
	if (!table)
		return 0;

	return table->cols;
}

int simplex_table_get(const struct simplex_table *table,
		      size_t i, size_t j, double *out)
{
	if (!table || !out)
		return -1;

	if (i >= table->rows || j >= table->cols)
		return -1;

	*out = table->data[i * table->cols + j];

	return 0;
}

int simplex_table_set(struct simplex_table *table,
		      size_t i, size_t j, double value)
{
	if (!table)
		return -1;

	if (i >= table->rows || j >= table->cols)
		return -1;

	table->data[i * table->cols + j] = value;

	return 0;
}

struct simplex_table *simplex_table_copy(const struct simplex_table *src)
{
	struct simplex_table *dst;

	if (!src)
		return NULL;

	dst = simplex_table_create(src->rows, src->cols);
	if (!dst)
		return NULL;

	if (simplex_table_init(dst, src->data) != 0) {
		simplex_table_destroy(&dst);
		return NULL;
	}

	return dst;
}

/* ---------- внутренние хелперы ---------- */

static int simplex_is_zero(double x)
{
	return x > -SIMPLEX_EPS && x < SIMPLEX_EPS;
}

static int simplex_is_negative(double x)
{
	return x < -SIMPLEX_EPS;
}

static int simplex_is_positive(double x)
{
	return x > SIMPLEX_EPS;
}

/*
 * Проверка совместности: несовместна, если (кроме строки ЦФ) есть
 * строка с отрицательным свободным членом и всеми неотрицательными
 * коэффициентами при свободных переменных.
 */
static int simplex_check_consistent(const struct simplex_table *table)	
{
	size_t i, j;
	double b, a;

	for (i = 1; i < table->rows; i++) {
		simplex_table_get(table, i, table->cols - 1, &b);
		if (!simplex_is_negative(b))
			continue;

		for (j = 0; j < table->cols - 1; j++) {
			simplex_table_get(table, i, j, &a);
			if (simplex_is_negative(a))
				break;
		}
		if (j == table->cols - 1)
			return 0;	/* все a >= 0, а b < 0 — несовместна */
	}

	return 1;
}

/*
 * Проверка ограниченности: неограниченна, если в строке ЦФ есть
 * отрицательный элемент, в столбце которого нет ни одного
 * положительного элемента среди строк ограничений.
 */
static int simplex_check_bounded(const struct simplex_table *table)
{
	size_t i, j;
	double a0j, aij;

	for (j = 0; j < table->cols - 1; j++) {
		simplex_table_get(table, 0, j, &a0j);
		if (!simplex_is_negative(a0j))
			continue;

		for (i = 1; i < table->rows; i++) {
			simplex_table_get(table, i, j, &aij);
			if (simplex_is_positive(aij))
				break;
		}
		if (i == table->rows)
			return 0;	/* столбец весь <= 0 — неограниченна */
	}

	return 1;
}

/*
 * Проверка оптимальности: в строке ЦФ нет отрицательных элементов
 * (свободный член не учитываем).
 */
static int simplex_check_optimal(const struct simplex_table *table)
{
	size_t j;
	double a0j;

	for (j = 0; j < table->cols - 1; j++) {
		simplex_table_get(table, 0, j, &a0j);
		if (simplex_is_negative(a0j))
			return 0;
	}

	return 1;
}

/*
 * Поиск разрешающего элемента.
 * @row_out, @col_out: сюда записать координаты (k, l).
 *
 * Return: 0 если найден, -1 если нет (не должно случаться при
 *         пройденных проверках на ограниченность/оптимальность).
 */
static int simplex_find_pivot(const struct simplex_table *table,
			      size_t *row_out, size_t *col_out)
{
	size_t i, j, l = 0, k = 0;
	double a0j, min_a0j;
	double aij, bi, ratio, min_ratio;
	int found_col = 0, found_row = 0;

	/* Разрешающий столбец: минимальный отрицательный a0j. */
	for (j = 0; j < table->cols - 1; j++) {
		simplex_table_get(table, 0, j, &a0j);
		if (!simplex_is_negative(a0j))
			continue;

		if (!found_col || a0j < min_a0j || (a0j == min_a0j && rand() % 2)) { 
			min_a0j = a0j;
			l = j;
			found_col = 1;
		}
	}

	if (!found_col)
		return -1;

	/* Разрешающая строка: минимум b_i / a_il > 0. */
	for (i = 1; i < table->rows; i++) {
		simplex_table_get(table, i, l, &aij);
		if (simplex_is_zero(aij))
			continue;

		simplex_table_get(table, i, table->cols - 1, &bi);
		ratio = bi / aij;
		if (!simplex_is_positive(ratio))
			continue;

		if (!found_row || ratio < min_ratio || (ratio == min_ratio && rand() % 2)) {
			min_ratio = ratio;
			k = i;
			found_row = 1;
		}
	}

	if (!found_row)
		return -1;

	*row_out = k;
	*col_out = l;

	return 0;
}

/*
 * Пересчёт таблицы по разрешающему элементу (k, l).
 * Использует временный буфер, чтобы не портить исходные значения
 * в процессе вычислений.
 *
 * Return: 0 при успехе, -1 при ошибке выделения памяти.
 */
static int simplex_pivot(struct simplex_table *table, size_t k, size_t l)
{
	size_t i, j;
	size_t rows = table->rows;
	size_t cols = table->cols;
	double *old;
	double akl, aki, ajl, aji;

	old = malloc(rows * cols * sizeof(*old));
	if (!old)
		return -1;

	memcpy(old, table->data, rows * cols * sizeof(*old));

	akl = old[k * cols + l];

	/* Все элементы новой таблицы по единой формуле. */
	for (i = 0; i < rows; i++) {
		for (j = 0; j < cols; j++) {
			if (i == k && j == l) {
				table->data[i * cols + j] = 1.0 / akl;
			} else if (i == k) {
				aki = old[k * cols + j];
				table->data[i * cols + j] = aki / akl;
			} else if (j == l) {
				ajl = old[i * cols + l];
				table->data[i * cols + j] = -ajl / akl;
			} else {
				aji = old[i * cols + j];
				aki = old[k * cols + j];
				ajl = old[i * cols + l];
				table->data[i * cols + j] =
					aji - aki * ajl / akl;
			}
		}
	}

	free(old);

	return 0;
}

/* ---------- публичная функция ---------- */

enum simplex_status simplex_solve(struct simplex_table *table)
{
	size_t k, l;
	unsigned int iter;
	int rc;

	if (!table || table->rows < 2 || table->cols < 2)
		return SIMPLEX_INVALID;

	for (iter = 0; iter < SIMPLEX_ITER_LIMIT; iter++) {
		if (!simplex_check_consistent(table))
			return SIMPLEX_INCONSISTENT;

		if (simplex_check_optimal(table))
			return SIMPLEX_OK;

		if (!simplex_check_bounded(table))
			return SIMPLEX_UNBOUNDED;

		rc = simplex_find_pivot(table, &k, &l);
		if (rc != 0)
			return SIMPLEX_INVALID;	/* не должно случаться */

		rc = simplex_pivot(table, k, l);
		if (rc != 0)
			return SIMPLEX_INVALID;	/* malloc сдох */
		
		simplex_table_print(table, stdout);
	}

	return SIMPLEX_MAX_ITER;
}
