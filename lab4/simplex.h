/* SPDX-License-Identifier: GPL-2.0 */
/*
 * simplex.h - библиотека для решения задачи линейного программирования
 *            симплекс-методом.
 *
 * Таблица хранится как плоский массив double размера rows * cols,
 * элементы уложены построчно: data[i * cols + j].
 */

#ifndef _SIMPLEX_H
#define _SIMPLEX_H

#include <stddef.h>
#include <stdio.h>

/*
 * Порог для сравнения вещественных чисел. Всё, что по модулю меньше,
 * считается нулём. Можно переопределить при сборке:
 *   cc -DSIMPLEX_EPS=1e-7 ...
 */
#ifndef SIMPLEX_EPS
#define SIMPLEX_EPS 1e-9
#endif

/*
 * Максимум итераций симплекс-метода. Защита от зацикливания.
 */
#ifndef SIMPLEX_ITER_LIMIT
#define SIMPLEX_ITER_LIMIT 10000
#endif

/**
 * enum simplex_status - результат решения задачи.
 * @SIMPLEX_OK: оптимальное решение найдено.
 * @SIMPLEX_UNBOUNDED: целевая функция не ограничена, решения нет.
 * @SIMPLEX_INCONSISTENT: система ограничений несовместна.
 * @SIMPLEX_MAX_ITER: превышен лимит итераций (возможное зацикливание).
 * @SIMPLEX_INVALID: некорректные аргументы (NULL, пустая таблица).
 */
enum simplex_status {
	SIMPLEX_OK = 0,
	SIMPLEX_UNBOUNDED,
	SIMPLEX_INCONSISTENT,
	SIMPLEX_MAX_ITER,
	SIMPLEX_INVALID,
};

/*
 * Непрозрачный тип. Пользователь библиотеки не должен лазить внутрь —
 * для этого есть функции доступа.
 */
struct simplex_table;

/**
 * simplex_table_create() - создать пустую таблицу.
 * @rows: число строк матрицы.
 * @cols: число столбцов матрицы.
 *
 * Выделяет память под структуру и матрицу rows x cols. Содержимое
 * не инициализируется — для заполнения используйте simplex_table_init().
 *
 * Return: указатель на таблицу или NULL при ошибке (нулевые размеры,
 *         нехватка памяти).
 */
struct simplex_table *simplex_table_create(size_t rows, size_t cols);

/**
 * simplex_table_destroy() - освободить таблицу.
 * @table: адрес указателя на таблицу.
 *
 * Освобождает матрицу и саму структуру, зануляет *table.
 * Безопасно вызывается с *table == NULL.
 */
void simplex_table_destroy(struct simplex_table **table);

/**
 * simplex_table_init() - заполнить таблицу значениями.
 * @table: таблица, созданная simplex_table_create().
 * @values: плоский массив длиной rows * cols, уложенный построчно.
 *
 * Копирует values в матрицу таблицы.
 *
 * Return: 0 при успехе, -1 если table или values равны NULL.
 */
int simplex_table_init(struct simplex_table *table, const double *values);

/**
 * simplex_table_print() - напечатать таблицу.
 * @table: таблица.
 * @stream: куда печатать (обычно stdout или stderr).
 *
 * Return: 0 при успехе, -1 если table или stream равны NULL.
 */
int simplex_table_print(const struct simplex_table *table, FILE *stream);

/**
 * simplex_table_rows() - число строк таблицы.
 * @table: таблица.
 *
 * Return: число строк или 0, если table == NULL.
 */
size_t simplex_table_rows(const struct simplex_table *table);

/**
 * simplex_table_cols() - число столбцов таблицы.
 * @table: таблица.
 *
 * Return: число столбцов или 0, если table == NULL.
 */
size_t simplex_table_cols(const struct simplex_table *table);

/**
 * simplex_table_get() - прочитать элемент.
 * @table: таблица.
 * @i: строка (0-based).
 * @j: столбец (0-based).
 * @out: сюда записать значение.
 *
 * Return: 0 при успехе, -1 если table == NULL, out == NULL
 *         или (i, j) вне границ.
 */
int simplex_table_get(const struct simplex_table *table,
		      size_t i, size_t j, double *out);

/**
 * simplex_table_set() - записать элемент.
 * @table: таблица.
 * @i: строка (0-based).
 * @j: столбец (0-based).
 * @value: значение.
 *
 * Return: 0 при успехе, -1 если table == NULL
 *         или (i, j) вне границ.
 */
int simplex_table_set(struct simplex_table *table,
		      size_t i, size_t j, double value);

/**
 * simplex_table_copy() - сделать копию таблицы.
 * @src: исходная таблица.
 *
 * Return: указатель на новую таблицу или NULL при ошибке.
 */
struct simplex_table *simplex_table_copy(const struct simplex_table *src);

/**
 * simplex_solve() - решить ЗЛП табличным симплекс-методом.
 * @table: симплекс-таблица размера (m + 1) x (n + 1), где строка 0 -
 *         целевая функция, строки 1..m - ограничения, столбцы 0..n-1 -
 *         свободные переменные, последний столбец - свободные члены.
 *
 * Таблица модифицируется на месте. При статусе SIMPLEX_OK содержит
 * оптимальную таблицу: свободные члены в последнем столбце - значения
 * базисных переменных, правый нижний угол - значение ЦФ (с учётом
 * знака, см. методичку).
 *
 * Return: статус решения.
 */
enum simplex_status simplex_solve(struct simplex_table *table);

#endif /* _SIMPLEX_H */