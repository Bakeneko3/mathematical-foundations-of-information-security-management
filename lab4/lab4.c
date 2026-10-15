/* SPDX-License-Identifier: GPL-2.0 */
/*
 * lab4.c - точка входа. Данные задаются прямо здесь.
 */

#include <stdio.h>

#include "simplex.h"

static const double source[4 * 5] = {
	/*  a00  a01  a02  a03  b0  */
	    -14, -10, -14, -11,  0,
	/*  a10  a11  a12  a13  b1  */
	      4,   2,   2,   3, 35,
    /*  a20  a21  a22  a23  b2  */
	      1,   1,   2,   3, 30,
    /*  a30  a31  a32  a33  b3  */
          3,   1,   2,   1, 40,
};

int main(void)
{
	struct simplex_table *t1;
    struct simplex_table *t2;
    int err;

	t1 = simplex_table_create(4, 5);
	if (!t1) {
		fprintf(stderr, "simplex_table_create failed\n");
		return 1;
	}

	if (simplex_table_init(t1, source) != 0) {
		fprintf(stderr, "simplex_table_init failed\n");
		simplex_table_destroy(&t1);
		return 1;
	}

	if (simplex_table_print(t1, stdout) != 0) {
		fprintf(stderr, "simplex_table_print failed\n");
		simplex_table_destroy(&t1);
		return 1;
	}

    t2 = simplex_table_copy(t1);
    if (!t2) {
        fprintf(stderr, "simplex_table_copy failed\n");
        simplex_table_destroy(&t1);
    }

    if (simplex_table_print(t2, stdout) != 0) {
		fprintf(stderr, "simplex_table_print failed\n");
		simplex_table_destroy(&t1);
        simplex_table_destroy(&t2);
		return 1;
    }

    printf("\nsolve:\n\n");
    err = simplex_solve(t1);
    printf("err: %d\n", err);

	simplex_table_destroy(&t1);
    simplex_table_destroy(&t2);
	return 0;
}

