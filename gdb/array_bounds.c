#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "show.h"

int main(int argc, char *argv[]) {
    double a[10], sum;
    int i;
    for (i = 0; i <= 10; i++)
        a[i] = sqrt(i);
    for (i = 0; i < 10; i++)
        sum += a[i];
    printResult("sum", sum);
    return EXIT_SUCCESS;
}

