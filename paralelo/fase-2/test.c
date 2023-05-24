#include <stdio.h>
#include <stdint.h>

void swap(int **a, int **b){
    uint64_t x = (uint64_t) *a;
    uint64_t y = (uint64_t) *b;
    x ^= y;
    y ^= x;
    x ^= y;
    
    *a = (int*) x;
    *b = (int*) y;
}

void otroSwap(int **a, int **b){
    uint64_t aux = *a;
    *a = *b;
    *b = aux;
}

void swapNormal(int *a, int *b){
    uint64_t x = (uint64_t) a;
    uint64_t y = (uint64_t) b;
    x ^= y;
    y ^= x;
    x ^= y;
    
    a = (int*) x;
    b = (int*) y;
}

int main(int argc, char **argv){
    int test, test1;
    int *a, *b;
    test = 0;
    test1 = 1;

    a = &test;
    b = &test1;

    printf("test: %p, test1: %p\na: %p, b: %p\n\n", &test, &test1, a, b);

    swap(&a, &b);
    printf("test: %p, test1: %p\na: %p, b: %p\n\n", &test, &test1, a, b);

    swapNormal(a, b);
    printf("test: %p, test1: %p\na: %p, b: %p\n\n", &test, &test1, a, b);

    otroSwap(&a, &b);
    printf("test: %p, test1: %p\na: %p, b: %p\n\n", &test, &test1, a, b);

    return 0;
}