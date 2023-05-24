#include <stdio.h>
#include <stdint.h>

void otroSwap(int **a, int **b){
    uint64_t aux = (uint64_t) *a;
    *a = *b;
    *b = (int*)aux;
}

void do_other_thing(int **a, int **b){
    otroSwap(a, b);
}
void do_something(int **a, int **b){

    do_other_thing(a, b);

}
int main(int argc, char **argv){
    int test, test1;
    int *a, *b;
    a = &test;
    b = &test1;
    printf("%p, %p\n", a, b);
    do_something(&a, &b);
    printf("%p, %p\n", a, b);


    return 0;
}