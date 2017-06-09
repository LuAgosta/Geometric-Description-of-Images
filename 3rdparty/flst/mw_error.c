#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

void mwerror(int code, int exit_code, char *fmt) {
  printf("%s\n", fmt);
  exit(exit_code);
}
