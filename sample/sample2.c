#include <stdio.h>
#include <stdlib.h>

#include <tepla/ec.h>

int main(void)
{
    Field f;

    field_init(f, "bn254_fpa");
    field_clear(f);

    return 0;
}
