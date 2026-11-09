/* The low-level routines must be exported from a Windows DLL (issue #91):
   this program fails to link if one of them is not. */
#include <stdio.h>
#include <cminpack.h>
#include <minpack.h>

typedef void (*fp)(void);

int main(void)
{
    fp f[] = {
        (fp)__cminpack_func__(dogleg), (fp)__cminpack_func__(qrfac),
        (fp)__cminpack_func__(qrsolv), (fp)__cminpack_func__(qform),
        (fp)__cminpack_func__(r1updt), (fp)__cminpack_func__(r1mpyq),
        (fp)__cminpack_func__(lmpar), (fp)__cminpack_func__(rwupdt),
        (fp)__minpack_func__(dogleg), (fp)__minpack_func__(qrfac),
        (fp)__minpack_func__(qrsolv), (fp)__minpack_func__(qform),
        (fp)__minpack_func__(r1updt), (fp)__minpack_func__(r1mpyq),
        (fp)__minpack_func__(lmpar), (fp)__minpack_func__(rwupdt)
    };
    int i;
    for (i = 0; i < 16; ++i) {
        if (f[i] == NULL) {
            printf("routine %d is null\n", i);
            return 1;
        }
    }
    return 0;
}
