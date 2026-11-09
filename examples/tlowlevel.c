/* The low-level routines must be exported from a Windows DLL (issue #91):
   this program fails to link if one of them is not. */
#include <stdio.h>
#include <cminpack.h>
#include <minpack.h>

typedef void (*fp)(void);

int main(void)
{
    fp f[] = { (fp)dogleg, (fp)qrfac, (fp)qrsolv, (fp)qform, (fp)r1updt,
               (fp)r1mpyq, (fp)lmpar, (fp)rwupdt,
               (fp)dogleg_, (fp)qrfac_, (fp)qrsolv_, (fp)qform_, (fp)r1updt_,
               (fp)r1mpyq_, (fp)lmpar_, (fp)rwupdt_ };
    int i;
    for (i = 0; i < 16; ++i) {
        if (f[i] == NULL) {
            printf("routine %d is null\n", i);
            return 1;
        }
    }
    return 0;
}
