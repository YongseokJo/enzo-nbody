#include <stdio.h>
void auto_show_version(FILE *fp) {
   fprintf (fp,"make[1]: Warning: File `DEPEND' has modification time 100 s in the future\n");
   fprintf (fp,"\n");
   fprintf (fp,"Git Branch   master\n");
   fprintf (fp,"Git Revision 02b55e6dc10d47d485a7c84ff1f0209f0a322e41\n");
   fprintf (fp,"\n");
   fprintf (fp,"make[1]: warning:  Clock skew detected.  Your build may be incomplete.\n");
}
