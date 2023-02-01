#include <stdio.h>
void auto_show_version(FILE *fp) {
   fprintf (fp,"make[1]: Warning: File `DEPEND' has modification time 80 s in the future\n");
   fprintf (fp,"\n");
   fprintf (fp,"Git Branch   master\n");
   fprintf (fp,"Git Revision d167f916f730ee8f4578336bcee6c0d44a8d1f97\n");
   fprintf (fp,"\n");
   fprintf (fp,"make[1]: warning:  Clock skew detected.  Your build may be incomplete.\n");
}
