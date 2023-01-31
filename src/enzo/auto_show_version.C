#include <stdio.h>
void auto_show_version(FILE *fp) {
   fprintf (fp,"make[1]: Warning: File `DEPEND' has modification time 96 s in the future\n");
   fprintf (fp,"\n");
   fprintf (fp,"Git Branch   master\n");
   fprintf (fp,"Git Revision 6a721658990deb6c0c1c48a87abd52d997c8801d\n");
   fprintf (fp,"\n");
   fprintf (fp,"make[1]: warning:  Clock skew detected.  Your build may be incomplete.\n");
}
