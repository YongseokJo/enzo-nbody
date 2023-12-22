#include <stdio.h>
void auto_show_version(FILE *fp) {
   fprintf (fp,"make[1]: Warning: File `Make.config.targets' has modification time 117 s in the future\n");
   fprintf (fp,"\n");
   fprintf (fp,"Git Branch   HEAD\n");
   fprintf (fp,"Git Revision 60a60f5375bbc3a626739fb9d4ceb41ad7467850\n");
   fprintf (fp,"\n");
   fprintf (fp,"make[1]: warning:  Clock skew detected.  Your build may be incomplete.\n");
}
