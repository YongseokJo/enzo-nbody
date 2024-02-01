// header file cuda_functions.h
// that includes the necessary cuda functions to communicate with GPU

extern "C" {
  void gpunb_devinit_(int *irank);
  void gpunb_open_(int *nbmax, int *irank);
  void gpunb_close_();
  void gpunb_send_(int *nj, double mj[], double xj[][3], double vj[][3]);
  void gpunb_regf_(int *ni, double h2[], double dtr[], double xi[][3],
                   double vi[][3], double acc[][3], double jrk[][3],
                   double pot[], int *lmax, int *nnbmax, int *list, int *m_flag);
  void gpunb_profile_(int *irank);
}

// end of header file