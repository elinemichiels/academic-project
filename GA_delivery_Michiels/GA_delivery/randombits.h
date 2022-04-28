#define USHRT_WIDTH 16
#define UINT_WIDTH 32

/*********************************
 * ran1 from "Numerical Recipes" *
 * note #undef's at end of file  *
**********************************/
#define IA 16807
#define IM 2147483647L
#define AM (1.0/IM)
#define IQ 127773L
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran1(long *idum);
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

void randomize(void);
float uniform(void);
unsigned char random_bit(void);
unsigned UINTran(void);
unsigned short USHRTran(void);
unsigned long ULONGran(unsigned char j);
unsigned char UCHARran(void);


