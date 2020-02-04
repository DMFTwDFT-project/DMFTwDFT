#include <sys/types.h>
#include <sys/times.h>
#include <sys/time.h>
#include <stdio.h>
#include <unistd.h>

void cputim_(dsec)
double *dsec;
{
      struct tms buffer;

      times(&buffer);
      *dsec = (double)buffer.tms_utime/sysconf(_SC_CLK_TCK);
/*      fprintf(stderr,"inside cputim: %lf\n",dsec);*/
}


void walltim_(dsec)
double * dsec;
{       struct timeval tv;
        struct timezone tz;
        gettimeofday(&tv, &tz);
        *dsec = (double) tv.tv_sec + tv.tv_usec/1000000.0;
}

void cputim(dsec)
double *dsec;
{ cputim_(dsec);
}

void walltim(dsec)
double * dsec;
{ walltim_(dsec);
}

