#ifndef COMMON_H_
#define COMMON_H_

#define HEADER_LENGTH 10
#define READ_LENGHT 150
#define READS_PER_CHUNK 15000

#ifndef SDIV
  #define SDIV(x,y)(((x)+(y)-1)/(y))
#endif

#endif // COMMON_H_
