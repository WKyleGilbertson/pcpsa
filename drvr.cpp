#include <cstdio>
#include <cstdint>
#include "G2INIT.h"
#include "NCO.h"

using namespace std;

int main(int argc, char* argv[]) {
  uint8_t PRN = 9;
  int8_t bin = 0;
  int8_t num_c = 0, num_s = 0, line=0;
  int16_t codephase = 0;
  uint16_t width = 500;
  uint32_t idx = 0, NCO_IDX = 0;
  float RefFreq = 4.092e6;
  //float RefFreq = 16.368e6;
  float SampleFreq = 16.368e6;
  float s = 0.0, c = 0.0;
  FILE * OUT = fopen("NCO.bin", "w");

  NCO CARRNCO(5, SampleFreq);
  CARRNCO.SetFrequency(RefFreq + (bin * width));

  if (argc >= 2)
    PRN = atoi(argv[1]);
  if (argc >= 3)
    codephase = atoi(argv[2]);
  if (argc >= 4)
    bin = atoi(argv[3]);
  printf("PRN: %d CP:%d bin: %d\n",
         PRN, codephase, bin);
  G2INIT sv(PRN, codephase);
  CARRNCO.LoadCACODE(sv.CODE);
  
  //  for (idx = 0; idx < 1 << 4; idx++) {
  //    //num = round(n() * 18);
  //    printf("%10.6f\n", n());
  //  }   
/*  for (idx = 0; idx<1023; idx++) {
    if (idx % 20 == 0) printf("\n%2d: ", ++line);
    printf("%2d ", sv.CODE[idx]);
  } */ 

  for (idx = 0; idx<32; idx++) {
    NCO_IDX = CARRNCO.clk();
    s = CARRNCO.sine(NCO_IDX)*127;
    c = CARRNCO.cosine(NCO_IDX)*127;
    num_s = (int8_t)round(s);
    num_c = (int8_t)round(c);
    fputc(num_s, OUT);
    if (idx < 100) {
      if (idx % 6 == 0) printf("\n");
      printf("[%4d,%4d] ", num_s, num_c);
    }
  }
  printf("\n");
  fclose(OUT);
  return 0;
}
