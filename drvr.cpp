#include <cstdio>
#include <cstdint>
#include <cmath>
#include <complex>
#include <vector> 
#include "FftComplex.hpp"
#include "G2INIT.h"
#include "NCO.h"

using namespace std;

void stuffVector(vector<complex<double>> &vec, FILE *fp) {
  int8_t buffer[32736];
  size_t idx = 0;

  size_t bytesRead = fread(buffer, 1, sizeof(buffer), fp);

  for (idx=0; idx<16368; idx++) {
    double re = (double)buffer[2*idx];
    double im = (double)buffer[2*idx + 1];  
    vec[idx] = complex<double>(re, im);
  }
  for (idx = 16368; idx < 16384; idx++) {
    vec[idx] = complex<double>(0, 0);
  }
}

int main(int argc, char* argv[]) {
  std::vector<std::complex<double>> originalData(16384);
  std::vector<std::complex<double>> data(16384);
  std::vector<std::complex<double>> codeVec(16384);
  uint8_t PRN = 21;
  int8_t bin = 0;
  int8_t num_c = 0, num_s = 0, line=0;
  int16_t codephase = 0;
  uint16_t width = 500;
  uint32_t idx = 0, NCO_IDX = 0;
  float RefFreq = 4.092e6;
  float SampleFreq = 16.368e6;
  float s = 0.0, c = 0.0;
  FILE * IN = fopen("IF.bin", "rb");
  FILE * OUT = fopen("NCO.bin", "wb");

  if (argc >= 2)
    PRN = atoi(argv[1]);
  if (argc >= 3)
    bin = atoi(argv[2]);
  printf("PRN: %d bin: %d\n", PRN, bin);

  G2INIT sv(PRN, codephase);
  //NCO CARRNCO(5, SampleFreq);
  NCO CARRNCO(10, SampleFreq);
//  CARRNCO.SetFrequency(RefFreq + (bin * width));
  CARRNCO.LoadCACODE(sv.CODE);
  for(idx = 0; idx < 16384; idx++) {
    int chipIdx = (idx/16) % 1023;
    double codeVal = (double)sv.CODE[chipIdx];
    codeVec[idx] = complex<double>(codeVal, 0.0);
}
Fft::transform(codeVec, false); // Forward FFT of the code vector, to be used in correlation in the frequency domain  
for (auto &val : codeVec) val = conj(val); // Take the complex conjugate of the code FFT for correlation

  //  for (idx = 0; idx < 1 << 4; idx++) {
  //    //num = round(n() * 18);
  //    printf("%10.6f\n", n());
  //  }   
/*  for (idx = 0; idx<1023; idx++) {
    if (idx % 20 == 0) printf("\n%2d: ", ++line);
    printf("%2d ", sv.CODE[idx]);
  } */ 

  stuffVector(originalData, IN);  
  fclose(IN);

  for (bin = -20; bin <= 20; bin++) {
    data = originalData;
    CARRNCO.SetFrequency(RefFreq + (bin * width));
  // 1. Process the vector data with the NCO
for (idx = 0; idx < data.size(); idx++) {
    // Get the next NCO index and look up the complex LO
    uint32_t NCO_IDX = CARRNCO.clk();
    complex<double> lo(CARRNCO.cosine(NCO_IDX), CARRNCO.sine(NCO_IDX));

    // Perform the complex multiply: (I + jQ) * (cos + jsin)
    data[idx] *= lo;
}
// Now data is baseband (or at your target offset) and ready for FFT
Fft::transform(data, false);// False means not inverse FFT, i.e. forward FFT 

for (idx = 0; idx < 16384; idx++) {
    // Perform the correlation in the frequency domain: element-wise multiply with the conjugate of the code FFT
    data[idx] *= codeVec[idx]; 
} 

Fft::transform(data, true); // Inverse FFT to get the correlation result in the time domain

double maxMag = 0;
int peakIndex = 0;
double sumPower = 0;

for (int idx =0; idx < data.size(); idx++) {
    double mag = std::abs(data[idx])/16384.0; 
    sumPower += mag*mag; // Sum of Squares (Power)

    if (mag > maxMag) {
        maxMag = mag;
        peakIndex = idx;
    }
}

double peakPower = maxMag * maxMag;
double avgNoisePower = (sumPower -peakPower) / (data.size() -1);
double snr_power = 10.0 * log10(peakPower / avgNoisePower);

printf("bin %3d Peak found at index: %5d with magnitude: %10.1f %4d SNR: %5.2f dB\n", 
        bin, peakIndex, maxMag, peakIndex / 16, snr_power);
  }
/*  for (idx = 0; idx<32; idx++) {
    NCO_IDX = CARRNCO.clk();
    s = CARRNCO.sine(NCO_IDX);
    c = CARRNCO.cosine(NCO_IDX);
    num_s = (int8_t)round(s);
    num_c = (int8_t)round(c);
    fputc(num_s, OUT);
    if (idx < 100) {
      if (idx % 5 == 0) printf("\n");
      //printf("[%4d,%4d] ", num_s, num_c);
      printf("[%6.3f,%6.3f] ", s, c);
    }
  } */
  printf("\n");
  fclose(OUT);
  return 0;
}
