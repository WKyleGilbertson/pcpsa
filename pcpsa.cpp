#include <cstdio>
#include <cstdint>
#include <cmath>
#include <complex>
#include <vector> 
#include <iostream>
#include <string>
#include <algorithm>
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

struct AcqResult {
    int bin;
    int peakIndex;
    double peakMagnitude;
    double snr;
};

class AcquisitionEngine {
  private:
   const size_t N = 16384;
   //const double SampleFreq = 16.368e6;
   double m_sampleFreq;
//   const double RefFreq = 4.092e6;
  // Pre-computed PRN code FFTs for each PRN (1-32)
  // Store them to avoid re-calculaing inside search loops
  std::vector<std::complex<double>> codeFfts[33];

  public:
    AcquisitionEngine(double sampleFreq) : m_sampleFreq(sampleFreq) {}

    void initPrn(int prn) {
      std::vector<std::complex<double>> codeVec(N);
      G2INIT sv(prn, 0); // Initialize with code phase 0 for FFT
      for (size_t idx = 0; idx < N; idx++) {
        int chipIdx = (idx / 16) % 1023;
        double codeVal = (double)sv.CODE[chipIdx];
        codeVec[idx] = std::complex<double>(codeVal, 0.0);
      }

      Fft::transform(codeVec, false); // Forward FFT of the code vector
      for (auto &val : codeVec) val = std::conj(val); // Conjugate for correlation
      codeFfts[prn] = codeVec; // Store the pre-computed FFT
  }
  AcqResult search(int prn, const std::vector<std::complex<double>> &rawData,
     double centerFreq, int binRange, float binWidth  ) {
    AcqResult bestResult = {0, 0, -99.0, 0.0};
    NCO nco(10, m_sampleFreq);
    std::vector<std::complex<double>> workspace(N);

    for(int bin = -binRange; bin<=binRange; bin++) {
      workspace = rawData; // Copy original data to workspace for processing
      nco.SetFrequency(centerFreq + (bin * binWidth));
       // 1. Carrier Wipeoff
       for (size_t idx = 0; idx < N; idx++) {
         uint32_t ncoIdx = nco.clk();
         std::complex<double> lo(nco.cosine(ncoIdx), nco.sine(ncoIdx));
         workspace[idx] *= lo; // Mix down to baseband
       }

       // 2. Correlation in Frequency Domain
    Fft::transform(workspace, false); // Forward FFT of the mixed signal
    for (size_t idx = 0; idx < N; idx++) {
        workspace[idx] *= codeFfts[prn][idx]; // Element-wise multiply with pre-computed code FFT 
    }
    Fft::transform(workspace, true); // Inverse FFT to get correlation in time domain

    //3. Peak Detection and SNR Calculation
    double maxMag = 0;
    int peakIndex = 0;
    double sumPower = 0;  
    for (size_t idx = 0; idx < N; idx++) {
        double mag = std::abs(workspace[idx]) / N; // Normalize by N
        sumPower += (mag * mag); // Accumulate power

        if (mag > maxMag) {
            maxMag = mag;
            peakIndex = idx;
        }
    }

   double noiseSum =0;
   int noiseCount = 0;
   for (size_t idx = 0; idx < N; idx++) {
        if (std::abs((int)idx - (int)peakIndex) > 16) {
          noiseSum += (std::abs(workspace[idx]) / N);
          noiseCount++;
        }
      }
    double cleanNoiseFloor = noiseSum / noiseCount; // Average noise magnitude
    double betterSNR = 20.0 * log10(maxMag / cleanNoiseFloor); // SNR in dB using magnitude ratio

    double peakPower = maxMag * maxMag;
    double avgNoisePower = (sumPower - peakPower) / (N - 1); // Exclude peak from noise power calculation
    double snr = 10.0 * log10(peakPower / avgNoisePower);
    if (snr > bestResult.snr) {
        //bestResult = {bin, peakIndex, maxMag, betterSNR}; // Update result if this bin has better SNR 
        bestResult = {bin, peakIndex, maxMag, snr}; // Update result if this bin has better SNR 
      }
    } 
    return bestResult;
  }
};

int main(int argc, char* argv[]) {
    // 1. Defaults
    float SampleFreq = 16.368e6;
    float RefFreq = 4.092e6;
    uint16_t binWidth = 500;
    int searchRange = 20;
    const char* filename = "IF.bin"; // Default file
    std::vector<int> prnsToSearch;

    // 2. Parse Arguments
    // Usage: ./program [filename] [prn1 prn2 ...]
    if (argc >= 2) {
        filename = argv[1];
    }

    if (argc >= 3) {
        // User provided specific PRNs
        for (int i = 2; i < argc; ++i) {
            prnsToSearch.push_back(std::stoi(argv[i]));
        }
    } else {
        // Default: Search all 32 PRNs
        for (int i = 1; i <= 32; ++i) prnsToSearch.push_back(i);
    }

    // 3. Initialize engine and data
    AcquisitionEngine engine(SampleFreq);
    std::vector<std::complex<double>> data(16384);
    
    FILE * IN = fopen(filename, "rb");
    if (!IN) {
        fprintf(stderr, "Error: Could not open file %s\n", filename);
        return 1;
    }
    stuffVector(data, IN);
    fclose(IN);

    // 4. Search Loop
    printf("Searching file: %s\n", filename);
    printf("------------------------------------------------------------------\n");
    printf("PRN | Bin | Freq Offset | Peak Index | Chip Phase | SNR (dB)\n");
    printf("------------------------------------------------------------------\n");

    for (int prn : prnsToSearch) {
        engine.initPrn(prn);
        AcqResult result = engine.search(prn, data, RefFreq, searchRange, binWidth);
        printf("%3d | %4d | %10.0f | %10d | %10.2f | %6.2f\n", 
               prn, result.bin, result.bin * (float)binWidth, 
               result.peakIndex, result.peakIndex / 16.0, result.snr);
    }

    return 0;
}
