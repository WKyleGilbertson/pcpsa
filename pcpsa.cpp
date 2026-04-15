#include <cstdio>
#include <iostream>
#include <vector>
#include <string>
#include "AcquisitionEngine.hpp"

// Updated to handle multiple blocks and your specific 16368-sample read
void stuffVector(std::vector<std::complex<double>> &vec, FILE *fp, int numMs)
{
  int8_t buffer[32736]; // 16368 samples * 2 bytes
  for (int ms = 0; ms < numMs; ms++)
  {
    size_t bytesRead = fread(buffer, 1, sizeof(buffer), fp);
    if (bytesRead < 32736)
      break;

    size_t offset = ms * 16384;
    for (size_t idx = 0; idx < 16368; idx++)
    {
      vec[offset + idx] = std::complex<double>((double)buffer[2 * idx], (double)buffer[2 * idx + 1]);
    }
    // Zero-pad the remaining 16 samples of this millisecond block
    for (size_t idx = 16368; idx < 16384; idx++)
    {
      vec[offset + idx] = std::complex<double>(0, 0);
    }
  }
}

int main(int argc, char *argv[])
{
  std::string filename = "IF.bin";
  std::vector<int> prnsToSearch;
  int numMs = 1; // Default to 1ms

  // Manual Argument Parsing
  for (int i = 1; i < argc; ++i)
  {
    std::string arg = argv[i];
    if (arg == "-f" && i + 1 < argc)
      filename = argv[++i];
    else if (arg == "-m" && i + 1 < argc)
      numMs = std::stoi(argv[++i]);
    //        else if (arg == "-p" && i + 1 < argc) prnsToSearch.push_back(std::stoi(argv[++i]));
    else if (arg == "-p" && i + 1 < argc) {
    std::string list = argv[++i];
    size_t start = 0;
    size_t end = list.find(',');

    while (end != std::string::npos) {
        std::string token = list.substr(start, end - start);
        if (!token.empty()) prnsToSearch.push_back(std::stoi(token));
        start = end + 1;
        end = list.find(',', start);
    }
    
    // This part captures the last item (or the only item if no comma exists)
    std::string lastToken = list.substr(start);
    if (!lastToken.empty()) {
        prnsToSearch.push_back(std::stoi(lastToken));
    }
   }
  }

  if (prnsToSearch.empty())
  {
    for (int i = 1; i <= 32; ++i)
      prnsToSearch.push_back(i);
    prnsToSearch.push_back(131);
    prnsToSearch.push_back(133);
    prnsToSearch.push_back(135);
  }

  AcquisitionEngine engine(16.368e6);
  std::vector<std::complex<double>> data(16384 * numMs);

  FILE *IN = fopen(filename.c_str(), "rb");
  if (!IN)
  {
    fprintf(stderr, "Error opening %s\n", filename.c_str());
    return 1;
  }
  stuffVector(data, IN, numMs);
  fclose(IN);

  printf("Searching %d ms of data from: %s\n", numMs, filename.c_str());
  printf("------------------------------------------------------------------\n");
  printf("PRN | Bin | Freq Offset | Peak Index | Chip Phase | SNR (dB)\n");
  printf("------------------------------------------------------------------\n");

  for (int prn : prnsToSearch)
  {
    engine.initPrn(prn);
    AcqResult result = engine.search(prn, data, 4.092e6, 20, 500);
    printf("%3d | %4d | %10.0f | %10d | %10.2f | %6.2f\n",
           prn, result.bin, result.bin * 500.0, result.peakIndex,
           result.peakIndex / 16.0, result.snr);
  }

  return 0;
}
