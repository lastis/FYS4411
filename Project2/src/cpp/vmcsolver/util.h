#include <cmath>
#include "../CPhys/CPhys.h"

namespace util {
    void blockingVar(int nBins, Vector& dataVec, 
            Vector& contStdArray, Vector& contBlockSizes){
        // Get the data from the data vector.
        double* data = dataVec.getArrayPointer();
        int dataLength = dataVec.getLength();
        // Make a list of block sizes. We
        contBlockSizes = Vector(nBins);
        contBlockSizes.linspace(1,dataLength);
        double* blockSizes = contBlockSizes.getArrayPointer();
        // Chose a block of a spesific block size, 
        // and calculate the sample with this. 
        contStdArray = Vector(nBins);
        double* stdArray = contStdArray.getArrayPointer();
        for (int bin = 0; bin < nBins; bin++) {
            int blockSize = blockSizes[bin];
            int nBlocks = dataLength/blockSize;
            int cnt = 0;
            double std = 0;
            // Calculate the mean with a given blocksize.
            for (int block = 0; block < nBlocks; block++) {
                double sample;
                double samples = 0;
                double samplesSq = 0;
                for (int i = 0; i < blockSize; i++) {
                   sample = data[cnt];
                   samples += sample; 
                   samplesSq += sample*sample;
                   cnt++;
                }
                std += samplesSq/blockSize - samples*samples/(blockSize*blockSize);
            }
            stdArray[bin] = std/nBlocks;
        }
    }
}
