#include <cmath>
#include "../CPhys/CPhys.h"
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

namespace util {
    Vector getMeanArray(int blockSize, Vector& dataVec){
        double* data = dataVec.getArrayPointer();
        int nBlocks = dataVec.getLength()/blockSize;
        Vector meanArray = Vector(nBlocks);
        double* means = meanArray.getArrayPointer();
        int cnt = 0;
        // Calculate the mean with a given blocksize.
        for (int block = 0; block < nBlocks; block++) {
            double sample;
            double mean = 0;
            for (int i = 0; i < blockSize; i++) {
               sample = data[cnt];
               mean += sample; 
               cnt++;
            }
            means[block] = mean/blockSize;

        }
        return meanArray;
    }

    /** \brief Normalize a vector
     *
     * \param vec Vector
     */         
    void normalize(Vector &vec){
        int N = vec.getLength();
        double* pvec = vec.getArrayPointer();
        double sum = 0;
        for (int i = 0; i < N; i++) 
        {
            sum += pvec[i];
        }
        for (int i = 0; i < N; i++) 
        {
            pvec[i] /= sum;
        }
    }

    /** \brief Load matrix from txt file. Warning, the matrix dimensions
     * must be correct.
     *
     * \param adress Adress to file. Warning, from the location of this
     * function.
     * \param fileName Name of file. 
     * \param mat Matrix to load into. 
     */         
    void loadFromFile(string adress, string fileName, Matrix& mat){
        int N = mat.getN();
        int M = mat.getM();
        double** pMat = mat.getArrayPointer();
        fstream myFile ((adress + fileName).c_str(), ios::in);
        double number;
        for (int i = 0; i < N; i++) 
        {
            for (int j = 0; j < M; j++) 
            {
                myFile >> number;
                pMat[i][j] = number;
            }
        }
        myFile.close();
    }

    /** \brief Write vector to file. Will overwrite.
     *
     * \param adress Adress to file. Warning, from the location of this
     * function.
     * \param fileName Name of file. 
     * \param mat Matrix to print. 
     */         
    void writeToFile(string adress, string fileName, Matrix mat){
        ofstream myFile;
        int N = mat.getN();
        int M = mat.getM();
        double** pMat = mat.getArrayPointer();
        // Dump results to the end of the file. 
        myFile.open((adress + fileName).c_str());
        for (int i = 0; i < N; i++) 
        {
            for (int j = 0; j < M; j++) 
            {
                myFile << pMat[i][j] << " ";
            }
            myFile << endl;
        }
        myFile << endl;
        myFile.close();
    }

    /** \brief Write vector to file. Will overwrite.
     *
     * \param adress Adress to file. Warning, from the location of this
     * function.
     * \param fileName Name of file. 
     * \param array Vector to print. 
     */         
    void writeToFile(string adress, string fileName, Vector array){
        ofstream myFile;
        double* pArray = array.getArrayPointer();
        // Dump results to the end of the file. 
        myFile.open((adress + fileName).c_str());
        for (int i = 0; i < array.getLength(); i++) {
            myFile << pArray[i] << " ";
        }
        myFile << endl;
        myFile.close();
    }

    /** \brief Appends vector to file.
     *
     * \param adress Adress to file. Warning, from the location of this
     * function.
     * \param fileName Name of file. 
     * \param array Vector to print. 
     */         
    void appendToFile(string adress, string fileName, Vector array){
        ofstream myFile;
        double* pArray = array.getArrayPointer();
        // Dump results to the end of the file. 
        myFile.open((adress + fileName).c_str(),std::ios_base::app);
        for (int i = 0; i < array.getLength(); i++) {
            myFile << pArray[i] << " ";
        }
        myFile << endl;
        myFile.close();
    }


    void blockingVar(int minBlockSize, Vector& dataVec, 
            Vector& contStdArray, Vector& contBlockSizes){
        // Get the data from the data vector.
        double* data = dataVec.getArrayPointer();
        int dataLength = dataVec.getLength();
        // From the minimum block size we can find out what the
        // maxium number of block "bins" there are. 
        int nBins = dataLength/minBlockSize;
        // Make a list of block sizes. 
        contBlockSizes = Vector(nBins);
        contBlockSizes.linspace(1,dataLength);
        double* blockSizes = contBlockSizes.getArrayPointer();
        // Make an array we can fill with variances. 
        contStdArray = Vector(nBins);
        double* stdArray = contStdArray.getArrayPointer();
        // This counter is used to get rid of unecessary data points.
        int zeroes = 0;
        // Chose a block of a spesific block size, 
        // and calculate the sample with this. 
        for (int bin = 0; bin < nBins; bin++) {
            int blockSize = blockSizes[bin];
            int nBlocks = dataLength/blockSize;
            // If we are not on the last "bin", we want to check
            // if the next bin has the same number of blocks that we are currently
            // on. If it's the same, we skip it. This ensures the maximum
            // number of bins possible.
            if (bin != nBins-1){
                int blockSizeNext = blockSizes[bin+1];
                int nBlocksNext = dataLength/blockSizeNext;
                if (nBlocks == nBlocksNext){
                    zeroes++;
                    continue;
                }
            }
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
                samples = samples/blockSize;
                samplesSq = samplesSq/blockSize;
                std += samplesSq - samples*samples;
            }
            stdArray[bin] = sqrt(std/nBlocks);
        }
        // This method of finding the maximum number of 
        // block sizes has left us with a lot of zeroes. 
        // This next code reduces that. (Added -1 at the end of
        // the array length because the last element always became zero).
        Vector newStdArray = Vector(nBins-zeroes-1);
        Vector newBlockSizes = Vector(nBins-zeroes-1);
        double* ptr1 = newStdArray.getArrayPointer();
        double* ptr2 = newBlockSizes.getArrayPointer();
        int cnt = 0;
        for (int bin = 0; bin < nBins; bin++) {
            if(stdArray[bin] == 0) continue;
            ptr1[cnt] = stdArray[bin];
            ptr2[cnt] = blockSizes[bin];
            cnt++;
        }
        contStdArray = newStdArray;
        contBlockSizes = newBlockSizes;
    }
}
