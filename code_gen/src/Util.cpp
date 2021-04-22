//
// Created by kazem on 4/23/17.
//

#include <atomic>
#include <fstream>
#include <iostream>
#include <sstream>
#include "Util.h"

namespace Sympiler {
namespace Internal {

namespace {
// We use 64K of memory to store unique counters for the purpose of
// making names unique. Using less memory increases the likelihood of
// hash collisions. This wouldn't break anything, but makes stmts
// slightly confusing to read because names that are actually unique
// will get suffixes that falsely hint that they are not.

    const int num_unique_name_counters = (1 << 14);
    std::atomic<int> unique_name_counters[num_unique_name_counters];

    int unique_count(size_t h) {
        h = h & (num_unique_name_counters - 1);
        return unique_name_counters[h]++;
    }

    void reset_counter(size_t h) {
        h = h & (num_unique_name_counters - 1);
        unique_name_counters[h] = 0;
    }
}
// There are three possible families of names returned by the methods below:
// 1) char pattern: (char that isn't '$') + number (e.g. v234)
// 2) string pattern: (string without '$') + '$' + number (e.g. fr#nk82$42)
// 3) a string that does not match the patterns above
// There are no collisions within each family, due to the unique_count
// done above, and there can be no collisions across families by
// construction.

void reset_counter(char prefix) {
    if (prefix == '$') prefix = '_';
    reset_counter((size_t) (prefix));
}

std::string unique_name(char prefix) {
    if (prefix == '$') prefix = '_';
    return prefix + std::to_string(unique_count((size_t) (prefix)));
}

std::string unique_name(const std::string &prefix) {
    std::string sanitized = prefix;

    // Does the input string look like something returned from unique_name(char)?
    bool matches_char_pattern = true;

    // Does the input string look like something returned from unique_name(string)?
    bool matches_string_pattern = true;

    // Rewrite '$' to '_'. This is a many-to-one mapping, but that's
    // OK, we're about to hash anyway. It just means that some names
    // will share the same counter.
    int num_dollars = 0;
    for (size_t i = 0; i < sanitized.size(); i++) {
        if (sanitized[i] == '$') {
            num_dollars++;
            sanitized[i] = '_';
        }
        if (i > 0 && !isdigit(sanitized[i])) {
            // Found a non-digit after the first char
            matches_char_pattern = false;
            if (num_dollars) {
                // Found a non-digit after a '$'
                matches_string_pattern = false;
            }
        }
    }
    matches_string_pattern &= num_dollars == 1;
    matches_char_pattern &= prefix.size() > 1;

    // Then add a suffix that's globally unique relative to the hash
    // of the sanitized name.
    int count = unique_count(std::hash<std::string>()(sanitized));
    if (count == 0) {
        // We can return the name as-is if there's no risk of it
        // looking like something unique_name has ever returned in the
        // past or will ever return in the future.
        if (!matches_char_pattern && !matches_string_pattern) {
            return prefix;
        }
    }

    return sanitized + "$" + std::to_string(count);
}

bool readCSCMatrixPattern(std::string fName, int &n, int& m, int &NNZ, int* &col,
                int* &row){
    /*This function reads the input matrix from "fName" file and
     * allocate memory for matrix A, L and U.
     * - The input file is a coordinate version and e
     * ach row of the file shows (col, row, nnz)
     * - The matrices are zero-indexed
     */

    std::ifstream inFile;
    inFile.open(fName);
    std::string line,banner, mtx, crd, arith, sym;
    /*  File format:
     *    %%MatrixMarket matrix coordinate real general/symmetric/...
     *    % ...
     *    % (optional comments)
     *    % ...
     *    #rows    #non-zero
     *    Triplet in the rest of lines: row    col    value
     */
    std::getline(inFile,line);
    for (unsigned i=0; i<line.length(); line[i]=tolower(line[i]),i++);
    std::istringstream iss(line);
    if (!(iss >> banner >> mtx >> crd >> arith >> sym)){
        std::cout<<"Invalid header (first line does not contain 5 tokens)\n";
        return false;
    }

    if(banner.compare("%%matrixmarket")) {
        std::cout<<"Invalid header (first token is not \"%%%%MatrixMarket\")\n";
        return false;
    }
    if(mtx.compare("matrix")) {
        std::cout<<"Not a matrix; this driver cannot handle that.\"\n";
        return false;
    }
    if(crd.compare("coordinate")) {
        std::cout<<"Not in coordinate format; this driver cannot handle that.\"\n";
        return false;
    }
    if(arith.compare("real") and arith.compare("integer")) {
        if(!arith.compare("complex")) {
            std::cout<<"Complex matrix; use zreadMM instead!\n";
            return false;
        }
        else if(!arith.compare("pattern")) {
            std::cout<<"Pattern matrix; values are needed!\n";
            return false;
        }
        else {
            std::cout<<"Unknown arithmetic\n";
            return false;
        }
    }
    while (!line.compare(0,1,"%"))
    {
        std::getline(inFile, line);
    }
    std::istringstream issDim(line);
    if (!(issDim >> n >> n >> NNZ)){
        std::cout<<"The matrix dimension is missing\n";
        return false;
    }
    if(n <= 0 || NNZ <= 0)
        return false;
    col = new int[n + 1]();
    // colL = new int[n + 1]; colU = new int[n + 1];
    row = new int[NNZ];
    // rowL = new int[factorSize]; rowU = new int[factorSize];
    // valL = new double[factorSize]; valU = new double[factorSize];
    if( !col || !row)
        return false;
    //Initializing the result vector
    int y, x, colCnt=0, nnzCnt=0;
    double value;

    col[0]=0;
    for (int i = 0; nnzCnt<NNZ; ) {//Reading from file row by row
        inFile>>x;x--;
        inFile>>y;y--;//zero indexing
        inFile>>value;
        if(y > n)
            return false;
        if(y==i){
            row[nnzCnt]=x;
            colCnt++; nnzCnt++;
        }
        else{//New col
            col[i+1]=col[i]+colCnt;
            i++;//next iteration
            colCnt=1;
            row[nnzCnt]=x;
            nnzCnt++;
        }

    }
    col[n]= col[n - 1] + colCnt;//last col

    return true;
}

}
}