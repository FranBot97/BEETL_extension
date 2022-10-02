/**
 ** Copyright (c) 2011-2014 Illumina, Inc.
 **
 ** This file is part of the BEETL software package,
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** Citation: Markus J. Bauer, Anthony J. Cox and Giovanna Rosone
 ** Lightweight BWT Construction for Very Large String Collections.
 ** Proceedings of CPM 2011, pp.219-231
 **
 **/

#ifndef TRANPOSEFASTA_INCLUDED
#define TRANPOSEFASTA_INCLUDED

#include "SequenceExtractor.hh"
#include "Tools.hh"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using std::string;
using std::vector;
typedef unsigned char uchar;

#define BUFFERSIZE 1024// 2^20
#define DIFFERENT_LEN 1 //if 0 all sequences must be of the same length
#define TERMINATE_CHAR '$'
#define DUMMY_CHAR '#'
#define SIZE_ALPHA 256
#define dataTypedimAlpha uchar
#define dataTypelenSeq uchar
#define output_array_file "position_array"
#define TEMP_DIR "BEETL-Temp"
//#define QS 0
//#define CYCLENUM 100


class SeqReaderFile;

class TransposeFasta
{
public:
    TransposeFasta();
    void init( SeqReaderFile *pReader, const string &input, const bool processQualities = true );
    ~TransposeFasta();

    bool convert( /*const string &input,*/ const string &output, bool generatedFilesAreTemporary = true );   //Input from Fasta file (converts Fasta File into cyc Files)
    bool convertLenOrder(const string &input, const string &output, bool generatedFilesAreTemporary = true ); //Used if sequences have different length
    bool convertRLO(const string &input, const string &output, bool generatedFilesAreTemporary = true);
    bool convertAcceptDiffLen( const string &input, const string &output, bool generatedFilesAreTemporary = true);
    bool inputCycFile( const string &cycPrefix );                                    //Input from cyc files
    bool convertFromCycFileToFastaOrFastq( const string &fileInputPrefix, const string &fileOutput, bool generatedFilesAreTemporary = true, SequenceExtractor *sequenceExtractor = NULL );      //Convert cyc files into Fasta or Fastq File
    bool hasProcessedQualities() const
    {
        return processQualities_;
    }

    int findInfoSeq( const string &input, const string &output, vector<int> &keepOrder);

    SequenceLength lengthRead;    //Length of each text
    LetterNumber lengthTexts;   //Total length of all texts without $-symbols

    SequenceNumber nSeq;   //number total of texts in filename1
    LetterNumber freq[256];  //contains the distribution of the symbols. It is useful only for testing. It depends on the #characters

private:
    SeqReaderFile *pReader_;
    uint cycleNum_;
    vector<FILE *> outputFiles_;
    vector<vector<uchar> > buf_;
    //    FILE* outputFiles_[CYCLENUM];

    //    uchar buf_[CYCLENUM][BUFFERSIZE];
    bool processQualities_;
};

#endif