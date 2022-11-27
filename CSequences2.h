/***************************************************************************************************
*  The PolyMoSim project is distributed under the following license:
*  
*  Copyright (c) 2006-2022, Christoph Mayer, Forschungsmuseum Alexander Koenig, Bonn, Germany
*  All rights reserved.
*  
*  Redistribution and use in source and binary forms, with or without
*  modification, are permitted provided that the following conditions are met:
*  1. Redistributions of source code (complete or in parts) must retain
*     the above copyright notice, this list of conditions and the following disclaimer.
*  2. Redistributions in binary form must reproduce the above copyright
*     notice, this list of conditions and the following disclaimer in the
*     documentation and/or other materials provided with the distribution.
*  3. All advertising materials mentioning features or any use of this software
*     e.g. in publications must display the following acknowledgement:
*     This product includes software developed by Christoph Mayer, Forschungsmuseum
*     Alexander Koenig, Bonn, Germany.
*  4. Neither the name of the organization nor the
*     names of its contributors may be used to endorse or promote products
*     derived from this software without specific prior written permission.
*  
*  THIS SOFTWARE IS PROVIDED BY CHRISTOPH MAYER ''AS IS'' AND ANY
*  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
*  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHTHOLDER OR ITS ORGANISATION BE LIABLE FOR ANY
*  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
*  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
*  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
*  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
*  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
*  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*  
*  IMPORTANT (needs to be included, if code is redistributed):
*  Please not that this license is not compatible with the GNU Public License (GPL)
*  due to paragraph 3 in the copyright. It is not allowed under any
*  circumstances to use the code of this software in projects distributed under the GPL.
*  Furthermore, it is not allowed to redistribute the code in projects which are
*  distributed under a license which is incompatible with one of the 4 paragraphs above.
*  
*  This project makes use of code coming from other projects. What follows is a complete
*  list of files which make use of external code. Please refer to the copyright within
*  these files.
*  
*  Files in tclap foler:         Copyright (c) 2003 Michael E. Smoot
*                                See copyright in tclap/COPYRIGHT file for details.	
*  discrete_gamma.c:             Copyright 1993-2004 by Ziheng Yang.
*                                See copyright in this file for details.
*  CRandom.h:                    Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura
*                                See copyright in this file for details.
***************************************************************************************************/

//----------------------------------------------------------------------------------------------
// Beginning of copyright statement
//
// Program:     For multiple programs by the author
// Version:     2.0 of this file
// Copyright:   Christoph Mayer, 2003-2011
// Address:     Christoph Mayer,
//              Forschungsmuseum Alexander Koenig
//              53123 Bonn
//              Germany
//
// The source code of SAMS is copyright protected and NOT released for public use.
// It is prohibited to use this source code for any purpose whatsoever without
// explicit permission from the author Christoph Mayer.
// If permission to use this source code has been given to individual persons,
// this copyright statement must be included in any new source code that uses source
// code from this file or project. In particular, it must be included in any modified
// version or copy of this source code.
// If permission to use this source code has been given to individual persons,
// it is still not allowed to distribute or publish any code of this file or project
// by these persons.
// Programs that use this code can only be distributed or made available to others in
// binary form and only if no copyright of third parities is infringed.
// Publications that present or describe methods that where implemented or programmed using
// source code or algorithmic ideas from this file or project must include the author
// Christoph Mayer of this source code as one of the authors of the publication.



//////////////////////////////////////////////////////////////////////
// CSequences.h: interface for the CSequences class.
//
//////////////////////////////////////////////////////////////////////

#ifndef  CSEQUENCES2_H
#define  CSEQUENCES2_H

#include <vector>
#include <map>
#include "faststring2.h"
#include "fast-realloc-vector.h"
#include "CSplit2.h"
#include "CSequence_Mol2_1.h"
#include <cassert>
#include <cmath>
#include <utility>
#include "basic-DNA-RNA-AA-routines.h"

//#include "basic-DNA-RNA-AA-routines.h"

#define PrintMessage_cerr(s1)                fputs(s1, stderr)
#define ErrorMessage(s1)                     fputs(s1, stderr)
#define flush_cerr()                         fflush(stderr)


//////////////////////////////////
// Global Functions
//////////////////////////////////

//template<class T>
inline void add_or_count(std::map<faststring, unsigned> &m, faststring &x)
{
  std::map<faststring,unsigned>::iterator it;

  it = m.find(x);
  if (it == m.end() )
  {
    m[x] = 1;
  }
  else
  {
    ++it->second;
  }
}


//////////////////////////////////
// Data types
//////////////////////////////////

//////////////////////////////////
// Class CSequences
//////////////////////////////////

class CSequences2
{
public:
  typedef char chartype;

private:
  unsigned                  taxaNum;
  unsigned	            posNum;
  char                      ambig_char;


  CSequence_Mol::DataTypesEnum	datatype;
  std::vector<unsigned>         originalPosNumbers;  // In case we excluded positions from the alignment, we want to remember the original positions.
  std::vector<CSequence_Mol*>   seqData;             // The vector of sequences


  bool equal_length_of_all_sequences()
  {
    size_t len;
    unsigned i;

    if (taxaNum==0)
      return true;

    len = seqData[0]->length();

    for (i=1; i<taxaNum; ++i)
    {
      if (len != seqData[i]->length())
	return false;
    }
    return true;
  }

public:
  CSequences2(CSequence_Mol::DataTypesEnum Dt):
  taxaNum(0), posNum(0), ambig_char('?'), datatype(Dt)
  {}

  ~CSequences2()
  {
    int i, n=taxaNum;

    for (i=0; i<n; ++i)
    {
      delete seqData[i];
    }

  }


   unsigned         GetTaxaNum() { return taxaNum;}
   unsigned         GetPosNum()  { return posNum;}

   CSequence_Mol::DataTypesEnum    get_datatype()
   {
     return datatype;
   }

   chartype         GetChar(unsigned TaxaIndex,
			    unsigned PosIndex) const
   {
     assert(TaxaIndex < taxaNum);
     assert(PosIndex  < posNum);
     return seqData[TaxaIndex]->get_pos(PosIndex);
   }

/*    void             SetChar(unsigned TaxaIndex, unsigned PosIndex, chartype mychar) */
/*    { */
/*      assert(TaxaIndex < taxaNum); */
/*      assert(PosIndex < posNum); */
/*      seqData[TaxaIndex]->set_pos(PosIndex, mychar); */
/*    } */

   void SetOriginalPosNumber(unsigned index, unsigned theOriginalPosNumber)
   {
     assert(index < posNum);
     originalPosNumbers[index] = theOriginalPosNumber;
   }

   unsigned GetOriginalPosNumber(unsigned index) const
   {
     assert(index < posNum);
     return originalPosNumbers[index];
   }

   CSequence_Mol *get_seq(unsigned i)
   {
     return seqData[i];
   }

   char get_ambiguity_character()
   {
     return ambig_char;
   }

   void unique_phylip_names(std::vector<faststring> &snames)
   {
     std::map<faststring, unsigned>  m_snames;
     std::map<faststring, unsigned>::iterator  f_it;
     unsigned                             digits = log10(taxaNum) + 1;
     unsigned                             i;
     
     for (i=0; i < taxaNum; ++i)
     {
       snames.push_back(seqData[i]->getPhylipName());
       f_it = m_snames.find(snames[i]);

       if (f_it != m_snames.end() ) // This name is not unique
       {
	 ++f_it->second;
       }
       else
       {
	 m_snames.insert(std::make_pair(snames[i], 1));
       }
     }
     // We have build the map - now we need to create unique names:
     f_it = m_snames.begin();
     while (f_it != m_snames.end() )
     {
       if (f_it->second > 1)
       {
	 unsigned internal_counter = 1;
	 unsigned j;
	 for (j=0; j<taxaNum; ++j)
	 {
	   if (snames[j] == f_it->first)
	   {  
	     snames[j]  = snames[j].substr(0,10-digits);

	     faststring nn = faststring(internal_counter, '0', digits);
	     nn.replace_char(' ', '0');
	     snames[j].append(nn);
	     ++internal_counter;
	   }
	 }
       }


       ++f_it;
     }

     

   } 


   unsigned PairwiseSequenceSymbolMatches(unsigned taxon1, unsigned taxon2,
					  char     numDistCharacters,
					  const signed char* pos_vec,
					  signed char* match_vec) const
   {
     unsigned  count    = 0;
     unsigned  pos;
     
     for (pos = 0; pos < posNum; ++pos)
     {
       if (pos_vec[pos]  &&
	   seqData[taxon1]->get_pos(pos) == seqData[taxon2]->get_pos(pos) &&
	   seqData[taxon1]->get_pos(pos) < numDistCharacters)
       {
	 ++count;
	 match_vec[pos] = 1;
       }
       else
	 match_vec[pos] = 0;
     }
     return count;
   }
   
   
   unsigned PairwiseSequenceSymbolMatches(
					  unsigned taxon1,
					  const faststring & ref_seq,
					  char               numDistCharacters,
					  const signed char* pos_vec,
					  signed char* match_vec) const
  

   // Ist hier alles OK?????????????????????????????????????????
   {
     unsigned  count    = 0;
     unsigned  pos;
     
     for (pos = 0; pos < posNum; ++pos)
     {
       if (pos_vec[pos]  &&  seqData[taxon1]->get_pos(pos) == ref_seq[pos]
	              &&  ref_seq[pos] < numDistCharacters)
       {
	 ++count;
	 match_vec[pos] = 1;
       }
       else
	 match_vec[pos] = 0;
     }
     return count;
   }
   

   void     ConsensusSequence(faststring&           conSeq,
			      const CSplit&         setOfTaxa,
			      char                  numDistCharacters,
			      unsigned              numSymbols,
			      double                consensusThreshold )
   {
     unsigned          pos;
     unsigned          taxon, maxindex, equalindex;
     unsigned          taxaCount = 0;
     unsigned          consensusSymMinimum;
     char              i;


     std::vector<unsigned>  counterSymbolsVec(numSymbols,0);
     
     for (pos=0; pos < posNum; ++pos)
     {
       // Initialise variable
       taxaCount = 0;
       for (i=0; i < numDistCharacters; ++i)
	 counterSymbolsVec[i] = 0;
       
       // Count number of occuring symbols for this position over all taxa
       for (taxon = 0; taxon < taxaNum; ++taxon)
       {
	 if (setOfTaxa.test(taxon))
	 {
	   ++taxaCount;
	   ++counterSymbolsVec[seqData[taxon]->get_pos(pos)];
	 }
       }
       
       consensusSymMinimum = (unsigned) std::ceil(consensusThreshold * taxaCount);
       
       maxindex = 0;
       equalindex = numDistCharacters;
       for (i = 1; i < numDistCharacters; ++i)
       {
	 if (counterSymbolsVec[i] >= counterSymbolsVec[maxindex])
	 {
	   if (counterSymbolsVec[i] == counterSymbolsVec[maxindex])
	     equalindex = i;
	   maxindex = i;
	 }
       }
       if (counterSymbolsVec[maxindex] >= consensusSymMinimum &&
	   maxindex != equalindex)
	 conSeq[pos] = maxindex;
       else // Default value, in case consensus cannot be determined
	 conSeq[pos] = numDistCharacters;
     }
   }
   

   void GetSymbolFrequenciesAtPosition(unsigned pos, const CSplit &taxaSet,
				       unsigned numSymbols,
				       std::vector<unsigned>& frequencies)
   {
     unsigned i;
     //  chartype sym;
     
     //  assert(frequencies.size() <= numDistCharacters);
     
     for (i=0; i<numSymbols; ++i)
       frequencies[i]=0;
     
     for (i=0; i<taxaNum; ++i)
     {
       if (taxaSet.test(i))
       {
	 ++frequencies[seqData[i]->get_pos(pos)];
       }
     }
   }

   int read_from_Fasta_File(CFile& infile,
			    CSequence_Mol::processing_flag pflag,
			    unsigned first_seq,
			    unsigned last_seq,
			    bool report_status_to_cerr=true
			    )
   {
     CSequence_Mol   *seq;
     unsigned        count_seq = 0;
     unsigned        count_seqs_stored = 0;


     // TODO
     (void) report_status_to_cerr;

     // Remove all sequences from this object - we want to read a new data set.
     clear(datatype);

     while (!infile.eof())
     {
       seq = new CSequence_Mol (datatype);
       ++count_seq;

       seq->readFastaSequence_generic(infile, pflag);

       // PrintMessage_cerr(seq.getName());PrintMessage_cerr("-\n");PrintMessage_cerr(seq.getFullName());PrintMessage_cerr("-\n");PrintMessage_cerr(seq.getDescription());PrintMessage_cerr("-\n");
       if ( infile.fail() && !infile.eof() )
       {
	 PrintMessage_cerr("\n\n");
	 faststring errormsg;
	 errormsg   =  "An error occurred while reading the input file. It might not be a valid fasta file.\n";
	 errormsg  +=  "File position: line ";
	 errormsg  +=  faststring(infile.line());
	 ErrorMessage(errormsg.c_str());
	 PrintMessage_cerr("\n");
	 flush_cerr();
	 return -25;
       }

       if ( count_seq < first_seq || count_seq > last_seq )
       {
	 if (false)
	 {
	   PrintMessage_cerr("Skipping sequence ");
	   PrintMessage_cerr(faststring(count_seq).c_str());
	   PrintMessage_cerr(": ");
	   PrintMessage_cerr(seq->getName());
	   PrintMessage_cerr("\n");
	   flush_cerr();
	 }
	 continue;
       }

       PrintMessage_cerr("Processing sequence ");
       PrintMessage_cerr(faststring(count_seq).c_str());
       PrintMessage_cerr(": ");
       PrintMessage_cerr(seq->getName());
       PrintMessage_cerr("\n");
       flush_cerr();
       
       seqData.push_back(seq);
       ++count_seqs_stored;
     } // End while

     // Check that all sequences have the same length - we require that they have the same length
     if (!equal_length_of_all_sequences() )
       return -1;

     // Check equal data type and ambig character:

     taxaNum = count_seqs_stored;

     if (taxaNum > 0)
     {
       posNum  = seqData[0]->length();
       ambig_char = seqData[0]->get_ambiguity_character();
     }

     return 0;
   }  // End this element function.


   int read_from_Phylip_File(CFile& infile)
   {
     CSequence_Mol   *seq;
     //     unsigned        count_seq;
     //     unsigned        global_sequence_from=0, global_sequence_to=-1u;
     //     bool            global_gaps_to_Ns=false;
     unsigned        all_lines_N;
     unsigned        curr_taxon;

     // The allowed symbols from the phylip manual:


     char ch;
     faststring line;
     std::vector<faststring> fsvec;
     faststring *p;

     clear(datatype); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     originalPosNumbers.clear();

     ch = infile.peek();

     // If the first char is a white space, we expect this
     // is the line with the number of taxa and seq. posiitons.
     if ( isspace(ch) )
     {
       infile.getline(line);
       split(fsvec, line);

       if (fsvec.size()==2 && fsvec[0].isAnUnsigned() && fsvec[1].isAnUnsigned() )
       {
	 taxaNum = fsvec[0].ToUnsigned();
	 posNum  = fsvec[1].ToUnsigned();

	 std::cout << "Read number of taxa     " << taxaNum << std::endl;
	 std::cout << "Read number of residues " << posNum  << std::endl;
       }
     }

     //     std::vector<faststring*> sqs;
     std::vector<faststring*> all_lines;

/*      if (taxaNum == 0 && posNum == 0) */
/*      { */
/*        faststring *pp = new faststring; */
/*        *pp = line; */
/*        all_lines.push_back(pp); */
/*      } */

     // We read the complete file into a vector of strings.
     p = new faststring();
     infile.getline(*p);
     p->removeSpacesFront();
     p->removeSpacesBack();

     while (!infile.fail() )
     {
       all_lines.push_back(p);
       p = new faststring();
       infile.getline(*p);
       p->removeSpacesFront();
       p->removeSpacesBack();
     }

     // The problem is: The phylip format is very difficult
     // to parse. Main problem is to distinguish between
     // sequential and interleaved format.
     // The interleaved format is distinguished as follows:
     // After each block there must be a blank line.
     // So if we find blank lines which are not the last line
     // the format is interleaved. Otherwise it is sequential.

     // For simplicity we check whether the last line is empty.
     // If it is, we remove it.
     
     unsigned i;
     all_lines_N = all_lines.size();

     //     std::cerr << all_lines[all_lines_N-1]->length() << std::endl;
          
     while (all_lines[all_lines_N-1]->length() == 0)
     {
       --all_lines_N;
       delete all_lines[all_lines_N];
       std::vector<faststring*>::iterator it = all_lines.end();
       --it;
       all_lines.erase(it);
     }

     unsigned number_blank_lines=0;

     for (i=0; i<all_lines_N; ++i)
     {
       if (all_lines[i]->length() == 0)
       {
	 ++number_blank_lines;
       }
     }
     
     if (all_lines_N == 0 || number_blank_lines == all_lines_N) // Nothing to do???
       return 0;
  
     if (number_blank_lines == 0) // sequential format
     {
       // This can only be interpreted, if the number
       // of taxa and positions has been specified in
       // the first line.
       // What makes it more complicated
       // is the fact, that the phylip format tolerates
       // numbers and other stuff at the beginning of pure sequence lines.

       faststring  seq_taxon_name;
       faststring  seq_in_one_line;
       //       unsigned    current_length = 0;
       unsigned    curr_line_num=0;

       faststring  *currLine;
       unsigned    residuse_added_this_line;

       if (taxaNum == 0 || posNum == 0)
       {
	 // delete all faststrings in all_lines
	 for (unsigned ind=0; ind < all_lines.size(); ++ind)
	   delete all_lines[ind];
	 return -3;
       }


       // For all taxa
       //
       for (curr_taxon = 0; curr_taxon < taxaNum; ++curr_taxon)
       {
	 residuse_added_this_line = 0;
	 if (curr_line_num >= all_lines_N)
	 {
	   std::cerr << "Parse error while reading the input file:\nUnexpected end of file in line: " << curr_line_num <<  std::endl;
	   std::cerr << "Trying to read " << taxaNum << " taxa but found only " << curr_taxon << "!" <<  std::endl;

	   exit(-45);
	 }

	 currLine       = all_lines[curr_line_num];
	 const char* pos = currLine->c_str();
	 //	 currLineLen     = currLine->length();
	 seq_taxon_name.clear();
	 seq_in_one_line.clear();
	 
	 // Copy the sequence name
	 for (i=0; i<10 && *pos != '\0'; ++i, ++pos)
	 {
	   seq_taxon_name.push_back(*pos);
	 }
	 if (i<10)
	   seq_taxon_name.append(10-i, ' ');
	 if (*pos) // we neeed to read the rest of the line
	 {
	   // Copy the sequence data:
	   {
	     // Find first character that shall be copied.
	     // The value of i should be 10, the position after the sequence name
	     char c;

	     c = currLine->operator[](i);
	     while (!is_phylip_aa_dna_symbol(c) && !is_allowed_non_ABC(c) )
	     {
	       ++i;
	       c = currLine->operator[](i);
	     }
	   }
	   //	   i = currLine->find_first_of(DNARNA_allowed_symbols, i);

	   //	   seq_in_one_line.assign(*currLine, i);
	   seq_in_one_line.assign(*currLine, (size_t)i, (size_t)std::string::npos); // Third argument for compatibility with std::string
	   residuse_added_this_line = currLine->length()-i;
	 }
	 // Read lines until we have the required number of residues:
	 while (residuse_added_this_line < posNum)
	 {
	   ++curr_line_num;
	   if (curr_line_num >= all_lines_N)
	   {
	     std::cerr << "Parse error while reading the input file. The file has been interpreted as being in sequential phylip format." << std::endl;
	     std::cerr << "Several problems can trigger this error: (i) Wrong number of residues in sequences compared to the number specified " << std::endl;
	     std::cerr << "in the file header. (ii) Wrong number of taxa (sequences) compared to the number specified in the file header." << std::endl;
	     std::cerr << "(iii) Missing blank line between blocks of data if this file should be in interleaved format." << std::endl;
	     exit(-44);
	   }
	   currLine       = all_lines[curr_line_num];
	   //	   const char* pos = currLine->c_str();
	   {
	     // Find first character that shall be copied.
	     char c;
	     i=0;
	     c = currLine->operator[](i);
	     while (!is_phylip_aa_dna_symbol(c) && !is_allowed_non_ABC(c) )
	     {
	       ++i;
	       c = currLine->operator[](i);
	     }
	   }
	   //	   i = currLine->find_first_of(DNARNA_allowed_symbols, i);


	   //	   seq_in_one_line.append(*currLine, i);
	   seq_in_one_line.append(*currLine, i, std::string::npos);	   // Third argument for compatibility with std::string
	   residuse_added_this_line += currLine->length()-i;
	 }
	 // If we are here, the sequence should have been copied to
	 // seq_in_one_line

	 seq = new CSequence_Mol (datatype);
	 seq->set_taxon_and_sequence(seq_taxon_name, seq_in_one_line);
	 seqData.push_back(seq);

	 // The next line should start with a new sequence:
	 ++curr_line_num;

       } // End for loop - for all taxa

       if (curr_line_num != all_lines_N)
       {
	 std::cerr << "Parse error while reading the input file: This can have several reasons:\n"
                      "(i) The number of taxa found in the first block does not match the number specified in the file header.\n"
	              "(ii) Since only one data block was found, this file is interpreted as sequential format. If this should be a file\n"
                      "in interleaved format, blank line must be added between blocks of data.\n"
                      "(iii) It could also be that the number of residues in the sequences is larger than specified in the header, so that addtional\n"
	              "lines of data are interpreted as additional taxa in the sequential format." << std::endl;
	 exit(-67);
       }

     }    // End sequential format
     else // Interleaved format
     {
       // In the first block, the first 10 characters are
       // the sequence name. Then the first base starts
       // the sequence.
       // After the first block, the sequence beginns/continues at
       // the first base symbol in each line.
       unsigned     curr_line_num=0;
       faststring   *currLine;

       //       unsigned current_length = 0;
       unsigned curr_taxon = 0;

       faststring *pseq_taxon_name;
       faststring *pseq_in_one_line;


       std::vector<faststring*> taxonnames;
       std::vector<faststring*> sequences;

       currLine        = all_lines[curr_line_num];
       const char* pos = currLine->c_str();

       // read first block - until we find a blank line:
       while (currLine->length() != 0) // We read all lines of the first block. Counter: curr_line_num 
       {
	 pseq_taxon_name  = new faststring;
	 pseq_in_one_line = new faststring;

	 // Copy the sequence name
	 for (i=0; i<10 && *pos != '\0'; ++i, ++pos)
	 {
	   pseq_taxon_name->push_back(*pos);
	 }
	 if (i<10)
	   pseq_taxon_name->append(10-i, ' ');

	 //************
	 if (*pos) // we neeed to read the rest of the line
	 {
	   // Copy the sequence data:
	   {
	     // Find first character that shall be copied.
	     // The value of i should be 10, the position after the sequence name
	     char c;
	     
	     c = currLine->operator[](i);
	     while (!is_phylip_aa_dna_symbol(c) && !is_allowed_non_ABC(c) )
	     {
	       ++i;
	       c = currLine->operator[](i);
	     }
	   }
	   //	   i = currLine->find_first_of(DNARNA_allowed_symbols, i);
	   //	   pseq_in_one_line->assign(*currLine, i);
	   pseq_in_one_line->assign(*currLine, i, std::string::npos); 	   // Third argument for compatibility with std::string
	 } // if (*pos) // we neeed to read the rest of the line
	 
	 // Another taxon has been read: Let us store it:

	 taxonnames.push_back(pseq_taxon_name);
	 sequences.push_back(pseq_in_one_line);	 

	 // Read next line with next taxon name
	 ++curr_taxon;
	 ++curr_line_num;
	 if (curr_line_num >= all_lines_N)
	   break;
	 // If we break, currLine will be set to its new value further down
	 currLine       = all_lines[curr_line_num];
       } // End while loop over all non blank lines in first block

       if (taxaNum == 0 && posNum == 0)
	 taxaNum = curr_taxon;

       // Check that the number of taxa that might had been specified agrees with the
       // number of taxa we found now.
       if (taxaNum != 0 && taxaNum != curr_taxon) // do we have a problem?
       {
	 std::cerr << "Parse error while reading the input file:\nThe number of taxa in the first block of the input file ("
		   << curr_taxon  << ") does not match the number of taxa specified in the header of the phylip file (" << taxaNum << ") " 
		   << std::endl;
	 exit(-55);
       }

       // This line should be blank: Skip blank lines between the blocks
       while (curr_line_num < all_lines_N && all_lines[curr_line_num]->length()==0)
       {
	 ++curr_line_num;
       }

       // Now we need to read all remaining blocks:
       while (true)  // (curr_line_num < all_lines_N)
       {
	 // Read the next block:
	 for (curr_taxon=0; curr_taxon < taxaNum; ++curr_taxon)
	 {
	   // The current line should be the first line of this block
	   if (curr_line_num >= all_lines_N)
	   {
	     // We are within a block and ran out of lines:
	     std::cerr << "Parse error while reading the input file:\nUnexpected end of file. More data has been expected." << std::endl;
	     exit(-46);
	   }
	   
	   currLine       = all_lines[curr_line_num];
	   {
	     // Find first character that shall be copied.
	     char c;
	     i=0;
	     c = currLine->operator[](i);
	     while (!is_phylip_aa_dna_symbol(c) && !is_allowed_non_ABC(c) )
	     {
	       ++i;
	       c = currLine->operator[](i);
	     }
	   }
	   // Copy the sequence on this line to our data.
	   //	   sequences[curr_taxon]->append(*currLine, i);
	   sequences[curr_taxon]->append(*currLine, i, std::string::npos); 	   // Third argument for compatibility with std::string
	   ++curr_line_num;
	 } // End of for loop which reads the next block
	 
	 // The next line(s) should be blank or we should be at the end of the list of lines
	 if (curr_line_num < all_lines_N && all_lines[curr_line_num]->length()!=0)
	 {
	   std::cerr << "Parse error while reading the input file:\nEither the number of lines of data in one of the data blocks is not identical to the number of taxa\n"
	                "or a blank line might be missing between different data blocks in the interleaved format." << std::endl;
	   std::cerr << "This problem was detected in line " << curr_line_num+1 << " of the file, but it could have occured earlier." << std::endl;
	   exit(-77);
	 }

	 while (curr_line_num < all_lines_N && all_lines[curr_line_num]->length()==0)
	 {
	   ++curr_line_num;
	 }
	 if (curr_line_num >= all_lines_N) // we parsed all lines of the data set.
	 {
	   break;
	 }
	 // This should be the first line of the next block
	 currLine = all_lines[curr_line_num];
	 // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
       } // End while loop read all blocks

       if (taxaNum > 0 && posNum == 0)
	 posNum = sequences[0]->length();

       // Now we should have read the sequences:
       for (curr_taxon=0; curr_taxon < taxaNum; ++curr_taxon)
       {
	 seq = new CSequence_Mol (datatype);
	 seq->set_taxon_and_sequence(*taxonnames[curr_taxon], *sequences[curr_taxon]);
	 seqData.push_back(seq);
	 delete taxonnames[curr_taxon];
	 delete sequences[curr_taxon];
       }
     } // End interleaved format
     

     // Do some final error checking

     if (taxaNum == 0)
     {
       // delete all faststrings in all_lines
       for (unsigned ind=0; ind < all_lines.size(); ++ind)
	 delete all_lines[ind];
       return 0;
     }

     unsigned len1 = posNum;
     unsigned len2 = seqData[0]->length();

     if (len1 == 0)
       len1 = seqData[0]->length();

     if (len1 != seqData[0]->length())
     {
       std::cerr << "Parse error while reading the input file:\nThe number of positions specified in the file header is not identical to the number of positons in the first sequence.\n"
	            "If this is an interleaved format it can also be that blank lines between blocks are missing." << std::endl;
       std::cerr << "The number of positions in the first sequence is: " << len2
		 << ", while the number specified in the header is " << posNum << "."<<  std::endl;
       exit(-66);
     }

     // Check that the number of positions is correct and equal for all sequences: 
     for (curr_taxon=1; curr_taxon < taxaNum; ++curr_taxon)
     {
       seq = seqData[curr_taxon]; 
       len2 = seq->length();

       if (len1 != len2)
       {
	 std::cerr << "Parse error while reading the input file:\nThe number of residues in sequence " << curr_taxon+1 << " differs from the number of residues in sequence 1. The sequence lengths for sequence 1 and sequence " << curr_taxon+1 
		   <<  " are respectively " << len1 << " and " << len2 << "." << std::endl;
       exit(-67);
       }
     }

     // Check that the data type is equal for all sequences and that 
     CSequence_Mol::DataTypesEnum dt_tmp;
     
     datatype = seqData[0]->get_datatype();
     for (curr_taxon=1; curr_taxon < taxaNum; ++curr_taxon)
     {
       dt_tmp = seqData[curr_taxon]->get_datatype();
	
       if (datatype != dt_tmp)
       {
	 std::cerr << "Warning: Not all sequences seem to have the same data type.\n";
	 datatype = CSequence_Mol::mixed;
       }
     }
    
     char ambig0, ambig_tmp;
     ambig0 = seqData[0]->get_ambiguity_character();
     for (curr_taxon=1; curr_taxon < taxaNum; ++curr_taxon)
     {
       ambig_tmp = seqData[curr_taxon]->get_ambiguity_character();
	
       if (ambig0 != ambig_tmp)
       {
	 std::cerr << "Warning: Not all sequences seem to have the same ambiguity character.\n";
       }
     }
  


 
/*      // PrintMessage_cerr(seq.getName());PrintMessage_cerr("-\n");PrintMessage_cerr(seq.getFullName());PrintMessage_cerr("-\n");PrintMessage_cerr(seq.getDescription());PrintMessage_cerr("-\n"); */
/*      if ( infile.fail() && !infile.eof() ) */
/*      { */
/*        PrintMessage_cerr("\n\n"); */
/*        faststring errormsg; */
/*        errormsg   =  "An error occurred while reading the input file. It might not be a valid fasta file.\n"; */
/*        errormsg  +=  "File position: line "; */
/*        errormsg  +=  faststring(infile.line()).c_str(); */
/*        ErrorMessage(errormsg.c_str()); */
/*        PrintMessage_cerr("\n"); */
/*        flush_cerr(); */
/*        return -25; */
/*      } */
     
/*      if ( count_seq < global_sequence_from || count_seq > global_sequence_to ) */
/*      { */
/*        if (false) */
/*        { */
/* 	 PrintMessage_cerr("Skipping sequence "); */
/* 	 PrintMessage_cerr(faststring(count_seq).c_str()); */
/* 	 PrintMessage_cerr(": "); */
/* 	 PrintMessage_cerr(seq->getName()); */
/* 	 PrintMessage_cerr("\n"); */
/* 	 flush_cerr(); */
/*        } */
/*        continue; */
/*      } */
     
/*      PrintMessage_cerr("Processing sequence "); */
/*      PrintMessage_cerr(faststring(count_seq).c_str()); */
/*      PrintMessage_cerr(": "); */
/*      PrintMessage_cerr(seq->getName()); */
/*      PrintMessage_cerr("\n"); */
/*      flush_cerr(); */
     
/*      seqData.push_back(seq); */
/*    } // End while */
   
     unsigned j;
     
     for (j=0; j<all_lines.size(); ++j)
       delete all_lines[j];
     
     return 0;
   }  // End this element function.


   void ExportSequences(std::ostream &os, char format, unsigned interleaved_len)
   {
     if (taxaNum == 0)
       return;

     std::vector<faststring> vofs;

     unsigned i;

     for (i=0; i< taxaNum; ++i)
     {
       vofs.push_back(faststring());
       seqData[i]->getSequence_fill_in_gaps_and_stars(vofs[i]);
     }

     if (format=='p') // Phylip
     {
       std::vector<faststring> uphynames;
       unique_phylip_names(uphynames);

       unsigned pos=0;
       unsigned i;

       os << "\t" << taxaNum << "\t" << posNum << std::endl;

       while (pos < posNum)
       {
	 for (i=0; i<taxaNum; ++i)
	 {
           os << uphynames[i] << " ";
           os << vofs[i].substr(pos, interleaved_len) << std::endl;
	 }
	 os << std::endl;
	 pos += interleaved_len;
       }
     }
     else if (format=='f')
     {
       unsigned pos=0;
       unsigned i;

       for (i=0; i<taxaNum; ++i)
       {
         os << ">" << seqData[i]->getFullName() << std::endl;
	 pos = 0;
	 while (pos < posNum)
	 {
           os << vofs[i].substr(pos, interleaved_len) << std::endl;
	   pos += interleaved_len;
	 }
       }
     }
     else if (format=='n')
     {


/*        Begin data; */
/*        Dimensions ntax=4 nchar=15; */
/*        Format datatype=dna symbols="ACTG" missing=? gap=-; */
/*        Matrix */
/* 	 Species1   atgctagctagctcg */
/* 	 Species2   atgcta??tag-tag */
/* 	 Species3   atgttagctag-tgg */
/* 	 Species4   atgttagctag-tag            */
/* 	 ; */
/*        End; */

       os << "Begin data;" << std::endl;
       os << "Dimensions ntax="<< taxaNum << " nchar="<< posNum << ";" << std::endl;
       os << "Format datatype="<< seqData[0]->type_as_string()
	 //  << " symbols=\" "<< <<  "\"" 
	 << " missing=" << ambig_char << " gap=-;" << std::endl;
       os << "matrix" << std::endl;
       
       unsigned pos=0;
       unsigned i;

       while (pos < posNum)
       {
	 for (i=0; i<taxaNum; ++i)
	 {
           os << seqData[i]->getName() << " ";
           os << vofs[i].substr(pos, interleaved_len) << std::endl;
	 }
	 pos += interleaved_len;
       }
       os << ";" << std::endl << "end;" << std::endl;

     }


   }

   //   void            GetBaseFrequencies(unsigned, std::vector<unsigned>);
   
   //   void            initSplits();
   //   void            computeSplits();
   //   CSplits         *splits;
   //   CBasefreq       *basefreq;

   double compute_site_pattern_stat()
   {
     unsigned i, j;
     
     const char **tax = new const char* [taxaNum];

     for (j=0; j < taxaNum; ++j)
     {
       tax[j] = seqData[j]->getSeqStr();
     }

     faststring tmp;
     std::map<faststring, unsigned> m;

     for (i=0; i<posNum; ++i)
     {
       tmp.clear();

       for (j=0; j < taxaNum; ++j)
       {
	 tmp.push_back(tax[j][i]);
       }
       add_or_count(m, tmp);
     }
     std::map<faststring, unsigned>::iterator it, it_end;
     it     = m.begin();
     it_end = m.end();

     double res=0;
     double x;

     while (it != it_end)
     {
       x = it->second;
       res += x*log(x);
       //       std::cout << x << " " << res << std::endl;
       ++it;
     }
     delete [] tax;

     x = posNum;
     res -= x*log(x);
     //     std::cout << x << " " << res << std::endl;

     return res;
   }

 
   void get_vector_of_site_pattern_frequencies(fastvector<unsigned> &spf)
   {
     unsigned i, j;
     
     const char **tax = new const char* [taxaNum];

     for (j=0; j < taxaNum; ++j)
     {
       tax[j] = seqData[j]->getSeqStr();
     }

     faststring tmp;
     std::map<faststring, unsigned> m;

     for (i=0; i<posNum; ++i)
     {
       tmp.clear();

       for (j=0; j < taxaNum; ++j)
       {
	 tmp.push_back(tax[j][i]);
       }
       add_or_count(m, tmp);
     }
     std::map<faststring, unsigned>::iterator it, it_end;
     it     = m.begin();
     it_end = m.end();

     while (it != it_end)
     {
       spf.push_back(it->second);
       ++it;
     }
     delete [] tax;
   }

   void clear(CSequence_Mol::DataTypesEnum     datatype_param=CSequence_Mol::unknown)
   {
     taxaNum    = 0;
     posNum     = 0;
     ambig_char = '?';
     datatype   = datatype_param;

     originalPosNumbers.clear();

     int i, n=seqData.size();

     // Recently resolved memory leak:
     for (i=0; i<n; ++i)
         delete seqData[i];

     seqData.clear();
   }


   void add_seq_to_alignment(CSequence_Mol::DataTypesEnum     datatype_param,
                             const faststring                  &fullname_param,
                             const faststring                  &seq_data_param,
			     char                              ambig_param='\0' // default auto
			    )
   {
     if (datatype != datatype_param)
     {
	 std::cerr << "Warning: Not all sequences seem to have the same data type.\n";
	 datatype = CSequence_Mol::mixed;
     }

     // Handle general ambig TODO XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     //     if (general_ambig == CSequence_Mol::unknown)

     CSequence_Mol   *seq;
     seq = new CSequence_Mol (datatype_param, ambig_param);
     seq->set_taxon_and_sequence(fullname_param, seq_data_param);


     ++taxaNum;
     if (posNum == 0)
       posNum = seq->length();
     else if (posNum != seq->length())
     {
       std::cerr << "Error: Sequence added has different length than existing sequences." << std::endl;
       delete seq;
       return;
     }

     seqData.push_back(seq);

     // Further autodetect stuff???????


   }


   void alignment_without_gap_N_positions(CSequences2 &seqs)
   {
     seqs.clear(datatype);

     std::vector<bool> gap_N_positions(posNum, false);

     size_t i, j;
     const char *the_seq;
     char c;


     if (datatype == CSequence_Mol::dna)
     {
       for (i=0; i< taxaNum; ++i)
       {
	 the_seq = seqData[i]->getSeqStr();
	 for(j=0; j<posNum; ++j)
	 {
	   c = the_seq[j];
	   if (is_DNA_iupac_ambig(c) || c == '-')
	   {
	     gap_N_positions[j] = true;
	   }
	 }
       }
     }

     if (datatype == CSequence_Mol::protein)
     {
       for (i=0; i< taxaNum; ++i)
       {
	 the_seq = seqData[i]->getSeqStr();
	 for(j=0; j<posNum; ++j)
	 {
	   c = the_seq[j];
	   if (is_aa_ambig(c) || c == '-')
	   {
	     gap_N_positions[j] = true;
	   }
	 }
       }
     }

     faststring new_seq;

     for (i=0; i< taxaNum; ++i)
     {
       new_seq.clear();
       the_seq = seqData[i]->getSeqStr();

       for(j=0; j<posNum; ++j)
       {
	 if (!gap_N_positions[j])
	   new_seq.push_back(the_seq[j]);
       }
       seqs.add_seq_to_alignment(datatype, seqData[i]->getFullName(), new_seq, seqData[i]->get_ambiguity_character() );
     }
   }

};



// inline CSequences::CSequences():posNum(0),datatype(0),splits(NULL){}


//inline string    CSequences::GetTaxonLabel(unsigned i)
//                             {return taxa.GetTaxonLabel(i);}
//inline void      CSequences::SetTaxonLabel(unsigned i, string s)
//                             {taxa.SetTaxonLabel(i,s);}






/* inline void PrintSpecial_TaxaSeq(ostream& os, */
/* 			  const std::vector<unsigned>& taxa_vec, */
/* 			  const std::vector<unsigned>& pos_vec, */
/* 			  const CTaxa       *ptaxa, */
/* 			  const CSequences2 *pseq, */
/* 			  const char*      symbolsList, */
/* 			  bool       printTaxonLabels, */
/* 			  unsigned   taxonLabelPrintWidth, */
/* 			  unsigned   maxTaxonLabelLength) */
/* { */
/*   std::vector<unsigned>::const_iterator taxa_vec_it; */
/*   std::vector<unsigned>::const_iterator taxa_vec_end = taxa_vec.end(); */

/*   std::vector<unsigned>::const_iterator pos_vec_it; */
/*   std::vector<unsigned>::const_iterator pos_vec_end = pos_vec.end(); */

/*   unsigned theTaxonNumber; */

/*   std::ios_base::fmtflags orig_ios_flag = os.flags(); */
/*   os.flags(orig_ios_flag | ios::left); */

/*   for (taxa_vec_it=taxa_vec.begin(); taxa_vec_it!=taxa_vec_end; ++taxa_vec_it) */
/*   { */
/*     theTaxonNumber = *taxa_vec_it; */

/*     if (printTaxonLabels) */
/*     { */
/*       os.width(maxTaxonLabelLength); */
/*       faststring tl = ptaxa->GetTaxonLabel(theTaxonNumber); */
/*       tl.shorten(maxTaxonLabelLength); */
/*       os << tl; */
/*       int morespaces = (int)taxonLabelPrintWidth-(int)maxTaxonLabelLength; */

/*       if ( morespaces > 0 ) */
/* 	{ */
/* 	  os.width(morespaces); */
/* 	  os << " "; */
/* 	} */
/*     } */

/*     for (pos_vec_it = pos_vec.begin(); pos_vec_it != pos_vec_end; ++pos_vec_it) */
/*     { */
/*       os << symbolsList[pseq->GetChar(theTaxonNumber,*pos_vec_it)]; */
/*     } */
/*     os << std::endl; */
/*   } */
/*   os.flags(orig_ios_flag); */
/* } */

/* inline  */
/* void PrintSpecial_PosNumbers(ostream& os, */
/* 			     std::vector<unsigned> pos_vec, // we get a real copy */
/* 			     const CSequences2 *pseq, */
/* 			     unsigned         taxonLabelPrintWidth, */
/* 			     bool             UseNexusComments) */
/* { */
/*   unsigned i; */
/*   unsigned numPos           = pos_vec.size(); */
/*   unsigned printNoZerosTill = numPos; */
/*   bool     printNoZero_set  = false; */

/*   // We determine the original position numbers: */
/*   for (i=0; i<numPos; ++i) */
/*     pos_vec[i] = pseq->GetOriginalPosNumber(pos_vec[i]); */

/*   unsigned maxDigits = ((unsigned)log10((double)pos_vec[numPos-1]))+1; */
/*   unsigned factor    = 1; */
/*   unsigned digit; */

/*   for (i=1; i<maxDigits; ++i) */
/*     factor *= 10; */

/*   for (; factor; factor /= 10, printNoZero_set = false) */
/*   { */
/*     if (UseNexusComments) */
/*     { */
/*       os.width(taxonLabelPrintWidth); */
/*       os << "[ "; */
/*     } */
/*     else */
/*     { */
/*       os.width(taxonLabelPrintWidth); */
/*       os << ' '; */
/*     } */
/*     for (i=0; i < numPos; ++i) */
/*     { */
/*       digit = pos_vec[i]/factor; */

/*       if ( i >= printNoZerosTill || digit > 0) */
/*       { */
/* 	if (!printNoZero_set) */
/* 	{ */
/* 	  printNoZerosTill = i; */
/* 	  printNoZero_set = true; */
/* 	} */
/* 	os << (char)(digit+'0'); */
/* 	pos_vec[i] -= factor*digit; */
/*       } */
/*       else */
/*       { */
/* 	os << ' '; */
/*       } */
/*     } */
/*     if (UseNexusComments) */
/*     { */
/*       os << "]" << std::endl; */
/*     } */
/*     else */
/*     { */
/*       os << std::endl; */
/*     } */
/*   } */
/* } */


/* inline void      CSequences::GetBaseFrequencies(unsigned ti, */
/* 						std::vector<unsigned> bf) */
/* { */
/*   unsigned pos; */
/*   for (pos = 0; pos < posNum; ++pos) */
/*     ++bf[seqData[ti][pos]] */
/* } */


#endif
