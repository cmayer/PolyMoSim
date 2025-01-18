/***************************************************************************************************
*  The PolyMoSim project is distributed under the following license:
*  
*  Copyright (c) 2006-2025, Christoph Mayer, Leibniz Institute for the Analysis of Biodiversity Change,
*  Bonn, Germany
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

#include "tclap/CmdLine.h"


// #define GLOBAL_VARS
#include "global-types-and-parameters.hpp"
#include "PolyMoSim.h"
#include <string>
#include <ctime>

using namespace TCLAP;
using namespace std;

string                      global_tree_file;
string                      global_model_file;

string                      global_preanalysis_file;
string                      global_postanalysis_file;
string                      global_log_file;
bool                        global_logging;
string                      global_output_filename;

string                      global_ancestral_sequence_file;
string                      global_siterateshist_file;
string                      global_siteratesdata_file;

unsigned                    global_num_repetitions;
bool                        global_use_GUI;
unsigned                    global_seed_random_generator;
enum_rg                     global_rg;
enum_outputformat           global_outputformat;

unsigned                    global_verbosity;


void ReportErrorMessage(const char *e)
{
  if (global_use_GUI)
  {

  }
  else
  {
    fputs(e, stderr); //    fprintf(stderr, e);
  }
}


void myPrint(FILE *of, const char *s1)
{
  if (of == stderr && global_use_GUI)
  {
  }
  else
  {
    fputs(s1, of);  //    fprintf(of, s1);
  }
}

void myPrint(FILE *of, const char *s1, const char *s2)
{
  if (of == stderr && global_use_GUI)
  {
  }
  else
  {
    fputs(s1, of);  //    fprintf(of, s1);
    fputs(s2, of);  //    fprintf(of, s2);
  }
}

void myPrint(FILE *of, const char *s1, const char *s2, const char *s3)
{
  if (of == stderr && global_use_GUI)
  {
  }
  else
  {
    fputs(s1, of);  //    fprintf(of, s1);
    fputs(s2, of);  //    fprintf(of, s2);
    fputs(s3, of);  //    fprintf(of, s3);
  }
}

void myPrint(FILE *of, const char *s1, const char *s2, const char *s3, const char *s4)
{
  if (of == stderr && global_use_GUI)
  {}
  else
  {
    fputs(s1, of);  //    fprintf(of, s1);
    fputs(s2, of);  //    fprintf(of, s2);
    fputs(s3, of);  //    fprintf(of, s3);
    fputs(s4, of);  //    fprintf(of, s4);
  }
}




void myPrint(ostream &os, const char *s1)
{
  if (os.rdbuf() == cerr.rdbuf() && global_use_GUI)
  {
  }
  else
  {
    os << s1;
  }
}

void myPrint(ostream &os, const char *s1, const char *s2)
{
  if (os.rdbuf() == cerr.rdbuf() && global_use_GUI)
  {
  }
  else
  {
    os << s1 << s2;
  }
}

void myPrint(ostream &os, const char *s1, const char *s2, const char *s3)
{
  if (os.rdbuf() == cerr.rdbuf() && global_use_GUI)
  {
  }
  else
  {
    os << s1 << s2 << s3;
  }
}

void myPrint(ostream &os, const char *s1, const char *s2, const char *s3, const char *s4)
{
  if (os.rdbuf() == cerr.rdbuf() && global_use_GUI)
  {}
  else
  {
    os << s1 << s2 << s3 << s4;
  }
}




// void PrintMessage(FILE *of, const char *s1)
// {
//   fprintf(of, s1);
// }

// void PrintMessage(FILE *of, const char *s1, const char *s2)
// {
//   fprintf(of, s1);
//   fprintf(of, s2); 
// }

// void PrintMessage(FILE *of, const char *s1, const char *s2, const char *s3)
// {
//   fprintf(of, s1);
//   fprintf(of, s2); 
//   fprintf(of, s3);
// }

// void PrintMessage(FILE *of, const char *s1, const char *s2, const char *s3, const char *s4)
// {
//   fprintf(of, s1);
//   fprintf(of, s2); 
//   fprintf(of, s3);
//   fprintf(of, s4);
// }





void init_param()
{
  global_num_repetitions       = 1;
  global_use_GUI               = false;
  global_seed_random_generator = (unsigned)time(NULL);
  global_rg                    = rg_MT19937;
  global_outputformat          = outputformat_fasta;
  global_log_file              = ""; // string("PolyMoSim_") + number2str((unsigned)time(NULL)) + ".log";
  global_logging               = false;

  global_verbosity             = 1; // Reason: 0: no output, 1: only report seed, which we should report by default. 
}


int read_and_init_parameters(int argc, char** argv)
{
  string errormsg;

  init_param();

  try
  {
    CmdLine cmd("This program simulates the evolution of molecular (nucleotide and amino acid) sequences.",
		' ', VERSION);


    vector<string> allowed;
    allowed.push_back("nexus");
    allowed.push_back("phylip");
    allowed.push_back("phylip_no_spaces");
    allowed.push_back("fasta");
    allowed.push_back("site_pattern_freq_absolute");
    allowed.push_back("site_pattern_freq_relative");
    allowed.push_back("site_pattern_freq_absolute_fill");
    allowed.push_back("site_pattern_freq_relative_fill");

    ValueArg<string> global_outputformat_Arg("f", "outputFormat",
       "Output format of sequence data. Default: fasta. The site_pattern* format list the site pattern frequencies instead of the alignment. The site_pattern*_fill formats list all patterns the non fill formats only the site patterns that occured in the simulated data set.",
       false, "fasta", allowed);
    cmd.add( global_outputformat_Arg );

    ValueArg<unsigned> global_seed_random_generator_Arg("s", "seed",
       "The seed value for the random number generator. Default: time.",
       false, global_seed_random_generator, "unsigned int");
    cmd.add( global_seed_random_generator_Arg );

//     SwitchArg global_use_GUI_Arg("", "gui",
// 	"Use graphical user interface.",
//        global_use_GUI);
//     cmd.add( global_use_GUI_Arg );

//     SwitchArg global_use_NOGUI_Arg("", "nogui",
// 	"Do not use graphical user interface.",
//        !global_use_GUI);
//     cmd.add( global_use_NOGUI_Arg );

    ValueArg<string> global_preanalysis_file_Arg("", "pre",
	"File included in output before each generated data set.",
	false, "", "string");
    cmd.add( global_preanalysis_file_Arg );

    ValueArg<string> global_postanalysis_file_Arg("", "post",
	"File included in output after each generated data set.",
	false, "", "string");
    cmd.add( global_postanalysis_file_Arg );

    ValueArg<string> global_ancestral_sequence_file_Arg("", "print_ancestral_seq",
	"With this option and by providing a file name, the ancestral sequence is printed to this file. Default: ancestral sequence is not printed",
	false, "", "string");
    cmd.add( global_ancestral_sequence_file_Arg );

    ValueArg<string> global_siterateshist_file_Arg("", "print_siterate_histogram",
	"Print site rates histogram for each model to given file. Default: No siterate information is printed.",
	false, "", "string");
    cmd.add( global_siterateshist_file_Arg );

    ValueArg<string> global_siteratesdata_file_Arg("", "print_siterate_data",
	"Print full site rates of all models to given file. Default: No siterate information is printed.",
	false, "", "string");
    cmd.add( global_siteratesdata_file_Arg );

    ValueArg<unsigned> global_num_repetitions_Arg("n", "nreps",
	"Number of independent data generated in simulation. "
	"Default: 1",
	false, global_num_repetitions, "unsigned");
    cmd.add( global_num_repetitions_Arg );

    ValueArg<string> global_log_file_Arg("l", "log",
	"File to write log information to.",
	false, global_log_file, "string");
    cmd.add( global_log_file_Arg );

    ValueArg<string> output_filename_Arg("o", "outfile",
	"Name of output file. If not specified, results are printed to standard output.",
	false, global_output_filename, "string");
    cmd.add( output_filename_Arg );

    ValueArg<string> global_tree_file_Arg("t", "treefile",
	"Tree file for simulation.",
	true, global_tree_file, "string");
    cmd.add( global_tree_file_Arg );

    ValueArg<string> global_model_file_Arg("m", "modelfile",
	"Model file for simulation.",
	true, global_model_file, "string");
    cmd.add( global_model_file_Arg );

    ValueArg<unsigned> global_verbosity_file_Arg("", "verbosity",
	"Adjust the level of additional information given to the user. Values from 0-5 are valid. Default 1. Set to 0 for less output.",
	false, global_verbosity, "integer");
    cmd.add( global_verbosity_file_Arg );

    cmd.parse( argc, argv );

    // Assigning parameters to variables:
    //    global_input_filename          = infilename_Arg.getValue();
    global_output_filename       = output_filename_Arg.getValue();

    global_model_file            = global_model_file_Arg.getValue();
    global_tree_file             = global_tree_file_Arg.getValue();
    global_log_file              = global_log_file_Arg.getValue();
    if (global_log_file != "")
      global_logging = true;

    global_ancestral_sequence_file =     global_ancestral_sequence_file_Arg.getValue();
    global_siterateshist_file      =     global_siterateshist_file_Arg.getValue();
    global_siteratesdata_file      =     global_siteratesdata_file_Arg.getValue();

    global_preanalysis_file      = global_preanalysis_file_Arg.getValue();
    global_postanalysis_file     = global_postanalysis_file_Arg.getValue();

    global_num_repetitions       = global_num_repetitions_Arg.getValue();
    //    global_use_GUI               = global_use_GUI_Arg.getValue();
    //    global_use_GUI               = global_use_NOGUI_Arg.getValue();
    global_seed_random_generator = global_seed_random_generator_Arg.getValue();
    //    global_rg                    = global_rg_Arg.getValue();
    global_verbosity             = global_verbosity_file_Arg.getValue();

    if (global_outputformat_Arg.getValue() == "nexus")
    {
      global_outputformat = outputformat_nexus;
    }
    else if ( global_outputformat_Arg.getValue() == "phylip" )
    {
      global_outputformat = outputformat_phylip;
    }
    else if ( global_outputformat_Arg.getValue() == "phylip_no_spaces" )
    {
      global_outputformat = outputformat_phylip_no_spaces;
    }
    else if ( global_outputformat_Arg.getValue() == "site_pattern_freq_absolute" )
    {
      global_outputformat = outputformat_pattern_absolute;
    }
    else if ( global_outputformat_Arg.getValue() == "site_pattern_freq_relative" )
    {
      global_outputformat = outputformat_pattern_relative;
    }
    else if ( global_outputformat_Arg.getValue() == "site_pattern_freq_absolute_fill" )
    {
      global_outputformat = outputformat_pattern_absolute_fill;
    }
    else if ( global_outputformat_Arg.getValue() == "site_pattern_freq_relative_fill" )
    {
      global_outputformat = outputformat_pattern_relative_fill;
    }
    else
    {
      global_outputformat = outputformat_fasta;
    }

  }
  catch (TCLAP::ArgException &e)
  {
    errormsg = string("Error: ") + e.error().c_str() + " for arg " + e.argId().c_str() + "\n";
    ReportErrorMessage(errormsg.c_str() );
    return (-1);
  }

  //Checking rest of parameters:
  // Do some immediate error checking
  if ( !FileExists(global_model_file.c_str() ) )
  {
    errormsg = string("Error: the model file \"") + global_model_file.c_str() + "\" does not exist.\n";
    ReportErrorMessage(errormsg.c_str());
    return -2;
  }

  if ( !FileExists(global_tree_file.c_str() ) )
  {
    errormsg = string("Error: the tree file \"") + global_tree_file.c_str() + "\" does not exist.\n";
    ReportErrorMessage(errormsg.c_str());
    return -2;
  }

  if ( !global_preanalysis_file.empty() && !FileExists(global_preanalysis_file.c_str() ) )
  {
    errormsg = string("Error: the pre-analysis file \"") + global_preanalysis_file.c_str() + "\" does not exist.\n";
    ReportErrorMessage(errormsg.c_str());
    return -2;
  }

  


  if ( !global_postanalysis_file.empty() && !FileExists(global_postanalysis_file.c_str() ) )
  {
    errormsg = string("Error: the post-analysis file \"") + global_postanalysis_file.c_str() + "\" does not exist.\n";
    ReportErrorMessage(errormsg.c_str());
    return -2;
  }


  return 0;
} // End: read_and_init_parameters(int argc, char** argv)



void print_analysis_parameters(FILE *of, const char *s)
{
  // TODO: not fully implwmwnted in this version.

  myPrint(of, s, "Parameters used in this simulation:\n");
  myPrint(of, s, "Name of model file:                            ", global_model_file.c_str(), "\n");
  myPrint(of, s, "Name of tree file:                             ", global_tree_file.c_str(), "\n");

  myPrint(of, s, "Name of pre-analysis file:                     ", global_preanalysis_file.c_str(), "\n");
  myPrint(of, s, "Name of post-analysis file:                    ", global_postanalysis_file.c_str(), "\n");
  myPrint(of, s, "Name of log file:                              ", global_log_file.c_str(), "\n");
  if (!global_output_filename.empty())
    myPrint(of, s, "Output filename:                               ", global_output_filename.c_str(), "\n");
  else
    myPrint(of, s, "Output filename:                               ", "Not specified. Output is printed to standard output.", "\n");

  myPrint(of, s, "Print ancestral sequence to file:              ");
  if (!global_ancestral_sequence_file.empty() )
    myPrint(of, global_ancestral_sequence_file.c_str(), "\n");
  else
    myPrint(of,  "-\n");

  myPrint(of, s, "Print siterates histogram data to file:        ");
  if (!global_siterateshist_file.empty() )
  {
    myPrint(of, global_siterateshist_file.c_str(), "\n");
  }
  else
    myPrint(of,  "-\n");

  myPrint(of, s, "Print siterates data to file:                  ");
  if (!global_siteratesdata_file.empty() )
  {
    myPrint(of, global_siteratesdata_file.c_str(), "\n");
  }
  else
    myPrint(of,  "-\n");

  myPrint(of, s, "Logging in effect:                             ", (global_logging) ? "yes\n" : "no\n");
  myPrint(of, s, "Number of repetition:                          ", faststring(global_num_repetitions).c_str(), "\n");
  myPrint(of, s, "Random number seed: (default if not specified: ", faststring(global_seed_random_generator).c_str(), "\n");
  
  switch(global_rg)
  {
     case rg_MT19937:  myPrint(of, s, randomNumberGeneratorNames[global_rg], "\n"); break;
  }
  
  switch(global_outputformat)
  {
     case outputformat_nexus:                  myPrint(of, s, "Output format: nexus\n"); break;
     case outputformat_fasta:                  myPrint(of, s, "Output format: fasta\n"); break;
     case outputformat_phylip:                 myPrint(of, s, "Output format: phylip\n"); break;
     case outputformat_phylip_no_spaces:       myPrint(of, s, "Output format: phylip_no_spaces\n"); break;
     case outputformat_pattern_absolute:       myPrint(of, s, "Output format: site patterns absolute\n"); break;
     case outputformat_pattern_relative:       myPrint(of, s, "Output format: site patterns relative\n"); break;
     case outputformat_pattern_absolute_fill:  myPrint(of, s, "Output format: site patterns absolute and fill in missing patterns.\n"); break;
     case outputformat_pattern_relative_fill:  myPrint(of, s, "Output format: site patterns relative and fill in missing patterns.\n"); break;
     case outputformat_logfile:                myPrint(of, s, "Output format: special logfile format\n"); break;
  }

 myPrint(of, s, "Level of verbosity: ", faststring(global_verbosity).c_str(), "\n"); 

}


void print_analysis_parameters(ofstream &of, const char *s)
{
  // TODO: not fully implwmwnted in this version.

  myPrint(of, s, "Parameters used in this simulation:\n");
  myPrint(of, s, "Name of model file:                            ", global_model_file.c_str(), "\n");
  myPrint(of, s, "Name of tree file:                             ", global_tree_file.c_str(), "\n");

  myPrint(of, s, "Name of pre-analysis file:                     ", global_preanalysis_file.c_str(), "\n");
  myPrint(of, s, "Name of post-analysis file:                    ", global_postanalysis_file.c_str(), "\n");
  myPrint(of, s, "Name of log file:                              ", global_log_file.c_str(), "\n");

  myPrint(of, s, "Print ancestral sequence to file:              ");
  if (!global_ancestral_sequence_file.empty())
    myPrint(of, global_ancestral_sequence_file.c_str(), "\n");
  else
    myPrint(of,  "-\n");

  myPrint(of, s, "Print siterates histogram information to file: ");
  if (!global_siterateshist_file.empty())
  {
    myPrint(of, global_siterateshist_file.c_str(), "_data.txt, ");
    myPrint(of, global_siterateshist_file.c_str(), "_hist.txt\n");
  }
  else
    myPrint(of,  "-\n");

  myPrint(of, s, "Print siterates data to file:                  ");
  if (!global_siteratesdata_file.empty() )
  {
    myPrint(of, global_siteratesdata_file.c_str(), "\n" );
  }
  else
    myPrint(of,  "-\n");


  myPrint(of, s, "Logging in effect:                             ", (global_logging) ? "yes\n" : "no\n");
  myPrint(of, s, "Number of repetition:                          ", faststring(global_num_repetitions).c_str(), "\n");
  myPrint(of, s, "Random number seed: (default if not specified: ", faststring(global_seed_random_generator).c_str(), "\n");
  
  switch(global_rg)
  {
     case rg_MT19937:  myPrint(of, s, randomNumberGeneratorNames[global_rg], "\n"); break;
  }
  
  switch(global_outputformat)
  {
     case outputformat_nexus:                  myPrint(of, s, "Output format: nexus\n"); break;
     case outputformat_fasta:                  myPrint(of, s, "Output format: fasta\n"); break;
     case outputformat_phylip:                 myPrint(of, s, "Output format: phylip\n"); break;
     case outputformat_phylip_no_spaces:       myPrint(of, s, "Output format: phylip_no_spaces\n"); break;
     case outputformat_pattern_absolute:       myPrint(of, s, "Output format: site patterns absolute\n"); break;
     case outputformat_pattern_relative:       myPrint(of, s, "Output format: site patterns relative\n"); break;

     case outputformat_pattern_absolute_fill:  myPrint(of, s, "Output format: site patterns absolute and fill in missing patterns\n"); break;
     case outputformat_pattern_relative_fill:  myPrint(of, s, "Output format: site patterns relative and fill in missing patterns\n"); break;

     case outputformat_logfile:                myPrint(of, s, "Output format: special logfile format\n"); break;
  }

  myPrint(of, s, "Level of verbosity: ", faststring(global_verbosity).c_str(), "\n"); 


}






