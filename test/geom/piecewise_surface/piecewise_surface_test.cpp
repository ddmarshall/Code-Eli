/*********************************************************************************
* Copyright (c) 2013 David D. Marshall <ddmarsha@calpoly.edu>
*
* All rights reserved. This program and the accompanying materials
* are made available under the terms of the Eclipse Public License v1.0
* which accompanies this distribution, and is available at
* http://www.eclipse.org/legal/epl-v10.html
*
* Contributors:
*    David D. Marshall - initial code and implementation
********************************************************************************/

#include <cstdlib>  // EXIT_SUCCESS
#include <ctime>    // time functions
#include <cstring>  // std::strlen, std::strncmp
#include <string>   // std::string
#include <fstream>  // std::fstream
#include <memory>   // std::auto_ptr
#include <ctime>    // strftime

#include <cpptest.h> // CppTest Framework

#include "piecewise_surface_test_suite.hpp"  // piecewise_surface_test_suite

enum TestType {testTypeText, testTypeCompiler, testTypeHTML};

int main(int argc, char *argv[])
{
  Test::Output *pout;
  TestType tt(testTypeCompiler);
  int rtn(EXIT_SUCCESS);

  // if user passed command line argument then parse it
  if (argc==1)
  {
    tt=testTypeCompiler;
  }
  else if (argc==2)
  {
    size_t lenargv(strlen(argv[1]));

    // using short codes
    if (lenargv==2)
    {
      if (argv[1][0]=='-')
      {
        switch(argv[1][1])
        {
          case('t'):
            tt=testTypeText;
            break;
          case('c'):
            tt=testTypeCompiler;
            break;
          case('h'):
            tt=testTypeHTML;
            break;
          default:
            std::cerr << "Unrecognized parameter: " << argv[1] << std::endl;
            rtn=EXIT_FAILURE;
            break;
        }
      }
      else
      {
        std::cerr << "Unrecognized parameter: " << argv[1] << std::endl;
        rtn=EXIT_FAILURE;
      }
    }
    else
    {
      if (lenargv==6)
      {
        if (std::strncmp(argv[1], "--text", lenargv)==0)
          tt=testTypeText;
        else if (std::strncmp(argv[1], "--html", lenargv)==0)
          tt=testTypeHTML;
        else if (std::strncmp(argv[1], "--help", lenargv)==0)
          rtn=EXIT_FAILURE;
        else
        {
          std::cerr << "Unrecognized parameter: " << argv[1] << std::endl;
          rtn=EXIT_FAILURE;
        }
      }
      else if (lenargv==10)
      {
        if (std::strncmp(argv[1], "--compiler", lenargv)==0)
          tt=testTypeCompiler;
        else
        {
          std::cerr << "Unrecognized parameter: " << argv[1] << std::endl;
          rtn=EXIT_FAILURE;
        }
      }
    }
  }
  else
  {
    std::cerr << "Can only pass one argument." << std::endl;
    rtn=EXIT_FAILURE;
  }
  if (rtn==EXIT_FAILURE)
  {
    std::cerr << "  Usage: " << argv[0] << " [--text|-t|--compiler|-c|--html|-h]" << std::endl;
    return rtn;
  }

  // allocate the correct output type
  switch (tt)
  {
    case(testTypeText):
    {
      Test::TextOutput *po = new Test::TextOutput(Test::TextOutput::Verbose);
      pout = static_cast<Test::Output *>(po);
      break;
    }
    case(testTypeCompiler):
    {
      Test::CompilerOutput *po = new Test::CompilerOutput(Test::CompilerOutput::GCC);
      pout = static_cast<Test::Output *>(po);
      break;
    }
    case(testTypeHTML):
    {
      Test::HtmlOutput *po = new Test::HtmlOutput;
      pout = static_cast<Test::Output *>(po);
      break;
    }
    default:
    {
      rtn=EXIT_FAILURE;
      break;
    }
  }

  if (rtn == EXIT_SUCCESS)
  {
    // get the current time information
    time_t rawtime;
    struct tm * timeinfo;
    std::string datetime;
    time (&rawtime);
    timeinfo = localtime(&rawtime);
    datetime = asctime(timeinfo);

    // add the test suites to the test runner
    // NOTE: This is where changes should be needed
    //
    Test::Suite ts;
    std::string ostr_filename("piecewise_surface_test_results.html");

    // add the cppack test suites
    ts.add(std::auto_ptr<Test::Suite>(new piecewise_surface_test_suite<float>()));
    ts.add(std::auto_ptr<Test::Suite>(new piecewise_surface_test_suite<double>()));
    ts.add(std::auto_ptr<Test::Suite>(new piecewise_surface_test_suite<long double>()));

    //
    // NOTE: End of section that should be changed

    // run the test
    rtn = (ts.run(*pout)) ? EXIT_SUCCESS : EXIT_FAILURE;

    // generate the html data if requested
    if (tt==testTypeHTML)
    {
      std::ofstream ostr(ostr_filename.c_str());
      Test::HtmlOutput *phtmlout=dynamic_cast<Test::HtmlOutput *>(pout);

      phtmlout->generate(ostr, true, datetime);

      ostr.close();
    }

    delete pout;
  }

  return rtn;
}
