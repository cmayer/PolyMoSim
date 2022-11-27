#include <iostream>
#include <string>
#include <vector>
#include "faststring2.h"

using namespace std;

typedef faststring mystring;


class cartesian_product
{
 private:
  

  void cartesian_product_for_2strings(const vector<mystring> &a, const vector<mystring> &b, vector<mystring> &res)
  {
    res.clear();
    for (unsigned i=0; i<a.size(); ++i)
    {
      for (unsigned j=0; j<b.size(); ++j)
      {
	mystring tmp = a[i] + b[j];
	res.push_back(tmp);
      }
    }
  }

  void nth_cartesian_product_of_symbolset(unsigned replicates, const mystring symbols, vector<mystring> &res)
  {
    res.clear();

    vector<mystring> start, tmp_res, tmp_res2;
    for (unsigned i=0; i < symbols.size(); ++i)
    {
      // string:
      //    mystring tmp = mystring(1u,symbols[i]);
      // faststring:
      mystring tmp = mystring(symbols[i], 1u);
      
      start.push_back(tmp);
    }
    
    //    cout << "DEBUG: start" << endl;
    //    print_container(cout, start, "", ",","");
    //    cout << "DEBUG: start -- END" << endl;
    
    tmp_res = start;
    for (unsigned j=1; j<replicates; ++j)
    {
      cartesian_product_for_2strings(tmp_res, start, tmp_res2);
      tmp_res = tmp_res2;
    }

    res = tmp_res2;
  }


 public:
  void operator()(unsigned replicates, const mystring symbols, vector<mystring> &res)
  {
    nth_cartesian_product_of_symbolset(replicates, symbols, res);
  }
};



