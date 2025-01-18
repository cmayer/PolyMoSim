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

#ifndef MATH_EXPRESSION_PARSER_H
#define MATH_EXPRESSION_PARSER_H


#include <iostream>
#include "faststring3.h"
#include "special_functions.h"
#include <vector>
#include <cstdlib>
#include <cstring>


// Operations:
// precedence   operators       associativity
// 2            * /             left to right
// 3            + -             left to right
// 4            =               right to left^
// 5            ^               right to left (right associative)
// 100          # $             - since unary

static const char * ALLOWED_FUNCTION_NAMES[] = {"+", "-", "*", "/", "^", "#", "$", "Dist_Gauss", "Dist_Gamma", "Dist_Uniform", "Heavy", "Gamma", "Dist_Beta", "Dist_Cauchy", "Dist_Uniform_Interval"};
static unsigned     NUMBER_OF_ARGS[]         = { 2,   2,   2,   2,   2,   1,   1,   3,            3,            1,              1,       1,       3,           3,            3                };
static unsigned     N_FUNC                   = 15;

// Determine the expected number of arguments for this function.
int get_number_of_args(faststring fname)
{
  unsigned i;
  const char *cfname = fname.c_str();

  for (i=0; i<N_FUNC; ++i)
  {
    if (std::strcmp(cfname, ALLOWED_FUNCTION_NAMES[i]) == 0)
    {
      return NUMBER_OF_ARGS[i];
    }
  }

  //  if (fname.size() == 1)

  return -1; // Nothing found
}


inline double op_add(double x, double y)
{
  return x+y;
}

inline double op_sub(double x, double y)
{
  return x-y;
}

inline double op_mult(double x, double y)
{
  return x*y;
}

inline double op_div(double x, double y)
{
  return x/y;
}

inline double op_pow(double x, double y)
{
  return pow(x,y);
}

inline double unary_minus(double x)
{
  return -x;
}

inline double unary_plus(double x)
{
  return x;
}

inline void print_vec(std::ostream &os, std::vector<faststring> &v, char delim= '\n')
{
  size_t i=0, n=v.size();

  for (; i<n; ++i)
  {
    os << v[i].c_str() << delim;
  }
}

inline void print_vec(std::ostream &os, std::vector<double*> &v, char delim= '\n')
{
  size_t i=0, n=v.size();

  for (; i<n; ++i)
  {
    os << *v[i] << delim;
  }
}

inline int operator_preced(const char c)
{
    switch(c)
    {
        case '!':  // not implemented
            return 4;
        case '*':  case '/':
            return 3;
        case '+': case '-':
            return 2;
        case '=':   // not implemented
            return 1;
    case '^':
            return 5;
        case '#':
            return 100;
        case '$':
            return 100;
    }
    return 0;
}

inline bool operator_is_left_assoc(const char c)
{
    switch(c)
    {
      // left to right OK:
    case '*': case '/': case '+': case '-':
      return true;
      // right but not left
    case '=': case '!': case '^':
            return false;
    }
    return false;
}

inline bool is_operator(char c)
{
  return (c == '+' || c == '-' || c == '/' || c == '*' || c == '^' || c == '#' || c == '$');
}

inline bool is_upper(char c)
{
  return (c>='A' && c <= 'Z');
}

inline bool is_function(char c)
{
  return is_upper(c);
}

inline bool is_lower(char c)
{
  return (c>='a' && c <= 'z');
}

inline bool is_identifier(char c)
{
  return is_lower(c) || isdigit(c);
}


class operation
{
 protected:
  double *res;

 public:
  operation(double *pres):res(pres)
  {}

  virtual void compute()=0;

  double *get_res()
  {
    return res;
  }

  virtual ~operation()
  {
    delete res;
  }

};

class operation1:public operation
{
 private:
  // Function pointer to operation associated to this object.
  double (*op)(double);
  double *x;

 public:
 operation1(double (*p_op)(double), double *pres, double*px):operation(pres),op(p_op), x(px)
  {}

  virtual void compute()
  {
    *res = op(*x);
  }
};


class operation2:public operation
{
 private:
  // Function pointer to operation associated to this object.
  double (*op)(double, double);
  double *x1, *x2;

 public:
 operation2(double (*p_op)(double, double), double *pres, double*px1, double *px2):operation(pres),op(p_op),x1(px1),x2(px2)
  {}

  virtual void compute()
  {
    *res = op(*x1, *x2);
  }
};

class operation3:public operation
{
 private:
  // Function pointer to operation associated to this object.
  double (*op)(double, double, double);
  double *x1, *x2, *x3;

 public:
 operation3(double (*p_op)(double, double, double), double *pres, double*px1, double *px2, double *px3):operation(pres),op(p_op),x1(px1),x2(px2),x3(px3)
  {}

  virtual void compute()
  {
    *res = op(*x1, *x2, *x3);
  }
};


// Copy construktor not available
// Assignment opertor not available

class Cexpression_evaluator
{
  operation** op_list;

  double *double_memory;
  //  double *final_result;
  std::vector<double *> x_positions;

  unsigned num_ops_reserved;
  unsigned num_double_memory_reserved;

  unsigned num_ops;
  unsigned num_double_memory;

 public:
 Cexpression_evaluator(): op_list(NULL), double_memory(NULL), // final_result(NULL),
    num_ops_reserved(0), num_double_memory_reserved(0), num_ops(0), num_double_memory(0)
  {}

  ~Cexpression_evaluator()
  {
    //    std::cerr << "Entering Cexpression_evaluator destructor" << std::endl;

    unsigned i;
    for (i=0; i<num_ops; ++i)
      delete op_list[i];

    if (op_list)
      delete [] op_list;


    if (double_memory)
      delete [] double_memory;
    //    std::cerr << "Exiting Cexpression_evaluator destructor" << std::endl;
  }

  Cexpression_evaluator(const Cexpression_evaluator&);
  Cexpression_evaluator& operator=(const Cexpression_evaluator&);


  void add_op(operation *p_o)
  {
    if (num_ops >= num_ops_reserved)
    {
      std::cout << "An internal error occurred when trying to add an operation to the array of operations.\n"
	        << "The array is already full." << std::endl;
      exit(-11);
    }

    op_list[num_ops] = p_o;
    ++num_ops;
  }

  double compute_for_indeterminate_x(double param_x)
  {
    unsigned i, n= (unsigned)x_positions.size();

    if (op_list == NULL || double_memory == NULL)
    {
      std::cerr << "Error: Calling  compute_for_indeterminate_x, but memory not initialized." << std::endl;
      exit(-1);
    }

    for (i=0; i<n; ++i)
      *(x_positions[i]) = param_x;

    if (num_ops != 0)
    {
      for (i=0; i<num_ops; ++i)
      {
	op_list[i]->compute();
      }
      return *(op_list[num_ops-1]->get_res());
    }
    else
    {
      //      std::cerr << "ww" << num_double_memory << std::endl;
      if (num_double_memory == 1)
	return double_memory[0];
      else
      {
	std::cerr << "Internal error in function compute_for_indeterminate_x. Expression cannot be evaluted." << std::endl;
	exit(-1);
      }
    }

  }


  void reserve(unsigned p_num_double_memory, unsigned p_num_ops)
  {
    if (op_list)
    {
      //      std::cerr << "Deleting op_list in reserve" << std::endl;
      delete [] op_list;
    }
    if (double_memory)
    {
      //      std::cerr << "Deleting double_memory in reserve" << std::endl;
      delete [] double_memory;
    }

    num_double_memory_reserved = p_num_double_memory;
    num_ops_reserved   = p_num_ops;

    double_memory = new double[num_double_memory_reserved];
    op_list    = new operation*[num_ops_reserved];
  }


  void set_double_memory_pos_x(unsigned i)
  {
    x_positions.push_back((double *)(double_memory + i) );
  }

  void print_x_positons()
  {
    unsigned i=0, n = (unsigned)x_positions.size();

    std::cerr << "Size: " << n << std::endl;

    for (i=0; i< n; ++i)
    {
      std::cerr << x_positions[i] << " " << *x_positions[i] << std::endl;
    }
  }

  // Returns the index at which value was inserted.
  unsigned append_double_memory(double x)
  {
    double_memory[num_double_memory] = x;
    ++num_double_memory;
    return num_double_memory-1;
  }


  double get_double_memory(unsigned i)
  {
    if (i < num_double_memory_reserved)
    {
      return double_memory[i];
    }
    else
    {
      std::cerr << "Critical error. I do not know how to proceed in get_set_double_memory."
		<< std::endl;
      std::cerr << "num_double_memory_reserved: " << num_double_memory_reserved << std::endl;
      std::cerr << "index i:          "           << i << std::endl;
      exit(-2);
    }
  }

  double* get_pointer_double_memory(size_t i)
  {
    if (i < num_double_memory_reserved)
      return (double_memory + i);
    else
    {
      std::cerr << "Critical error. I do not know how to proceed in get_pointer_double_memory."
		<< std::endl;
      std::cerr << "num_double_memory_reserved: " << num_double_memory_reserved << std::endl;
      std::cerr << "index i:          " << i << std::endl;
      exit(-2);
    }
  }
  
  size_t get_actual_num_double_memory()
  {
    return num_double_memory;
  }

  size_t get_actual_num_ops()
  {
    return num_ops;
  }

  size_t get_reserved_num_double_memory()
  {
    return num_double_memory_reserved;
  }

  size_t get_reserved_num_ops()
  {
    return num_ops_reserved;
  }


/*   expression_evaluator(size_t reserve_ops, size_t reserve_double_memory) */
/*   { */
    
/*     op_list = new operation[reserve_ops]; */
/*     double_memory = new double[reserve_ops]; */

/*   } */

};

void read_float(faststring &input, size_t &i, faststring &id)
{
  size_t i_start = i, n = input.size();
  short  dot_count = 0;
  short  e_count   = 0;

  if (i>=n)
    return;

  while ( isdigit(input[i]) || (input[i]=='.' && dot_count<2) ||
	                       (toupper(input[i]) =='E' && dot_count<2) )
  {
    if (input[i]=='.')
    {
      ++dot_count;
    }
    if (toupper(input[i])=='E')
    {
      ++e_count;
      ++i;
      if (i==n)
	break;
    }
    else
    {
      ++i;
    }
    if (i==n)
      break;
  }
  id = input.substr(i_start, i-i_start);
}

void read_indentifier(faststring &input, size_t &i, faststring &id)
{
  size_t i_start = i;
  
  while ( (input[i]>='a' && input[i]<='z') ||  (input[i]>='0' && input[i]<='9') || input[i]=='_')
  {
    ++i;
  }

  id = input.substr(i_start, i-i_start);
}

bool read_function(faststring &input, size_t &i, faststring &id)
{
  size_t i_start = i;
  
  while ((input[i] >= 'A' && input[i] <= 'Z') || input[i]=='_' || (input[i] >= 'a' && input[i] <= 'z') )
  {
    ++i;
  }

  id = input.substr(i_start, i-i_start);

  if (get_number_of_args(id)==-1)
    return false;
  else
    return true;
}

class math_expression_parser
{
  faststring                  expression;
  std::vector<faststring>     infix;
  std::vector<faststring>     postfix;
  Cexpression_evaluator       exp_eval;

 public:
  math_expression_parser(faststring exp):expression(exp)
  {
  }

  math_expression_parser(const math_expression_parser&);
  math_expression_parser& operator=(const math_expression_parser&);

  ~math_expression_parser()
  {
    //    std::cout << "Entering math_expression_parser destructor" << std::endl;
  }

  void print_orig(std::ostream &os)
  {
    os << "Original math expression:" << std::endl;
    os << expression << std::endl;
  }

  // Make sure you know what you do if you use this function:
  Cexpression_evaluator* get_exp_eval_pointer()
  {
    return &exp_eval;
  }

  void print_postfix(std::ostream &os)
  {
    os << "Postfix expression (token by token):" << std::endl;
    print_vec(os, postfix);
  }

  void print_infix(std::ostream &os)
  {
    os << "Infix expression (token by token):" << std::endl;
    print_vec(os, infix);
  }

  void tokenize()
  {
    size_t  i=0, n = expression.size();
    char c;
    faststring tmp;
    
    expression.skip_spaces(i);
    c = expression[i];
    
    while (i < n)
    {
      // read one char
      expression.skip_spaces(i);
      c = expression[i];
      
      if (isdigit(c))                // float
      {
	read_float(expression, i, tmp);
	infix.push_back(tmp);
      }
      else if (islower(c))  // identifier
      {
	read_indentifier(expression, i, tmp);
	infix.push_back(tmp);
      }
      else  if (c >= 'A' && c <='Z') // Function name
      {
	if (!read_function(expression, i, tmp) )
	{
	  // TODO: Better error handling
	  printf("Error: Unknown function found in expression.\n");
	  exit(-1);
	}
	infix.push_back(tmp);
      }
      else  if (is_operator(c) || c == '(' || c == ')' || c == ',')
      {
	infix.push_back(faststring(c));
	++i;
      }
    }
  }

  void detect_unary_signs_in_infix()
  {
    size_t i, n = infix.size();
    
    // Unary +: #
    // Unary -: $
    
    if (n == 0)
      return;
    
    if (infix[0]=="-")
      infix[0] = "$";
    else if (infix[0]=="+")
      infix[0] = "#";
    
    i = 1;
    while (i<n)
    {
      if (     infix[i]=="-" && ( infix[i-1] == "(" || is_operator(infix[i-1][0]) || is_upper(infix[i-1][0]) || infix[i-1][0] == ',' ) )
	infix[i] = "$";
      else if (infix[i]=="+" && ( infix[i-1] == "(" || is_operator(infix[i-1][0]) || is_upper(infix[i-1][0]) || infix[i-1][0] == ',' ) )
	infix[i] = "#";
      ++i;
    }
  }


  // shunting_yard algorithm
  bool infix2postfix()
  {
    detect_unary_signs_in_infix();
    
    size_t  i=0, n = infix.size();
    char c;
 
    std::vector<faststring> stack;     // operator stack
    
    faststring     sc;          // used for record stack element
    faststring     tmp;

    while(i < n)
    {
      c = infix[i][0];
      // If the token is a number (identifier), then add it to the output queue.
      if (isdigit(c))
      {
	postfix.push_back(infix[i]);
      }
      else if (c >= 'a' && c <='z') // Identifier, i.e. variable names
      {
	postfix.push_back(infix[i]);
      }
      else
	// If the token is a function token, then push it onto the stack.
	if(c >= 'A' && c <='Z')
	{
	  stack.push_back(infix[i]);
	}
      // If the token is a function argument separator (e.g., a comma):
	else if(c == ',')
	{
	  bool pe = false;
	  while(stack.size() > 0)
	  {
	    sc = stack.back();
	    if(sc == "(")
	    {
	      pe = true;
	      break;
	    }
	    else
	    {
	      // Until the token at the top of the stack is a left parenthesis,
	      // pop operators off the stack onto the output queue.
	      postfix.push_back(sc);
	      stack.pop_back();
	    }
	  }
	  // If no left parentheses are encountered, either the separator was misplaced
	  // or parentheses were mismatched.
	  if (!pe)
	  {
	    printf("Error: separator or parentheses mismatched\n");
	    return false;
	  }
	}
      // If the token is an operator, op1, then:
	else if(is_operator(c))
	{
	  while(stack.size() > 0) 
	  {
	    sc = stack.back();
	    // While there is an operator token, o2, at the top of the stack
	    // op1 is left-associative and its precedence is less than or equal to that of op2,
	    // or op1 is right-associative and its precedence is less than that of op2,
	    if ( is_operator(sc[0]) &&
		 ( (operator_is_left_assoc(c) && (operator_preced(c) <= operator_preced(sc[0]))) ||
		  (!operator_is_left_assoc(c) && (operator_preced(c) <  operator_preced(sc[0]))))
		 )
	    {
	      // Pop o2 off the stack, onto the output queue;
	      postfix.push_back(sc);
	      stack.pop_back();
	    }
	    else
	    {
	      break;
	    }
	  }
	  // push op1 onto the stack.
	  stack.push_back(faststring(c));
	}
      // If the token is a left parenthesis, then push it onto the stack.
	else if(c == '(')
	{
	  stack.push_back("(");
	}
      // If the token is a right parenthesis:
	else if(c == ')')
	{
	  bool pe = false;
	  // Until the token at the top of the stack is a left parenthesis,
	  // pop operators off the stack onto the output queue
	  while(stack.size() > 0)
	  {
	    sc = stack.back();
	    if (sc[0] == '(')
	    {
	      pe = true;
	      break;
	    }
	    else
	    {
	      postfix.push_back(sc);
	      stack.pop_back();
	    }
	  }
	  // If the stack runs out without finding a left parenthesis, then there are mismatched parentheses.
	  if(!pe) 
	  {
	    printf("Error: parentheses mismatched\n");
	    return false;
	  }
	  // Pop the left parenthesis from the stack, but not onto the output queue.
	  stack.pop_back();
	  // If the token at the top of the stack is a function token, pop it onto the output queue.
	  if(stack.size() > 0)
	  {
	    sc = stack.back();
	    if(sc[0]>= 'A' && sc[0]<='Z')
	    {
	      postfix.push_back(sc);
	      stack.pop_back();
	    }
	  }
	}
	else
	{
	  printf("Unknown token %c\n", c);
	  return false; // Unknown token
	}
      ++i;
    } // END while (i<n)
    
    // There are no more tokens to read
    // While there are still operator tokens in the stack:
    while(stack.size() > 0)
    {
      sc = stack.back();
      if(sc[0] == '(' || sc[0] == ')')
      {
	printf("Error: parentheses mismatched\n");
	return false;
      }
      postfix.push_back(sc);
      stack.pop_back();
    }
    return true;
  } // END infix_to_postfix()


  // OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
  bool postfix2Cexpression_evaluator()
  {
    unsigned count_identifier = 0;
    unsigned count_op_or_func  = 0;


    //    printf("Individual operations:\n");
    
    if (postfix.size() == 0)
    {
      printf("Error: Can't determine operations without having the postfix stack.\n");
    }
 
    int i, n= (int)postfix.size();
    char   c;
    std::vector<faststring> stack_tokens;
    std::vector<double*>    stack_double_memory;
    
    // First pass - count the number of identifiers and operations
    i=0;
    while (i < n)
    {
      // Read the next token from input.
      c = postfix[i][0];  // First char of token determines its kind

      if (is_identifier(c) )
      {
	++count_identifier;
      }

      if (is_operator(c) || is_function(c))
      {
	++count_op_or_func;
      }
      //      printf("Found: %c -- %u %u\n", c, count_identifier, count_op_or_func);
      ++i;
    }

    //    printf("Found: Identifier: %u Operations: %u\n", count_identifier, count_op_or_func);
    exp_eval.reserve(count_identifier+count_op_or_func, count_op_or_func);
    
    //     printf("Actual:   %u %u\n", exp_eval.get_actual_num_double_memory() ,   exp_eval.get_actual_num_ops()); 
    //     printf("Reserved: %u %u\n", exp_eval.get_reserved_num_double_memory() , exp_eval.get_reserved_num_ops()); 

    unsigned identifier_pos;

    // Second pass: assign identifiers, numbers and results of operations to elemens in the double_memory field:
    i=0;
    while (i < n)
    {
      // Read the next token from input.
      c = postfix[i][0];  // First char of token determines its kind

      if (is_identifier(c) )
      {
	if (isdigit(c))
	{
	  identifier_pos = exp_eval.append_double_memory(postfix[i].ToDouble());
	  stack_double_memory.push_back(exp_eval.get_pointer_double_memory(identifier_pos));
	}
	else if (postfix[i] == "x")
	{
	  identifier_pos = exp_eval.append_double_memory(-7);
	  exp_eval.set_double_memory_pos_x(identifier_pos);

	  stack_double_memory.push_back(exp_eval.get_pointer_double_memory(identifier_pos));
	}
	else
	{
	  std::cout << "Error: Unknown identifier found in expression: "
		    << postfix[i] << std::endl;
	  exit(-1);
	}
      }

      if (is_operator(c) || is_function(c))
      {
	// Reserve the next memory position for this operation:
	identifier_pos = exp_eval.append_double_memory(-13);
	stack_double_memory.push_back(exp_eval.get_pointer_double_memory(identifier_pos));
      }
      ++i;
    }

    //    std::cout << "Postfix" << std::endl;
    //    print_vec(std::cout, postfix);
    //    std::cout << "stack_double_memory" << std::endl;
    //    print_vec(std::cout, stack_double_memory);

    std::vector<char> already_used_dm(stack_double_memory.size(), ' ');

    // Third pass: build list of operations
    i = 0;
    while (i < n)
    {
      // Read the next token from input.
      c = postfix[i][0];  // First char of token determines its kind

      // If the token is an operator  (operator here includes both operators, and functions).
      if (is_operator(c) || is_function(c))
      {
	unsigned num_arguments = get_number_of_args( postfix[i]);

	if (num_arguments == 1)
	{
	  int index1 = (int)i-1;
	  while (already_used_dm[index1] == 'u')
	  {
	    --index1;
	  }
	  double *p1  = stack_double_memory[index1];
	  double *res = stack_double_memory[i];
	  already_used_dm[index1] = 'u';

	  if (postfix[i] == "#")
	  {
	    operation1 *p_op1 = new operation1(unary_plus, res, p1);
	    exp_eval.add_op(p_op1);
	  }
	  else
	  if (postfix[i] == "$")
	  {
	    operation1 *p_op1 = new operation1(unary_minus, res, p1);
	    exp_eval.add_op(p_op1);
	  }
	  if (postfix[i] == "Dist_Uniform")
	  {
	    operation1 *p_op1 = new operation1(distribution_uniform, res, p1);
	    exp_eval.add_op(p_op1);
	  }
	  else
	  if (postfix[i] == "Heavi")
	  {
	    operation1 *p_op1 = new operation1(heaviside, res, p1);
	    exp_eval.add_op(p_op1);
	  }
	  else
	  if (postfix[i] == "Gamma")
	  {
	    operation1 *p_op1 = new operation1(local_gamma, res, p1);
	    exp_eval.add_op(p_op1);
	  }
	  else
	  {
	    std::cerr << "Error: Unknown function found in expression: " << postfix[i] << std::endl;
	    exit(-14);
	  }
	}
	else
	if (num_arguments == 2)
	{
	  unsigned index1 = i-1;
	  while (already_used_dm[index1] == 'u')
	  {
	    --index1;
	  }
	  double *p2  = stack_double_memory[index1];
	  already_used_dm[index1] = 'u';

	  --index1;
	  while (already_used_dm[index1] == 'u')
	  {
	    --index1;
	  }
	  double *p1  = stack_double_memory[index1];
	  already_used_dm[index1] = 'u';
	  double *res = stack_double_memory[i];

	  if (postfix[i] == "+")
	  {
	    operation2 *p_op2 = new operation2(op_add, res, p1, p2);
	    exp_eval.add_op(p_op2);
	  }
	  else
	  if (postfix[i] == "-")
	  {
	    operation2 *p_op2 = new operation2(op_sub, res, p1, p2);
	    exp_eval.add_op(p_op2);
	  }
	  else
	  if (postfix[i] == "*")
	  {
	    operation2 *p_op2 = new operation2(op_mult, res, p1, p2);
	    exp_eval.add_op(p_op2);
	  }
	  else
	  if (postfix[i] == "/")
	  {
	    operation2 *p_op2 = new operation2(op_div, res, p1, p2);
	    exp_eval.add_op(p_op2);
	  }
	  else
	  if (postfix[i] == "^")
	  {
	    operation2 *p_op2 = new operation2(op_pow, res, p1, p2);
	    exp_eval.add_op(p_op2);
	  }
	  else
	  {
	    std::cerr << "Error: Unknown function found in expression: " << postfix[i] << std::endl;
	    exit(-14);
	  }
	}
	else
	if (num_arguments == 3)
	{
	  unsigned index1 = i-1;
	  while (already_used_dm[index1] == 'u')
	  {
	    --index1;
	  }
	  double *p3  = stack_double_memory[index1];
	  already_used_dm[index1] = 'u';

	  --index1;
	  while (already_used_dm[index1] == 'u')
	  {
	    --index1;
	  }
	  double *p2  = stack_double_memory[index1];
	  already_used_dm[index1] = 'u';

	  --index1;
	  while (already_used_dm[index1] == 'u')
	  {
	    --index1;
	  }
	  double *p1  = stack_double_memory[index1];
	  already_used_dm[index1] = 'u';

	  double *res = stack_double_memory[i];

	  if (postfix[i] == "Dist_Gauss")
	  {
	    operation3 *p_op3 = new operation3(distribution_gauss, res, p1, p2, p3);
	    exp_eval.add_op(p_op3);
	  }
	  else
	  if (postfix[i] == "Dist_Gamma")
	  {
	    operation3 *p_op3 = new operation3(distribution_gamma, res, p1, p2, p3);
	    exp_eval.add_op(p_op3);
	  }
	  else
	  if (postfix[i] == "Dist_Beta")
	  {
	    operation3 *p_op3 = new operation3(distribution_beta, res, p1, p2, p3);
	    exp_eval.add_op(p_op3);
	  }
	  else
	  if (postfix[i] == "Dist_Cauchy")
	  {
	    operation3 *p_op3 = new operation3(distribution_cauchy, res, p1, p2, p3);
	    exp_eval.add_op(p_op3);
	  }
	  else
	  if (postfix[i] == "Dist_Uniform_Interval")
	  {
	    operation3 *p_op3 = new operation3(distribution_uniform_interval, res, p1, p2, p3);
	    exp_eval.add_op(p_op3);
	  }
	  else
	  {
	    std::cerr << "Error: Unknown function found in expression: " << postfix[i] << std::endl;
	    exit(-14);
	  }
	}
	else // Currently not supported
	{
	    std::cerr << "Error: Unknown function found in expression: " << postfix[i] << std::endl;
	    exit(-14);
	}
      } // End if (is_operator(c) || is_function(c))
      // If this is not an opertor, we do nothing.

      ++i;
    } // End while (i < n)
      

      // If there are more values in the stack
      // (Error) The user input has too many values.
      return false;
  }

  double eval(double x)
  {
    return exp_eval.compute_for_indeterminate_x(x);
  }
  

  // This part is has not been finished yet. I do not even know the following is correct.
  // In several sources people give infix 2 prefix algorithms which are indeed infix 2 post fix.

// Algorithm ConvertInfixtoPrefix
// Purpose: Convert and infix expression into a prefix expression. Begin // Create operand and operator stacks as empty stacks. Create OperandStack
// Create OperatorStack
// 1) While input expression still remains, read and process the next token.
// 2) while( not an empty input expression ) read next token from the input expression
// 3) Test if token is an operand or operator
// 3a)
//    if ( token is an operand ) // Push operand onto the operand stack. OperandStack.Push (token)
// 3b)
//    else // Token must be an operator.
// 3c)
//       if ( token is '(' or OperatorStack.IsEmpty() or OperatorHierarchy(token) > OperatorHierarchy(OperatorStack.Top()) ) // Push left parentheses onto the operator stack
//            OperatorStack.Push ( token )
// 3d)
//    else if( token is ')' ) // Continue to pop operator and operand stacks, building
//       prefix expressions until left parentheses is found.
//
//    Each prefix expression is push back onto the operand stack as either a left or right operand for the next operator. 
//    while( OperatorStack.Top() not equal '(' ) OperatorStack.Pop(operator) OperandStack.Pop(RightOperand) OperandStack.Pop(LeftOperand) operand = operator + LeftOperand + RightOperand OperandStack.Push(operand) endwhile
// // Pop the left parthenses from the operator stack. OperatorStack.Pop(operator)
// else if( operator hierarchy of token is less than or equal to hierarchy of top of the operator stack )
// // Continue to pop operator and operand stack, building prefix // expressions until the stack is empty or until an operator at // the top of the operator stack has a lower hierarchy than that // of the token. while( !OperatorStack.IsEmpty() and . OperatorHierarchy(token) lessThen Or Equal to OperatorHierarchy(OperatorStack.Top()) ) OperatorStack.Pop(operator) OperandStack.Pop(RightOperand) OperandStack.Pop(LeftOperand) operand = operator + LeftOperand + RightOperand OperandStack.Push(operand)
// endwhile // Push the lower precedence operator onto the stack OperatorStack.Push(token)
// endif
// endif endif endif endwhile // If the stack is not empty, continue to pop operator and operand stacks building // prefix expressions until the operator stack is empty. while( !OperatorStack.IsEmpty() ) OperatorStack.Pop(operator) OperandStack.Pop(RightOperand) OperandStack.Pop(LeftOperand) operand = operator + LeftOperand + RightOperand
// OperandStack.Push(operand) endwhile
// // Save the prefix expression at the top of the operand stack followed by popping // the operand stack.
// print OperandStack.Top()
// OperandStack.Pop()
// End

/*
void infix2prefix(faststring infix, faststring prostfix)
{
  unsigned i=0, n=infix.size();

  vector<char> operator_stack;
  vector<char> operands_stack;

  // 1)
  while (i<n)
  {
    // Skips spaces //2)
    infix.skip_spaces(i);

    // is operand // 3a)
    if ( isdigit(infix[i]) )
    {
      operands_stack.push_back(infix[i]);
    }

    // 3b
    if ( is_in_symbol_list__(infix[i],"+-/()*") )
    {
    // 3c
      if (infix[i] == '(' || operator_stack.empty() || precendence(infix[i]) > precendence( operator_stack.get_last()  ) )
	operator_stack.push_back('(');
      else
    // 3d
	if (infix[i] == ')')
	{
	  while (operator_stack.get_last() != '(' )
	  {
	    postfix.push_back(operator_stack.get_last_and_pop());
	  }
	  postfix.push_back(operator_stack.get_last_and_pop());
	}



      while
	postfix.push_back(infix[i]);
      ++i;
    }


    if( infix[i] == '(' )
    {
      stack.push_back(infix[i]);
      ++i;
    }
    
    

      if(infix[i] is ')' )
      {
         while( element at top of stack != '(' )
         {
            p = element at top of stack;
            increment p to next character in postfix exp<b></b>ression;
         }
         increment i to next character in infix exp<b></b>ression;
      }
 
      if( i is an operator )
      {
         if( stack is empty )
             stack_push(i);
         else
         {
 
            while priority(element at top of stack) >= priority(i)
            {
               p = element at top of stack;
               increment p to next character in postfix exp<b></b>ression;
            }
            stack_push(i);
         }
         increment i to next character in infix exp<b></b>ression;
      }
   }
   while stack is not empty
   {
      p = stack_pop()
      increment p to next character in postfix exp<b></b>ression;
   }
   p = '\0';
}
*/

  

};

#endif
