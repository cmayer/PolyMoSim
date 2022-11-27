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

//-----------------------------------------------------------------------------
// Author:    Christoph Mayer
//            Gropiusweg 13
//            44801 Bochum
//
// Copyright: Christoph Mayer
//
// Description:
//      A fast vector class for objects which do not require their constructors
//         or destructors to be called.
//      This class uses malloc, realloc and free to allocate and
//         free dynamic memory for the "vector". Especially the efficient
//         realloc has no C++ equivalent - the major drawback of new and
//         delete.
//      Disadvantage: No constructors and destructors are called for
//         objects in the vector. This limits is applicability.
//         The fastvector-nd.h uses new and delete and does not use
//         memcpy. It is therefore well suited for all classes.
//         Unfortunately, its less efficient.
//



#ifndef FASTVECTOR_H
#define FASTVECTOR_H

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cctype>
// #include <iomanip>

#define _MINBUFFER_CAPACITY 4

template<typename T>
class fastvector
{
  // Possibly even more efficient:
  //    Use pointer _end, _capacity_end instead of _len, _capacity


 private:
  T           *_buf;
  unsigned    _len;
  unsigned    _capacity;   // Capacity of container, i.e number of bytes availible
                           // for writing. If a capacity is requested by the user,
                           // an additional byte is reserved in order to be able to append
                           // the additional \0 which is appended by some function calls.
                           // Consequently, the capacity reported to user is _capacity-1.
 protected:


 public:
  ~fastvector ()
  {
    if (_buf)
      //      delete [] _buf;
      std::free(_buf);
  }

  fastvector ():_buf(NULL),  _len(0), _capacity(0)
  {
  }

  fastvector (const T *v, unsigned len)
  {
    //    std::cout << "New fastvector 2" << std::endl;
    _len      = len;
    _capacity = len;

    if (len)
    {
      _buf      = (T *) malloc(_capacity*sizeof(T)); 
      //_buf      = new T [_capacity]; 
      memcpy(_buf, v, len*sizeof(T));
    }
    else
    {
      _buf = NULL;
    }
  }

  fastvector (const T *v_begin, const T *v_end)
  {
    //    std::cout << "New fastvector 3" << std::endl;
    _len      = (v_end - v_begin);
    _capacity = _len;

    if (_len)
    {
      _buf      = (T *)malloc(_capacity*sizeof(T)); 
      // _buf      = new T [_capacity]; 
      memcpy(_buf, v_begin, _len*sizeof(T));
    }
    else
    {
      _buf = NULL;
    }
  }

  fastvector (const fastvector &a):_len(a._len),_capacity(a._len)
  {
    if (_len)
    {
      _buf = (T *)malloc(_capacity*sizeof(T));
      memcpy(_buf, a._buf, _len*sizeof(T));
    }
    else
    {
      _buf = NULL;
    }      
  }


  void push_back(const T &c)
  {
    if (_capacity == _len)
    {
      if (_buf == NULL)
      {
	_capacity = _MINBUFFER_CAPACITY;
	_buf = (T *)malloc(_capacity*sizeof(T));
      }
      else
      {
	T *tmp;

	_capacity <<= 1;
	tmp = (T *)realloc(_buf, _capacity*sizeof(T));

	if (tmp == NULL)
	{
	  std::cout << "Program aborded due to failed realloc!" << std::endl; 
	  exit(0);
	}
	_buf = tmp;
      }
    }
    _buf[_len] = c;
    ++_len;
  }


  T *begin() const
  {
    return _buf;
  }


  T *end() const
  {
    return _buf+_len;
  }


  T *rbegin() const
  {
    return _buf + (_len - 1);
  }


  T *rend() const
  {
    return _buf - 1;
  }


  void assign(const T *v, int len)
  {
    _len = len;
    reserve(len);
    memcpy(_buf, v, len*sizeof(T));
  }


  void assign(const T *v_begin, const T *v_end)
  {
    _len = v_end - v_begin;
    reserve(_len);
    memcpy(_buf, v_begin, _len*sizeof(T));
  }

  void assign(const fastvector &a)
  {
    _len = a._len;
    reserve(_len);
    memcpy(_buf, a._buf, _len*sizeof(T));
  }


  void reserve(unsigned s)
  {
    //    std::cout << "reserve" << std::endl;
    if (_capacity < s)
    {
      if (_buf == NULL)
      {
	//_buf = new T [s];
	_buf = (T *)malloc(s*sizeof(T));
      }
      else
      {
	T *tmp;

	tmp = (T *)realloc(_buf, s*sizeof(T));
	if (tmp == NULL)
	{
	  std::cout << "Program aborded due to failed realloc!" << std::endl; 
	  exit(0);
	}
	_buf = tmp;
      }
      _capacity = s;
    }
  }


  unsigned size() const
  {
    return _len;
  }


  unsigned capacity() const
  {
    return _capacity; // One char is reserved for the \0 at end of string
  }


  void clear()
  {
    //        std::cout << "clear" << std::endl;
    _len = 0;
  }


  void reset()
  {
    //    std::cout << "reset" << std::endl;
    if (_buf)
    {
      // delete [] _buf;
      free(_buf);
      _buf = NULL;
      _capacity = 0;
      _len = 0;
    }
  }


  bool check_pos(unsigned pos)
  {
    return (pos < _len);
  }


  void set_unckecked(unsigned pos, const T &c)
  {
    _buf[pos] = c;
  }


  bool set(unsigned pos, const T &c)
  {
    if ( check_pos(pos) )
    {
      _buf[pos] = c;
      return true;
    }
    else
    {
      return false;
    }
  }


  void get_unckecked(unsigned pos, T &c)
  {
    c = _buf[pos];
  }


  T get_unckecked(unsigned pos)
  {
    return _buf[pos];
  }



  bool get(unsigned pos, T &c)
  {
    if ( check_pos(pos) )
    {
      c = _buf[pos];
      return true;
    }
    else
    {
      return false;
    }
  }

  fastvector &operator=(const fastvector &a)
    {
      _len = a._len;
      if (_len)
      {
	reserve(_len);
	memcpy(_buf, a._buf, _len*sizeof(T));
      }
      return *this;
    }

  void swap(fastvector &a)
  {
    T         *_tmpbuf      = _buf;             _buf = a._buf;           a._buf = _tmpbuf; 
    unsigned   _tmplen      = _len;             _len = a._len;           a._len = _tmplen; 
    unsigned   _tmpcapacity = _capacity;   _capacity = a._capacity; a._capacity = _tmpcapacity;     
  }

};


#endif
