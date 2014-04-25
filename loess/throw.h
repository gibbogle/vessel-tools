#ifndef THROW_EXCEPT_H
#define THROW_EXCEPT_H

#include <stdexcept>
#include <sstream>
#include "varkitversion.h"
#ifdef __GNUC__
#include <execinfo.h> /* compile with -rdynamic */
#include <cstdlib>


#define THROW(a) do{std::ostringstream _os;\
	_os << "\nRevision:" <<   VARKIT_REVISION << "\nFile : "<<  __FILE__ << "\nLine  : "<< __LINE__ << "\nMethod : " << __FUNCTION__ << "\nWhat : " << a << std::endl;\
	 void   *_array[30];\
	std::size_t  _size = ::backtrace(_array, 30);\
	char    **_trace_strings = ::backtrace_symbols(_array, _size);\
	if(_trace_strings!=NULL)\
	    {\
	    _os << "StackTrace:\n";\
	    for (std::size_t _i = 0; _i < _size; ++_i)\
		   {\
		   _os << "\t" << _trace_strings[_i] << std::endl;\
		   }\
	    _os << std::endl;\
	    std::free(_trace_strings);\
	    }\
	throw std::runtime_error(_os.str());\
	}while(0)
#else

#define THROW(a) do{std::ostringstream _os;\
	_os << "\nRevision:" <<   VARKIT_REVISION << "\nFile : "<<  __FILE__ << "\nLine  : "<< __LINE__ << "\nMethod : " << __FUNCTION__ << "\nWhat : " << a << std::endl;\
	throw std::runtime_error(_os.str());\
	}while(0)

#endif

#endif
