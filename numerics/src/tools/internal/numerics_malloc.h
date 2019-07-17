/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/
#ifndef NUMERICS_MALLOC_H
#define NUMERICS_MALLOC_H
#ifdef NUMERICS_MALLOC
#include "numerics_verbose.h"
#include <stdlib.h>

static inline void *numerics_malloc(size_t size, const char* file, const int line,
                                    const char* function)
{
  void *p = malloc(size);
  if(!p && size > 0)
  {
    perror(function);
    numerics_printf("%s:%d",file,line);
    numerics_error(function, "malloc failed");
  }
  numerics_printf_verbose(1, "%s:%d in %s malloc %d bytes at %p\n", file, line, function, size, p);
  return(p);
}

static inline void *numerics_calloc(size_t nmemb, size_t size, const char* file,
                                    const int line, const char* function)
{
  void *p = calloc(nmemb, size);
  if(!p && nmemb*size > 0)
  {
    perror(function);
    numerics_printf("%s:%d",file,line);
    numerics_error(function, "calloc failed");
  }
  numerics_printf_verbose(1, "%s:%d in %s calloc %d bytes at %p\n", file, line, function, size*nmemb, p);
  return(p);
}

static inline void *numerics_realloc(void *ptr, size_t size, const char* file,
                                     const int line, const char* function)
{
  void *p = realloc(ptr, size);
  if(!p && size > 0)
  {
    perror(function);
    numerics_printf("%s:%d",file,line);
    numerics_error(function, "realloc failed");
  }
  numerics_printf_verbose(1, "%s:%d in %s realloc %d bytes at %p\n", file, line, function, size, p);
  return(p);
}


static inline void numerics_free(void *ptr, const char* file,
                                 const int line, const char* function)
{
  char mess[256];
  numerics_printf_verbose(1, "%s:%d in %s free %p\n", file, line, function, ptr);
  if (ptr)
    free(ptr);
  else
  {
    snprintf(mess, 256, "%s:%d", file, line);
    numerics_warning(function, mess);
    numerics_error(function, "free ptr twice");
  }
}
#undef malloc
#undef calloc
#undef free
#undef realloc
#define malloc(X) numerics_malloc(X,__FILE__,__LINE__,__FUNCTION__)
#define calloc(X,Y) numerics_calloc(X,Y,__FILE__,__LINE__,__FUNCTION__)
#define realloc(X,Y) numerics_realloc(X,Y,__FILE__,__LINE__,__FUNCTION__)
#define free(X) numerics_free(X,__FILE__,__LINE__,__FUNCTION__)
#endif
#endif
