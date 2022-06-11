
/*
 * Software License
 *
 * Copyright (c) 2001 Joshua M. Deutsch. All rights 
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * 3. The end-user documentation included with the redistribution, if
 *    any, must include the following acknowlegement:  
 *       "This product includes software developed Joshua M. Deutsch
 *        Dept of Physics University of California, Santa Cruz"
 *    Alternately, this acknowlegement may appear in the software itself,
 *    if and wherever such third-party acknowlegements normally appear.
 *
 * 4. The name "GESSES"
 *    must not be used to endorse or promote products derived
 *    from this software without prior written permission. For written 
 *    permission, please contact josh@physics.ucsc.edu.
 *
 * 5. Products derived from this software may not be called "gesses"
 *    nor may "gesses" appear in their names without prior written
 *    permission of Joshua M. Deutsch.
 *
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESSED OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED.  IN NO EVENT SHALL JOSHUA M. DEUTSCH OR THE UNIVERSITY
 * OF CALIFORNIA BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 * USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 * OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 * ====================================================================
 */
#ifndef DEFS_H

#define DEFS_H
#include <sys/types.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifndef __USE_GNU
#define __USE_GNU
#endif

typedef unsigned char uchar;

#define memcmp_(a,b,l) memcmp(a,b,l*sizeof(a[0]))

#define memcpy_(dest,src,l) memcpy((dest),(src),(l)*sizeof((src)[0]))
#define memmove_(dest,src,l) memmove((dest),(src),(l)*sizeof((src)[0]))
#define memset_(dest,int_const,l) memset((dest),(int_const),(l)*sizeof((dest)[0]))

#define free_(pointer)  do {if (pointer) free(pointer);} while(0)

#define calloc_(pointer,l) pointer = (typeof(pointer)) calloc((l),sizeof((pointer)[0])) 
#define malloc_(pointer,l) (pointer) = (typeof(pointer)) malloc((l)*sizeof((pointer)[0])) 

#define realloc_(pointer,l) pointer = (typeof(pointer)) realloc ((pointer), (l)*sizeof((pointer)[0]))

#define size_2d(pointer,i,j) do{   	                   \
   int _j;                                                 \
   calloc_ (pointer,(i));                                  \
   for(_j=0; _j < (i); _j++) calloc_((pointer)[_j],(j));       \
} while (0)


#define size_3d(pointer,i,j,k) do{   	                   \
   int _k;                                                 \
   calloc_ (pointer,(i));                                  \
   for(_k=0; _k < i; _k++) size_2d((pointer)[_k],j,k);     \
} while (0)

#define free_2d(pointer,i) do{                             \
   int _j;                                                 \
   for(_j=0; _j < i; _j++) free((pointer)[_j]);            \
   free (pointer);                                         \
} while (0)

#define free_3d(pointer,j,k) do{                             \
   int _l;                                                   \
   for(_l=0; _l < k; _l++) free_2d((pointer)[_l],j);         \
   free(pointer);                                            \
} while (0)

#define cp_3d(a,b,nx,ny,nz)  do{                             \
   int _i_,_j_;                                              \
   for(_i_=0;_i_<(nx);_i_++){                                \
      for(_j_=0;_j_<(ny);_j_++){                             \
	 memcpy_(a[_i_][_j_],b[_i_][_j_],(nz));              \
      }                                                      \
   }                                                         \
} while(0)
   
#define clone_3d(a,b,nx,ny,nz) do {                          \
   size_3d(a,(nx),(ny),(nz));                                \
   cp_3d(a,b,(nx),(ny),(nz));                                \
}  while(0)

#define alloca_(pointer,l) (pointer) = (typeof(pointer)) alloca((l)*sizeof((pointer)[0])); 

#define alloca_0(pointer,l) do {                                    \
   (pointer) = (typeof(pointer)) alloca((l)*sizeof((pointer)[0]));  \
   memset(pointer,0, l * sizeof((pointer)[0]));                     \
} while (0)
   
    

#define alloca_2d(pointer,i,j) do{                         \
   int _j;                                                 \
   alloca_ ((pointer),(i));                                \
   for(_j=0; _j < i; _j++) alloca_((pointer)[_j],j);       \
} while (0)

    

#define alloca_2d_0(pointer,i,j) do{                         \
   int _j;                                                   \
   alloca_0 ((pointer),(i));                                 \
   for(_j=0; _j < i; _j++) alloca_0((pointer)[_j],j);        \
} while (0)




void * concat(void * a1, size_t n1, void * a2, size_t n2);
void reverse_int_array(int * a, int length);
#define concat_(array_1,n_1,array_2,n2) ((typeof(array_1)) concat(array_1,n_1*sizeof(array_1[0]),array_2,n2*sizeof(array_2[0])))

#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define QUART(x) ((x)*(x)*(x)*(x))
#define QUINT(x) ((x)*(x)*(x)*(x)*(x))

#define SWAP(a,b) do {  	\
   typeof(a) tmp = (a); 	\
   (a)=b;  			\
   (b)=tmp;  			\
} while(0) 


#define ran() (random()/(RAND_MAX+1.0))

#define randint(lower,upper) ((int)((lower) + ((upper)-(lower))*ran()))

#define MN(a,b) ((a) < (b)?(a):(b))
#define MX(a,b) ((a) > (b)?(a):(b))

#define is_power_2(n) ((((n)&(n-1)) == 0))

#define newarr_(_pointer_,n) do {                                                         \
   void * _a_ = calloc((n)*sizeof((_pointer_)[0]) + sizeof(unsigned long),1);             \
   (_pointer_) = (typeof(_pointer_))(_a_ + sizeof(unsigned long));                        \
   *((unsigned long *) _a_) =  (n);                                                       \
} while(0)

#define new2darr_(_pointer_,i,j) do{   	                     \
   int _j;                                                   \
   newarr_ (_pointer_,(i));                                  \
   for(_j=0; _j < (i); _j++) newarr_((_pointer_)[_j],j);     \
} while (0)

#define new3darr_(_pointer_,i,j,k) do{   	                     \
   int _j;                                                           \
   newarr_ (_pointer_,(i));                                          \
   for(_j=0; _j < (i); _j++) new2darr_((_pointer_)[_j],j,k);         \
} while (0)

#define new4darr_(_pointer_,i,j,k,l) do{   	                     \
   int _j;                                                           \
   newarr_ (_pointer_,(i));                                          \
   for(_j=0; _j < (i); _j++) new3darr_((_pointer_)[_j],j,k,l);       \
} while (0)

#define new5darr_(_pointer_,i,j,k,l,m) do{   	                     \
   int _j;                                                           \
   newarr_ (_pointer_,(i));                                          \
   for(_j=0; _j < (i); _j++) new4darr_((_pointer_)[_j],j,k,l,m);     \
} while (0)




#define freearr_(_pointer_)   do {                                 \
   if (_pointer_)                                                  \
   {                                                               \
      free_(((void *) (_pointer_)) - sizeof(unsigned long));       \
   }                                                               \
   /* _pointer_=0;  */                                              \
} while(0)


#define size_(_pointer_) (((unsigned long *) ((void *) (_pointer_) - sizeof(unsigned long) ))[0])

#define resizearr_(_pointer_,_num_) do {                                                      \
   if (_pointer_) {                                                                           \
   int _n_ = size_(_pointer_);                                                                \
   void * _beg_ = ((void *) (_pointer_)) - sizeof(unsigned long);                             \
   _beg_ = realloc(_beg_,(_num_)*sizeof((_pointer_)[0]) + sizeof(unsigned long));             \
   *((unsigned long *) _beg_) =  (_num_);                                                     \
   (_pointer_) = (typeof(_pointer_))(_beg_ + sizeof(unsigned long));                          \
   if (_num_ > _n_) memset_((_pointer_)+_n_,0,(_num_-_n_));                                   \
   } else {                                                                                   \
   newarr_(_pointer_,_num_);                                                                  \
   }                                                                                          \
} while(0)

#define incrarr_(_pointer_) do {            \
   if (_pointer_) {                         \
   int __num_ = (size_(_pointer_)+1);       \
   resizearr_(_pointer_,__num_);            \
   } else {                                 \
       newarr_(_pointer_,1);                \
   }                                        \
} while(0)

#define incrarr_n_(_pointer_,n) do {        \
   if (_pointer_) {                         \
   int __num_ = (size_(_pointer_)+n);       \
   resizearr_(_pointer_,__num_);            \
   } else {                                 \
       newarr_(_pointer_,n);                \
   }                                        \
} while(0)



#define appendarr_(_pointer_,_value_) do {         \
   incrarr_((_pointer_));                          \
   _pointer_[size_(_pointer_)-1] = _value_;        \
} while(0)


#define append_n_arr_(_pointer1_,_pointer2_) do {                                               \
   incrarr_n_((_pointer1_),size_(_pointer2_));                                                  \
   memcpy_((_pointer1_+size_(_pointer1_)-size_(_pointer2_)),_pointer2_,size_(_pointer2_));      \
} while(0)


#define cpy_pntr_to_arr_(_dest_,_source_,_n_) do{             \
   newarr_(_dest_,(_n_));                                     \
   memcpy_(_dest_,_source_,(_n_));                            \
} while(0)

#define free2darr_(pointer) do{                 \
   if (pointer)                                 \
   {                                            \
   int _i_;                                     \
   for(_i_=0;_i_ < size_(pointer);_i_++)        \
      freearr_(pointer[_i_]);                   \
   freearr_(pointer);                           \
   }                                            \
   /* pointer = 0;*/                            \
} while(0)

#define free3darr_(pointer) do{                 \
   if (pointer)                                 \
   {                                            \
   int _i_;                                     \
   for(_i_=0;_i_ < size_(pointer);_i_++)        \
      free2darr_(pointer[_i_]);                 \
   freearr_(pointer);                           \
   }                                            \
   /* pointer = 0;*/                            \
} while(0)

#define free4darr_(pointer) do{                 \
   if (pointer)                                 \
   {                                            \
   int _i_;                                     \
   for(_i_=0;_i_ < size_(pointer);_i_++)        \
      free3darr_(pointer[_i_]);                 \
   freearr_(pointer);                           \
   }                                            \
   /* pointer = 0;*/                            \
} while(0)

#define free5darr_(pointer) do{                 \
   if (pointer)                                 \
   {                                            \
   int _i_;                                     \
   for(_i_=0;_i_ < size_(pointer);_i_++)        \
      free4darr_(pointer[_i_]);                 \
   freearr_(pointer);                           \
   }                                            \
   /* pointer = 0;*/                            \
} while(0)




#define conv_pntr_to_arr_(_dest_,_pntr_,_n_) do{              \
   cpy_pntr_to_arr_(_dest_,_pntr_,_n_);                       \
   free(_pntr_);                                              \
} while(0)

#define memcpyarr_(_dest_,_source_) cpy_pntr_to_arr_(_dest_,_source_,size_(_source_))

#define setarr_(dest,int_const) memset_((dest),(int_const),(size_(dest)))

#define cp2darr_(_dest_,_source_) do {                   \
   if (_source_)                                         \
   {                                                     \
       int _i_;                                          \
       newarr_(_dest_,size_(_source_));                  \
       for(_i_=0; _i_ < size_(_source_);_i_++)           \
       {                                                 \
          if (_source_[_i_])                             \
          {                                              \
              memcpyarr_(_dest_[_i_],_source_[_i_]);     \
          }                                              \
       }                                                 \
   }                                                     \
} while(0)

#define minarr_(_arr_,_index_) do {                     \
   int _i_;                                             \
   _index_=0;                                           \
   typeof (_arr_[0]) _min_ = _arr_[0];                  \
   for(_i_=1; _i_ < size_(_arr_);_i_++)                 \
   {                                                    \
      if (_min_ > _arr_[_i_])                           \
      {                                                 \
         _min_ = _arr_[_i_];                            \
         _index_ = _i_;                                 \
      }                                                 \
   }                                                    \
} while(0)

#define avearr_(_arr_, _ave_) do {                      \
   int _i_;                                             \
   _ave_=0.0;                                           \
   for(_i_=0; _i_ < size_(_arr_);_i_++)                 \
   {                                                    \
      _ave_ += _arr_[_i_];                              \
   }                                                    \
   _ave_ /= ((double)size_(_arr_));                     \
} while(0)

#define sumarr_(_arr_, _sum_) do {                      \
   int _i_;                                             \
   _sum_=0.0;                                           \
   for(_i_=0; _i_ < size_(_arr_);_i_++)                 \
   {                                                    \
      _sum_ += _arr_[_i_];                              \
   }                                                    \
} while(0)

#define sum2darr_(_arr_, _sum_) do {                    \
   int _i_;                                             \
   _sum_=0.0;                                           \
   for(_i_=0; _i_ < size_(_arr_);_i_++)                 \
   {                                                    \
      _sum_ += sumarr_(_arr_[_i_]);                     \
   }                                                    \
} while(0)



#define MIN(a,b) ((a) < (b)?(a):(b))
#define MAX(a,b) ((a) > (b)?(a):(b))

int * consolidate_int(int * array_1,int n1,int *  array_2,int  n2, int * output_length);
int * shrink_arr_int(int *arry, int indx);//remove elt at indx
void * shrinkarr(void *arry, int indx, size_t size);//remove elt at indx
int sizearr(void * arr);//for gdb
int all_ones(int * arr);
int all_zeros(int * arr);
int * d2i_arr(double * arr);
double avearr(double * arr);
double sumarr(double * arr);
double sum2darr(double ** arr);
void printarr_(FILE * file, double * arr);
void printarr(double * arr);
void print2darr_(FILE * file, double ** arr);
void print2darr(double ** arr);


#endif// DEFS_H
