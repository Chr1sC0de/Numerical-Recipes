#pragma once

#include <fstream>
#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <fcntl.h>
#include <string.h>
#include <ctype.h>

namespace NumericalRecipes {

    //---------------------------------------------------------------------------------------s
    template<class T>
    inline T SQR(const T a) {return a*a;}
    //---------------------------------------------------------------------------------------
    template<class T>
    inline const T &MAX(const T &a, const T &b){
        return b > a ? (b) : (a);
    }
    //---------------------------------------------------------------------------------------
    inline float MAX(const double &a, const float &b){
        return b > a ? (b) : float(a);
    }
    //---------------------------------------------------------------------------------------
    inline float MAX(const float &a, const double &b){
        return b > a ? float(b) : (a);
    }
    //---------------------------------------------------------------------------------------
    template<class T>
    inline const T &MIN(const T &a, const T &b){
        return b < a ? (b) : (a);
    }
    //---------------------------------------------------------------------------------------
    inline float MIN(const double &a, const float &b){
        return b < a ? (b) : float(a);
    }
    //---------------------------------------------------------------------------------------
    inline float MIN(const float &a, const double &b){
        return b < a ? float(b) : (a);
    }
    //---------------------------------------------------------------------------------------
    template<class T>
    inline T SIGN(const T &a, const T &b){
        return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
    }
    //---------------------------------------------------------------------------------------
    inline float SIGN(const float &a, const double &b){
        return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
    }
    //---------------------------------------------------------------------------------------
    inline float SIGN(const double &a, const float &b){
        return (float)(b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a));
    }
    //---------------------------------------------------------------------------------------
    template<class T>
    inline void SWAP(T &a, T &b){
        T dum=a; a=b; b=dum;
    }
    //---------------------------------------------------------------------------------------
    template <class T>
    class Vector_t {
        private:
            // Size of array, indices 0..nn-1.
            int nn;
            // Pointer to data array.
            T *v;
        public:
            // Make T available.
            typedef T value_type;
            // Default constructor.
            Vector_t(): nn(0), v(NULL){}
            // Construct vector of size n.
            explicit Vector_t(int n): nn(n), v(n>0 ? new T[n]: NULL){}
            // Initialize to constant value a.
            Vector_t(int n, const T &a)
                : nn(n), v(n>0 ? new T[n]: NULL)
            {
                for(int i=0; i<n; i++){v[i] = a;}
            }
            // Initialize to values in C-style array a.
            Vector_t(int n, const T *a)
            :
                nn(n), v(n>0 ? new T[n]: NULL)
            {
                for(int i=0; i<n; i++){v[i] = *a;}
            }
            // Copy constructor.
            Vector_t(const Vector_t<T> &rhs)
            :
                nn(rhs.nn), v(nn>0 ? new T[nn]: NULL)
            {
                for (int i=0; i<nn; i++)
                    v[i] = rhs[i];
            }
            // Assignment operator.
            Vector_t & operator=(const Vector_t<T> &rhs){
                if (this != &rhs)
                {
                    if (nn != rhs.nn) {
                        if (v != NULL){
                            delete [] (v);
                        }
                        nn = rhs.nn;
                        v = nn>0 ? new T[nn] : NULL;
                    }
                    for (int i=0; i<nn; i++)
                        v[i]=rhs[i];
                }
                return *this;
            }
            // Return element number i.
            inline T & operator[](const int i){
                #ifdef _CHECKBOUNDS_
                    if (i<0 || i>=nn) {
                        throw("NRvector subscript out of bounds");
                    }
                #endif
                    return v[i];
            }
            // const version.
            inline const T & operator[](const int i) const{
                #ifdef _CHECKBOUNDS_
                    if (i<0 || i>=nn) {
                        throw("NRvector subscript out of bounds");
                    }
                #endif
                    return v[i];
            }
            // Return size of vector.
            inline int size() const {return nn;}
            // Resize, losing contents.
            void resize(int newn){
                if (newn != nn) {
                    if (v != NULL) delete[] (v);
                    nn = newn;
                    v = nn > 0 ? new T[nn] : NULL;
                }
            }
            // Resize and assign a to every element.
            void assign(int newn, const T &a){
                if (newn != nn) {
                    if (v != NULL) delete[] (v);
                    nn = newn;
                    v = nn > 0 ? new T[nn] : NULL;
                }
                for (int i=0;i<nn;i++) v[i] = a;
            }
            // Destructor.
            ~Vector_t(){
                if (v != NULL) delete[] (v);
            }
    };
    //---------------------------------------------------------------------------------------
    template <class T>
    class Matrix2d_t {
        private:
            // Number of rows and columns. Index range is 0..nn-1, 0..mm-1.
            int nn;
            int mm;
            // Storage for data.
            T **v;
        public:
            // Make T available.
            typedef T value_type;

            // Default constructor.
            Matrix2d_t(): nn(0), mm(0), v(NULL) {}
            // Construct n x m matrix.
            Matrix2d_t(int n, int m)
            :
                nn(n), mm(m),
                v(n>0 ? new T*[n] : NULL)
            {
                int i,nel=m*n;
                if (v)
                   v[0] = nel>0 ? new T[nel] : NULL;
                for (i=1;i<n;i++)
                    v[i] = v[i-1] + m;
            }
            // Initialize to constant value a.
            Matrix2d_t(int n, int m, const T &a)
            :
                nn(n), mm(m),
                v(n>0 ? new T*[n] : NULL)
            {
                int i,j,nel=m*n;
                if (v)
                    v[0] = nel>0 ? new T[nel] : NULL;
                for (i=1; i< n; i++)
                    v[i] = v[i-1] + m;
                for (i=0; i< n; i++)
                    for (j=0; j<m; j++)
                        v[i][j] = a;
            }
            // Initialize to values in C-style array a.
            Matrix2d_t(int n, int m, const T *a)
            :
                nn(n), mm(m),
                v(n>0 ? new T*[n] : NULL)
            {
                int i,j,nel=m*n;
                if
                    (v) v[0] = nel>0 ? new T[nel] : NULL;
                for (i=1; i< n; i++)
                    v[i] = v[i-1] + m;
                for (i=0; i< n; i++)
                    for (j=0; j<m; j++)
                        v[i][j] = *a++;
            }
            // Copy constructor.
            Matrix2d_t(const Matrix2d_t &rhs)
            :
                nn(rhs.nn), mm(rhs.mm),
                v(nn>0 ? new T*[nn] : NULL)
            {
                int i,j,nel=mm*nn;
                if (v)
                    v[0] = nel>0 ? new T[nel] : NULL;
                for (i=1; i< nn; i++)
                    v[i] = v[i-1] + mm;
                for (i=0; i< nn; i++)
                    for (j=0; j<mm; j++)
                        v[i][j] = rhs[i][j];
            }
            // Assignment operator.
            Matrix2d_t & operator=(const Matrix2d_t &rhs){
                if (this != &rhs) {
                    int i,j,nel;
                    if (nn != rhs.nn || mm != rhs.mm) {
                        if (v != NULL) {
                            delete[] (v[0]);
                            delete[] (v);
                        }
                        nn=rhs.nn;
                        mm=rhs.mm;
                        v = nn>0 ? new T*[nn] : NULL;
                        nel = mm*nn;
                        if (v) v[0] = nel>0 ? new T[nel] : NULL;
                        for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
                    }
                    for (i=0; i< nn; i++) for (j=0; j<mm; j++) v[i][j] = rhs[i][j];
                }
                return *this;
            }
            // Subscripting: pointer to row i.
            inline T* operator[](const int i){
                #ifdef _CHECKBOUNDS_
                if (i<0 || i>=nn) {
                    throw("NRmatrix subscript out of bounds");
                }
                #endif
                    return v[i];
            }
            // const version.
            inline const T* operator[](const int i) const{
                #ifdef _CHECKBOUNDS_
                if (i<0 || i>=nn) {
                    throw("NRmatrix subscript out of bounds");
                }
                #endif
                    return v[i];
            }
            // Return number of rows.
            inline int nrows() const {
                return nn;
            }
            // Return number of columns.
            inline int ncols() const{
                return mm;
            }
            // Resize, losing contents.
            void resize(int newn, int newm){
                int i,nel;
                if (newn != nn || newm != mm) {
                    if (v != NULL) {
                        delete[] (v[0]);
                        delete[] (v);
                    }
                    nn = newn;
                    mm = newm;
                    v = nn>0 ? new T*[nn] : NULL;
                    nel = mm*nn;
                    if (v) v[0] = nel>0 ? new T[nel] : NULL;
                    for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
                }
            }
            // Resize and assign a to every element.
            void assign(int newn, int newm, const T &a){
                int i,j,nel;
                if (newn != nn || newm != mm) {
                    if (v != NULL) {
                        delete[] (v[0]);
                        delete[] (v);
                    }
                    nn = newn;
                    mm = newm;
                    v = nn>0 ? new T*[nn] : NULL;
                    nel = mm*nn;
                    if (v) v[0] = nel>0 ? new T[nel] : NULL;
                    for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
                }
                for (i=0; i< nn; i++) for (j=0; j<mm; j++) v[i][j] = a;
            }
            // Destructor.
            ~Matrix2d_t(){
                if (v != NULL) {
                    delete[] (v[0]);
                    delete[] (v);
                }
            }
    };


}