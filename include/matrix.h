#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "ofts.h"


/**
 * \file matrix.h
 * \brief An extension of the vector class of C++ for 2-dim matrix manipulation.
 *        See the .tpp file for additional comments and briefs of the routines.
 * \author BLB
 */


//----------------------------------------------------------------------------------------
// matrix<T> class
//----------------------------------------------------------------------------------------
template<typename T>
class matrix;



template <typename T>
class matrix
{

private:

    int m_size_r;      //rows
    int m_size_c;      //columns
    vector<T> m_coef;  //coefficients

public:

    //------------------------------------------------------------------------------------
    //Create
    //------------------------------------------------------------------------------------
    /**
     *  \brief Default creator for matrix class. Never use.
     **/
    matrix<T>();
    /**
     *  \brief Default creator for matrix class, with size m_size_r x m_size_c.
     **/
    matrix<T>(const int t_size_r, const int t_size_c);

    //------------------------------------------------------------------------------------
    //Copy
    //------------------------------------------------------------------------------------
    /**
     *  \brief Create a matrix equal to the matrix b. Requires the routines:
     *          - ccopy
     **/
    matrix<T>(matrix const& b);
    /**
     *  \brief Copy from a given matrix object (only the coefficients).
     */
    matrix<T>& ccopy(matrix<T> const& b);
    /**
     *  \brief Linked copy from a given matrix object (exact same object is obtained).
     */
    matrix<T>& lcopy(matrix<T> const& b);

    //------------------------------------------------------------------------------------
    //Delete
    //------------------------------------------------------------------------------------
    /**
     *  \brief Default destructor.
     */
    ~matrix<T>();

    //------------------------------------------------------------------------------------
    //Setters
    //------------------------------------------------------------------------------------
    /**
     *  \brief Set a given coefficient to the position (i, j) in the matrix<T>
     */
    void set_coef(T const & value, int i, int j);
    /**
     *  \brief Add a given coefficient to the position (i, j) in the matrix<T>
     */
    void add_coef(T const & value, int i, int j);
    /**
     *  \brief Set a given coefficient of type <U> to the position (i, j) in the matrix<T>
     */
    template <typename U> void set_coef(U const & value, int i, int j);
    /**
     *  \brief Set a given coefficient of type <cdouble> to the position (i, j) in the matrix<T>
     */
    void set_coef(cdouble const & value, int i, int j);
    /**
     *  \brief Zeroing of the matrix. Requires zero().
     **/
    void zero();
    //------------------------------------------------------------------------------------
    //Getters
    //------------------------------------------------------------------------------------
    /**
     *  \brief Get the const reference to the coefficient (i,j)
     */
    const T& get_ref_coef(int i, int j) const;
    /**
     *  \brief Get the coefficient (i,j)
     */
    T get_coef(int i, int j) const;
    /**
     *  \brief Get a ptr to the coefficient (i,j)
     */
    T* get_ptr_first_coef(int i, int j) const;
    /**
     *  \brief Get the size (num = 1 for m_size_r, num = 2 for m_size_c).
     */
    int get_size(int num) const;

    //------------------------------------------------------------------------------------
    //Operators
    //------------------------------------------------------------------------------------
    /**
     *  \brief Create a matrix equal to the matrix b. Requires the routines:
     *          - ccopy
     **/
    matrix<T>& operator  = (matrix<T> const& b);

    //------------------------------------------------------------------------------------
    //Operations
    //------------------------------------------------------------------------------------
    /**
     *  \brief Derivation of a wrt to the variable ni, set at position (i,j) in the matrix. Requires der.
     **/
    void der(T const &a, int ni, int i, int j);
    /**
     *  \brief Derivation of a wrt to the variable ni, set at position (i,j) in the matrix,
     *         at order k of the coefficients. Requires der (at order k).
     **/
    void der(T const &a, int ni, int i, int j, int k);
    /**
     *  \brief Derivation of a wrt to time, set at position (i,j) in the matrix. Frequency is n. Requires dot.
     **/
    void dot(T const &a, double n, int i, int j);
    /**
     *  \brief Derivation of a wrt to time of the whole matrix. Frequency is n. Requires dot.
     **/
    void dot(matrix<T> const &a, double n);
    /**
     *  \brief Derivation of a wrt to the variable ni, set at position (i,j) in the matrix,
     *         at order k of the coefficients. Requires tfts_der (at order k).
     **/
    void tfts_der(T const &a, int ni, int i, int j, int k);

    //------------------------------------------------------------------------------------
    //Operators
    //------------------------------------------------------------------------------------
    /**
     *  \brief Operator, to be able to use matrix[i].
     **/
    T operator [](int i) const    {return m_coef[i];}
    /**
     *  \brief Operator, to be able to use matrix[i].
     **/
    T & operator [](int i) {return m_coef[i];}
    /**
     *  \brief Operator, to be able to use matrix[i][j].
     **/
    T operator ()(int i, int j) const    {return m_coef[m_size_c*i+j];}
    /**
     *  \brief Operator, to be able to use matrix[i][j].
     **/
    T & operator ()(int i, int j) {return m_coef[m_size_c*i+j];}
};

//----------------------------------------------------------------------------------------
//Functions used with T = Ofts< Ofsc >, OFS format by default (see .tpp for brief)
//----------------------------------------------------------------------------------------
template <typename T>  void smvprod_t(matrix<T> const& a, vector<T> const& vIn, vector<T>& vOut);
template <typename T>  void smvprod_t(matrix<T> const& a, vector<T> const& vIn, vector<T>& vOut, int const& m);
template <typename T , typename U> void add_coef(vector<U> const& a, vector<T>& vOut);
template <typename T , typename U> void sub_coef(vector<U> const& a, vector<T>& vOut);

//----------------------------------------------------------------------------------------
//Functions used with T = Ofts< Ofsc >, TFS format (see .tpp for brief)
//----------------------------------------------------------------------------------------
template <typename T , typename U> void tfts_smvprod_u(matrix<U> const& a, vector<T> const& vIn, vector<T>& vOut, int constm);
template <typename T>  void tfts_smvprod_t(matrix<T> const& a, vector<T> const& vIn, vector<T>& vOut, int constm);
template <typename T , typename U> void tfts_subCoef(vector<U> const& a, vector<T>& vOut);

//----------------------------------------------------------------------------------------
//TFS <--> OFS format (see .tpp for brief)
//----------------------------------------------------------------------------------------
inline void tfs_from_ofs_inline(matrix< Ofsc >& a);
inline void tfs_to_ofs_inline(matrix< Ofsc >& a);
inline void tfs_from_ofs_inline(vector< Ofts< Ofsc > >& a, int m);
inline void tfs_to_ofs_inline(vector< Ofts< Ofsc > >& a, int m);
inline void tfs_to_ofs_inline(matrix< Ofts< Ofsc > >& a);

//----------------------------------------------------------------------------------------
//Functions used with T = Ofts< Ofsc >, to use only in tests (see .tpp for brief)
//----------------------------------------------------------------------------------------
template <typename T , typename U> void add_coef(vector<U> const& a, vector<T> const& vIn, vector<T>& vOut);
template <typename T , typename U> void sub_coef(vector<U> const& a, vector<T> const& vIn, vector<T>& vOut);


//----------------------------------------------------------------------------------------
//Functions used with T = Ofs<U> (see .tpp for brief)
//----------------------------------------------------------------------------------------
void smvprod_ofs(matrix<Ofsc> const& a, vector<Ofsc> const& vIn, vector<Ofsc>& vOut);

//----------------------------------------------------------------------------------------
//Functions used with U = cdouble (see .tpp for brief)
//----------------------------------------------------------------------------------------
template <typename U> void mvprod_u(matrix< Ofs<U> > const& a, vector<U> const& vIn, vector<U>& vOut, double const& t);
template <typename U> void smvprod_u(matrix< Ofs<U> > const& a, vector<U> const& vIn, vector<U>& vOut, double const& t);
template <typename U> void vvsum_u(vector< Ofs<U> > & a, vector<U>& vOut, double const& t);
template <typename U> void vvsub_u(vector< Ofs<U> > & a, vector<U>& vOut, double const& t);
template <typename U> void vvsub_u(vector< Ofs<U> > & a, vector<U> const& vIn, vector<U>& vOut, double const& t);
template <typename U> void smmprod_u(matrix< Ofs<U> > const& a, matrix< Ofs<U> > const& b, matrix< Ofs<U> >& mOut);
template <typename U> void smtmprod_u(matrix< Ofs<U> > const& a, matrix< Ofs<U> > const& b, matrix< Ofs<U> >& mOut);

//----------------------------------------------------------------------------------------
// Read & Write
//----------------------------------------------------------------------------------------
/**
 * \brief Writes a given matrix W of type \c Ofts<Ofsc >  in a binary files of the form "filename+i+j.bin", with i = 0, m_size_r(W)-1
 *        and j = 0, m_size_c(W)-1.
 **/
void write_mofts_bin(matrix<Ofts<Ofsc > > &W, string filename);

/**
 * \brief Reads a given matrix W of type \c Ofts<Ofsc >  in a binary files of the form "filename+i+j.bin", with i = 0, m_size_r(W)-1
 *        and j = 0, m_size_c(W)-1.
 **/
inline void read_mofts_bin(matrix<Ofts<Ofsc > > &W, string filename);

//----------------------------------------------------------------------------------------
//Include the implementation .tpp
//----------------------------------------------------------------------------------------
#include "matrix.tpp"

//########################################################################################
// end of matrix<T> class
//########################################################################################

#endif // MATRIX_H_INCLUDED
