//########################################################################################
// Implementation of the matrix template class
//########################################################################################

/**
 * \file matrix.tpp
 * \brief Some extension of the vector class of C++ for matrix manipulation
 * \author BLB
 * \date May 2015
 * \version 1.0
 */

//----------------------------------------------------------------------------------------
//Create
//----------------------------------------------------------------------------------------
/**
 *  \brief Default creator for matrix class. Never use.
 **/
template <typename T> matrix<T>::matrix()
{
    m_size_r = 0;//NV;
    m_size_c = 0;//REDUCED_NV;
    m_coef = *(new vector<T>());
}

/**
 *  \brief Default creator for matrix class, with size m_size_r x m_size_c.
 **/
template <typename T> matrix<T>::matrix(const int t_size_r, const int t_size_c)
{
    m_size_r = t_size_r;
    m_size_c = t_size_c;
    m_coef = *(new vector<T>(m_size_r*m_size_c));
}

//----------------------------------------------------------------------------------------
//Copy
//----------------------------------------------------------------------------------------
/**
 *  \brief Create a matrix equal to the matrix b. Requires the routines:
 *          - ccopy
 **/
template <typename T> matrix<T>::matrix(matrix const& b)
{
    m_size_r = b.m_size_r;
    m_size_c = b.m_size_c;
    m_coef  = *(new vector<T>(m_size_r*m_size_c));
    for(int i = 0 ; i< m_size_r*m_size_c; i++) m_coef[i].ccopy(b.m_coef[i]);
}

/**
 *  \brief Create a matrix equal to the matrix b. Requires the routine:
 *          - ccopy
 **/
template <typename T> matrix<T>& matrix<T>::operator = (matrix<T> const& b)
{
    if(this != &b)
    {
        m_size_r = b.m_size_r;
        m_size_c = b.m_size_c;
        m_coef  = *(new vector<T>(m_size_r*m_size_c));
        for(int i = 0 ; i < m_size_r*m_size_c; i++) m_coef[i].ccopy(b.m_coef[i]);
    }
    return *this; //same object if returned
}

/**
 *  \brief Copy the matrix b in this. Requires the routine:
 *          - ccopy
 **/
template <typename T> matrix<T>& matrix<T>::ccopy(matrix<T> const& b)
{
    if(m_size_r != b.m_size_r || m_size_c != b.m_size_c)
    {
        cout << "Erreur in ccopy for matrix: sizes do not match. Nothing is done." << endl;
        return *this;
    }
    else
    {

    for(int i = 0 ; i < m_size_r*m_size_c; i++) m_coef[i].ccopy(b.m_coef[i]);
    return *this; //same object if returned
    }
}

/**
 *  \brief Link the matrix b with this. Requires the routine:
 *          - lcopy
 **/
template <typename T> matrix<T>& matrix<T>::lcopy(matrix<T> const& b)
{
    m_size_r = b.m_size_r;
    m_size_c = b.m_size_c;
    for(int i = 0 ; i < m_size_r*m_size_c; i++) m_coef[i].lcopy(b.m_coef[i]);
    return *this;
}

//----------------------------------------------------------------------------------------
//Delete
//----------------------------------------------------------------------------------------
/**
 *  \brief Default destructor.
 *         Do we need to implement something here? Probably not, not pointer used.
 **/
template <typename T> matrix<T>::~matrix<T>()
{
}

//----------------------------------------------------------------------------------------
//Getters
//----------------------------------------------------------------------------------------
/**
 *  \brief Get the const reference to the coefficient (i,j)
 */
template <typename T> const T& matrix<T>::get_ref_coef(int i, int j) const
{
    if( i >= m_size_r || j >= m_size_c)
    {
        cout << "Error in matrix<T>::get_coef: indices are out of scope. First coefficient is returned." << endl;
        return m_coef[0];
    }
    else return m_coef[i*m_size_c + j];
}

/**
 *  \brief Gets the coefficient at the position (i,j)
 **/
template <typename T> T matrix<T>::get_coef(int i, int j) const
{
    if( i >= m_size_r || j >= m_size_c)
    {
        cout << "Error in matrix<T>::get_coef: indices are out of scope. First coefficient is returned." << endl;
        return m_coef[0];
    }
    else return m_coef[i*m_size_c + j];
}


/**
 *  \brief Get a ptr to the coefficient (i,j)
 */
template <typename T> T* matrix<T>::get_ptr_first_coef(int i, int j) const
{
    if( i >= m_size_r || j >=m_size_c)
    {
        cout << "Error in matrix<T>::get_coef: indices are out of scope. First coefficient is returned." << endl;
        return m_coef[0].get_ptr();
    }
    else return m_coef[i*m_size_c + j].get_ptr();
}

/**
 *  \brief Get the size (num = 1 for m_size_r, num = 2 for m_size_c).
 */
template <typename T> int matrix<T>::get_size(int num) const
{
    if(num == 1) return m_size_r;
    else if(num == 2) return m_size_c;
    else cout << "Error in matrix<T>::get_size: required number must be 1 or 2." << endl; return 0;

}

//----------------------------------------------------------------------------------------
//Setters
//----------------------------------------------------------------------------------------
/**
 *  \brief Sets the coefficient T at the position (i,j). Requires ccopy.
 **/
template <typename T> void matrix<T>::set_coef(T const & value, int i, int j)
{
    this->get_ptr_first_coef(i,j)->ccopy(value);
}

/**
 *  \brief Sets the subcoefficient value at the order zero of the coefficient (i,j). Requires set_coef.
 **/
template <typename T> template <typename U> void matrix<T>::set_coef(U const & value, int i, int j)
{
    this->get_ptr_first_coef(i,j)->set_coef(value, (int const) 0);
}

/**
 *  \brief Adds the coefficient T at the position (i,j). Requires smult.
 **/
template <typename T> void matrix<T>::add_coef(T const & value, int i, int j)
{
   this->get_coef(i,j).smult(value, 1.0);
}

/**
 *  \brief Zeroing of the matrix. Requires zero().
 **/
template <typename T> void matrix<T>::zero()
{
    for(int i = 0 ; i < m_size_r*m_size_c; i++) m_coef[i].zero();
}

//----------------------------------------------------------------------------------------
//Operations
//----------------------------------------------------------------------------------------
/**
 *  \brief Derivation of a wrt to the variable ni, set at position (i,j) in the matrix. Requires der.
 **/
template <typename T> void matrix<T>::der(T const &a, int ni, int i, int j)
{
    if(i >= m_size_r || j >= m_size_c)
    {
        cout << "Error in matrix<T>::der: indices are out of scope. Nothing is done." << endl;
    }
    else m_coef[i*m_size_c +j].der(a, ni);
}

/**
 *  \brief Derivation of a wrt to the variable ni, set at position (i,j) in the matrix, at order k of the coefficients. Requires der (at order k).
 **/
template <typename T> void matrix<T>::der(T const &a, int ni, int i, int j, int k)
{
    if(i >= m_size_r || j >= m_size_c)
    {
        cout << "Error in matrix<T>::der: indices are out of scope. Nothing is done." << endl;
    }
    else m_coef[i*m_size_c +j].der(a, ni, k);
}

/**
 *  \brief Derivation of a wrt to time, set at position (i,j) in the matrix. Frequency is n. Requires dot.
 **/
template <typename T> void matrix<T>::dot(T const &a, double n, int i, int j)
{
    if(i >= m_size_r || j >= m_size_c)
    {
        cout << "Error in matrix<T>::dot: indices are out of scope. Nothing is done." << endl;
    }
    else get_ptr_first_coef(i,j)->dot(a, n);
}

/**
 *  \brief Derivation of a wrt to time of the whole matrix. Frequency is n. Requires dot.
 **/
template <typename T> void matrix<T>::dot(matrix<T> const &a, double n)
{
    if(a.m_size_r != m_size_r || a.m_size_c != m_size_c)
    {
        cout << "Error in matrix<T>::dot: indices are out of scope. Nothing is done." << endl;
    }
    else
    {
        for(int i = 0; i< m_size_r; i++)
            for(int j = 0; j< m_size_c; j++) this->dot(a.get_coef(i,j), n, i, j);
    }
}

//----------------------------------------------------------------------------------------
//Functions used with T = Ofts< Ofsc >
//----------------------------------------------------------------------------------------
/**
 *  \brief Matrix-vector product: vOut += a x vIn. Used with T = Ofts< Ofsc >. Requires ofts_sprod.
 **/
template <typename T>  void smvprod_t(matrix<T> const& a, vector<T> const& vIn, vector<T>& vOut)
{
    if(a.get_size(2) != vIn.size() || a.get_size(1) != vOut.size() )
    {
        cout << "Error in mvprod (matrix.tpp): lengths do not match. Nothing is done." << endl;
    }
    else
    {
        for(int i =0; i < a.get_size(1) ; i++)
        {
            for(int j =0; j < a.get_size(2); j++)
            {
                ( (T*)&vOut[i])->ofts_sprod(a.get_coef(i,j), vIn[j]);
            }
        }
    }
}

/**
 *  \brief Matrix-vector product at order m: vOut += a x vIn. Used with T = Ofts< Ofsc >. Requires ofts_sprod.
 **/
template <typename T>  void smvprod_t(matrix<T> const& a, vector<T> const& vIn, vector<T>& vOut, int const m)
{
    if( a.get_size(2) != (int) vIn.size() || a.get_size(1) != (int) vOut.size() )
    {
        cout << "Error in mvprod (matrix.tpp): lengths do not match. Nothing is done." << endl;
    }
    else
    {
        for(int i =0; i < a.get_size(1) ; i++)
        {
            for(int j =0; j < a.get_size(2); j++)
            {
                vOut[i].ofts_sprod(*a.get_ptr_first_coef(i,j), vIn[j], m);
            }
        }
    }
}

/**
 *  \brief Matrix-vector product: vOut += a x vIn. Used with U = Ofsc, T = Ofts< Ofsc >. Requires ofts_smult_t.
 **/
template <typename T , typename U> void smvprod_u(matrix<U> const& a, vector<T> const& vIn, vector<T>& vOut)
{
    if( a.get_size(2) != (int) vIn.size() || a.get_size(1) != (int) vOut.size() )
    {
        cout << "Error in mvprod (matrix.tpp): lengths do not match. Nothing is done." << endl;
    }
    else
    {
        for(int i =0; i < a.get_size(1) ; i++)
        {
            for(int j =0; j < a.get_size(2); j++)
            {
                ( (T*)&vOut[i])->ofts_smult_t(vIn[j], (U&) a.get_ref_coef(i,j));//m_coef[i*a.m_size_c + j]);
            }
        }
    }
}

/**
 *  \brief Matrix-vector product at order m: vOut += a x vIn. Used with U = Ofsc, T = Ofts< Ofsc >. Requires ofts_smult_t.
 **/
template <typename T , typename U> void smvprod_u(matrix<U> const& a, vector<T> const& vIn, vector<T>& vOut, int const m)
{
    if( a.get_size(2) != (int) vIn.size() || a.get_size(1) != (int) vOut.size() )
    {
        cout << "Error in mvprod (matrix.tpp): lengths do not match. Nothing is done." << endl;
    }
    else
    {
        for(int i =0; i < a.get_size(1) ; i++)
        {
            for(int j =0; j < a.get_size(2); j++)
            {
                ( (T*)&vOut[i])->ofts_smult_t(vIn[j], (U&) a.get_ref_coef(i,j), m);
            }
        }
    }
}

/**
 *  \brief Add the algebraic vector a to the order zero of vOut. Used with T = Ofts< Ofsc >. Requires add_coef from vector class.
 **/
template <typename T , typename U> void add_coef(vector<U> const& a, vector<T>& vOut)
{
    if(a.size() != vOut.size() )
    {
        cout << "Error in add_coef (matrix.tpp): lengths do not match. Nothing is done." << endl;
    }
    else
    {
        for(int i =0; i < (int) vOut.size() ; i++)
        {
                ( (T*)&vOut[i])->add_coef(a[i], 0, 0);
        }
    }
}

/**
 *  \brief Sub the algebraic vector a to the order zero of vOut. Used with T = Ofts< Ofsc >. Requires add_coef from vector class.
 **/
template <typename T , typename U> void sub_coef(vector<U> const& a, vector<T>& vOut)
{
    if(a.size() != vOut.size() )
    {
        cout << "Error in add_coef (matrix.tpp): lengths do not match. Nothing is done." << endl;
    }
    else
    {
        for(int i =0; i < (int) vOut.size() ; i++)
        {
                ( (T*)&vOut[i])->add_coef(-a[i], 0, 0);
        }
    }
}


//----------------------------------------------------------------------------------------
//Functions used with T = Ofts< Ofsc >, TFS format
//----------------------------------------------------------------------------------------
/**
 *  \brief Matrix-vector product at order m: vOut += a x vIn. Used with U = Ofsc, T = Ofts< Ofsc >, TFS format. Requires ofts_smult_t.
 **/
template <typename T , typename U> void tfts_smvprod_u(matrix<U> const& a, vector<T> const& vIn, vector<T>& vOut, int const m)
{
    if( a.get_size(2) != (int) vIn.size() || a.get_size(1) != (int) vOut.size() )
    {
        cout << "Error in mvprod (matrix.tpp): lengths do not match. Nothing is done." << endl;
    }
    else
    {
        for(int i =0; i < a.get_size(1) ; i++)
        {
            for(int j =0; j < a.get_size(2); j++)
            {
                ( (T*)&vOut[i])->tfts_smult_t(vIn[j], a.get_ref_coef(i,j), m);
            }
        }
    }
}

/**
 *  \brief Matrix-vector product at order m: vOut += a x vIn. Used with T = Ofts< Ofsc >, TFS format. Requires ofts_sprod.
 **/
template <typename T>  void tfts_smvprod_t(matrix<T> const& a, vector<T> const& vIn, vector<T>& vOut, int const m)
{
    if( a.get_size(2) != (int) vIn.size() || a.get_size(1) != (int) vOut.size() )
    {
        cout << "Error in mvprod (matrix.tpp): lengths do not match. Nothing is done." << endl;
    }
    else
    {
        for(int i =0; i < a.get_size(1) ; i++)
        {
            for(int j =0; j < a.get_size(2); j++)
            {
                vOut[i].tfts_sprod(*a.get_ptr_first_coef(i,j), vIn[j], m);
            }
        }
    }
}

/**
 *  \brief Sub the algebraic vector a to the order zero of vOut. Used with T = Ofts< Ofsc >. Requires add_coef from vector class.
 *         Note the routine version for TFS format is identical to the one for OFS format because the underlying routines are used the same way.
 **/
template <typename T , typename U> void tfts_subCoef(vector<U> const& a, vector<T>& vOut)
{
    if(a.size() != vOut.size() )
    {
        cout << "Error in add_coef (matrix.tpp): lengths do not match. Nothing is done." << endl;
    }
    else
    {
        for(int i =0; i < (int) vOut.size() ; i++)
        {
                ( (T*)&vOut[i])->add_coef(-a[i], 0, 0);
        }
    }
}

/**
 *  \brief Add the algebraic vector a to the order zero of vOut. Used with T = Ofts< Ofsc >. Requires add_coef from vector class.
 *         Note the routine version for TFS format is identical to the one for OFS format because the underlying routines are used the same way.
 **/
template <typename T , typename U> void tfts_add_coef(vector<U> const& a, vector<T>& vOut)
{
    if(a.size() != vOut.size() )
    {
        cout << "Error in add_coef (matrix.tpp): lengths do not match. Nothing is done." << endl;
    }
    else
    {
        for(int i =0; i < (int) vOut.size() ; i++)
        {
                ( (T*)&vOut[i])->add_coef(a[i], 0, 0);
        }
    }
}

/**
 *  \brief Derivation of a wrt to the variable ni, set at position (i,j) in the matrix, at order k of the coefficients. Requires tfts_der (at order k).
 **/
template <typename T> void matrix<T>::tfts_der(T const &a, int ni, int i, int j, int k)
{
    if(i >= m_size_r || j >= m_size_c)
    {
        cout << "Error in matrix<T>::der: indices are out of scope. Nothing is done." << endl;
    }
    else m_coef[i*m_size_c +j].tfts_der(a, ni, k);
}


//----------------------------------------------------------------------------------------
//TFS <--> OFS format
//----------------------------------------------------------------------------------------
/**
 *  \brief  Inline from frequency domain to time domain, for matrix< Ofsc > object.
 */
inline void tfs_from_ofs_inline(matrix< Ofsc >& a)
{
    Ofsc temp(a.get_ptr_first_coef(0,0)->get_order());
    for(int i =0; i < a.get_size(1) ; i++)
    {
        for(int j =0; j < a.get_size(2); j++)
        {
            a.get_ptr_first_coef(i,j)->tfs_from_ofs_inline(temp);
        }
    }
}

/**
 *  \brief  Inline from time domain to frequency domain, for matrix< Ofsc > object.
 */
inline void tfs_to_ofs_inline(matrix< Ofsc >& a)
{
    for(int i =0; i < a.get_size(1) ; i++)
    {
        for(int j =0; j < a.get_size(2); j++)
        {
            a.get_ptr_first_coef(i,j)->tfs_to_ofs_inline();
        }
    }
}

/**
 *  \brief  Inline from frequency domain to time domain, for vector< Ofsc > object.
 */
inline void tfs_from_ofs_inline(vector< Ofsc >& a)
{
    Ofsc temp(a[0].get_order());
    for(int i =0; i < (int) a.size() ; i++)
    {
            a[i].tfs_from_ofs_inline(temp);
    }
}

/**
 *  \brief  Inline from time domain to frequency domain, for vector< Ofsc > object.
 */
inline void tfs_to_ofs_inline(vector< Ofsc >& a)
{
    for(int i =0; i < (int) a.size() ; i++)
    {
            a[i].tfs_to_ofs_inline();
    }
}

/**
 *  \brief  Inline from frequency domain to time domain, for vector< Ofts< Ofsc >  > object.
 */
inline void tfs_from_ofs_inline(vector< Ofts< Ofsc > >& a, int m)
{
    for(int i =0; i < (int) a.size() ; i++)
    {
            a[i].tfs_from_ofs_inline(m);
    }
}

/**
 *  \brief  Inline from time domain to frequency domain, for vector< Ofts< Ofsc > > object.
 */
inline void tfs_to_ofs_inline(vector< Ofts< Ofsc > >& a, int m)
{
    for(int i =0; i < (int) a.size() ; i++)
    {
            a[i].tfs_to_ofs_inline(m);
    }
}

/**
 *  \brief  Inline from time domain to frequency domain, for vector< Ofts< Ofsc > > object.
 */
inline void tfs_to_ofs_inline(vector< Ofts< Ofsc > >& a)
{
    for(int m = 0; m <= a[0].get_order(); m++)
    {
        for(int i =0; i < (int) a.size() ; i++)
        {
            a[i].tfs_to_ofs_inline(m);
        }
    }
}

/**
 *  \brief  Inline from time domain to frequency domain, for matrix< Ofts< Ofsc > > object.
 */
inline void tfs_to_ofs_inline(matrix< Ofts< Ofsc > >& a)
{
    for(int i =0; i < a.get_size(1) ; i++)
    {
        for(int j =0; j < a.get_size(2); j++)
        {
            a.get_ptr_first_coef(i,j)->tfs_to_ofs_inline();
        }
    }
}


//----------------------------------------------------------------------------------------
//Functions used with T = Ofts< Ofsc >, to use only in tests
//----------------------------------------------------------------------------------------
/**
 *  \brief Add the algebraic vector a to the order zero of vIn and set the result in vOut. Used with T = Ofts< Ofsc >. Requires add_coef from vector class.
 *         Careful: this routine leads to the copy of an entire Ofts structure. Very heavy! Use only in tests.
 **/
template <typename T , typename U> void add_coef(vector<U> const& a, vector<T> const& vIn, vector<T>& vOut)
{
    if(a.get_size() != vOut.size() )
    {
        cout << "Error in add_coef (matrix.tpp): lengths do not match. Nothing is done." << endl;
    }
    else
    {
        for(int i =0; i < vOut.size() ; i++)
        {
                ((T*)&vOut[i])->ccopy(vIn[i]);
                ((T*)&vOut[i])->add_coef(a.get_coef(i), 0, 0);
        }
    }
}

/**
 *  \brief Sub the algebraic vector a to the order zero of vIn and set the result in vOut. Used with T = Ofts< Ofsc >. Requires add_coef from vector class.
*         Careful: this routine leads to the copy of an entire Ofts structure. Very heavy! Use only in tests.
 **/
template <typename T , typename U> void sub_coef(vector<U> const& a, vector<T> const& vIn, vector<T>& vOut)
{
    if(a.size() != vOut.size() )
    {
        cout << "Error in add_coef (matrix.tpp): lengths do not match. Nothing is done." << endl;
    }
    else
    {
        for(int i =0; i < (int) vOut.size() ; i++)
        {
                ((T*)&vOut[i])->ccopy(vIn[i]);
                ((T*)&vOut[i])->add_coef(-a[i], 0, 0);
        }
    }
}


//----------------------------------------------------------------------------------------
//Functions used with T = Ofs<U>
//----------------------------------------------------------------------------------------
/**
 *  \brief Matrix-vector product: vOut += a x vIn. Used with T = Ofsc. Requires sprod.
 **/
inline void smvprod_ofs(matrix<Ofsc> const& a, vector<Ofsc> const& vIn, vector<Ofsc>& vOut)
{
    if((unsigned int)  a.get_size(2) != vIn.size() || (unsigned int) a.get_size(1) != vOut.size() )
    {
        cout << "Error in smvprod_ofs (matrix.tpp): lengths do not match. Nothing is done." << endl;
    }
    else
    {
        for(int i =0; i < a.get_size(1) ; i++)
        {
            for(int j =0; j < a.get_size(2); j++)
            {
                ( (Ofsc*)&vOut[i])->ofs_sprod(a.get_coef(i,j), vIn[j]);
            }
        }
    }
}

//----------------------------------------------------------------------------------------
//Functions used with U = cdouble
//----------------------------------------------------------------------------------------
/**
 *  \brief Matrix-vector product. Used with T = Ofs<U>, U = cdouble. Requires evaluate.
 *         Applied for the scalar version of the change of coordinates (see COC).
 **/
template <typename U> void mvprod_u(matrix< Ofs<U> > const& a, vector<U> const& vIn, vector<U>& vOut, double const& t)
{
    //Zeroing the target
    for(int i = 0; i< (int) vOut.size(); i++) vOut[i] = 0+0.0*I;
    //Loop on coefficients
    for(int i =0; i < a.get_size(1) ; i++)
    {
        for(int j =0; j < a.get_size(2); j++)
        {
            vOut[i] += vIn[j]*a.get_ptr_first_coef(i,j)->evaluate(t);
        }
    }
}

/**
 *  \brief Matrix-vector product (with sum). Used with T = Ofs<U>, U = cdouble. Requires evaluate.
 *         Applied for the scalar version of the change of coordinates (see COC).
 **/
template <typename U> void smvprod_u(matrix< Ofs<U> > const& a, vector<U> const& vIn, vector<U>& vOut, double const& t)
{
    //Loop on coefficients
    for(int i =0; i < a.get_size(1) ; i++)
    {
        for(int j =0; j < a.get_size(2); j++)
        {
            vOut[i] += vIn[j]*a.get_ptr_first_coef(i,j)->evaluate(t);
        }
    }
}

/**
 *  \brief Matrix-matrix product a x b. Used with T = Ofs<U>, U = cdouble. Requires ofs_sprod.
 *         Applied for the scalar version of the change of coordinates (see COC).
 **/
template <typename U> void smmprod_u(matrix< Ofs<U> > const& a, matrix< Ofs<U> > const& b, matrix< Ofs<U> >& mOut)
{
    //Loop on coefficients
    for(int i =0; i < mOut.get_size(1) ; i++)
    {
        for(int j =0; j < mOut.get_size(2); j++)
        {
            for(int k =0; k < a.get_size(2); k++) mOut.get_ptr_first_coef(i,j)->ofs_sprod(*a.get_ptr_first_coef(i,k), *b.get_ptr_first_coef(k,j));
        }
    }
}

/**
 *  \brief Matrix-matrix product a^T x b. Used with T = Ofs<U>, U = cdouble. Requires ofs_sprod.
 *         Applied for the scalar version of the change of coordinates (see COC).
 **/
template <typename U> void smtmprod_u(matrix< Ofs<U> > const& a, matrix< Ofs<U> > const& b, matrix< Ofs<U> >& mOut)
{
    //Loop on coefficients
    for(int i =0; i < mOut.get_size(1) ; i++)
    {
        for(int j =0; j < mOut.get_size(2); j++)
        {
            for(int k =0; k < a.get_size(2); k++) mOut.get_ptr_first_coef(i,j)->ofs_sprod(*a.get_ptr_first_coef(k,i), *b.get_ptr_first_coef(k,j));
        }
    }
}

/**
 *  \brief vOut += a(t). Used with T = Ofs<U>, U = cdouble. Requires ofs_sprod.
 *         Applied for the scalar version of the change of coordinates (see COC).
 **/
template <typename U> void vvsum_u(vector< Ofs<U> >& a, vector<U>& vOut, double const& t)
{
    for(int i =0; i < (int) vOut.size() ; i++)
    {
        vOut[i] += a[i].evaluate(t);
    }
}

/**
 *  \brief vOut -= a(t). Used with T = Ofs<U>, U = cdouble. Requires ofs_sprod.
 *         Applied for the scalar version of the change of coordinates (see COC).
 **/
template <typename U> void vvsub_u(vector< Ofs<U> >& a, vector<U>& vOut, double const& t)
{
    for(int i =0; i < vOut.size() ; i++) vOut[i] -= a[i].evaluate(t);
}

/**
 *  \brief vOut = vIn - a(t). Used with T = Ofs<U>, U = cdouble. Requires ofs_sprod.
 *         Applied for the scalar version of the change of coordinates (see COC).
 **/
template <typename U> void vvsub_u(vector< Ofs<U> >& a, vector<U> const& vIn, vector<U>& vOut, double const& t)
{
    for(int i =0; i < (int) vOut.size() ; i++) vOut[i] = vIn[i]-a[i].evaluate(t);
}

//----------------------------------------------------------------------------------------
// Read & Write
//----------------------------------------------------------------------------------------
/**
 * \brief Writes a given matrix W of type \c Ofts<Ofsc >  in a binary files of the form "filename+i+j.bin", with i = 0, m_size_r(W)-1
 *        and j = 0, m_size_c(W)-1.
 **/
inline void write_mofts_bin(matrix<Ofts<Ofsc > > &W, string filename)
{
    string ss1, ss2;
    //Loop on all coefficients
    for(int i = 0; i < W.get_size(1); i++)
    {
        for(int j = 0; j < W.get_size(2); j++)
        {
            ss1 = static_cast<ostringstream*>( &(ostringstream() << i) )->str();
            ss2 = static_cast<ostringstream*>( &(ostringstream() << j) )->str();
            write_ofts_bin(*W.get_ptr_first_coef(i,j), (filename+"["+ss1+"]"+"["+ss2+"].bin"));
        }
    }
}

/**
 * \brief Reads a given matrix W of type \c Ofts<Ofsc >  in a binary files of the form "filename+i+j.bin", with i = 0, m_size_r(W)-1
 *        and j = 0, m_size_c(W)-1.
 **/
inline void read_mofts_bin(matrix<Ofts<Ofsc > > &W, string filename)
{
    string ss1, ss2;
    //Loop on all coefficients
    for(int i = 0; i < W.get_size(1); i++)
    {
        for(int j = 0; j < W.get_size(2); j++)
        {
            ss1 = static_cast<ostringstream*>( &(ostringstream() << i) )->str();
            ss2 = static_cast<ostringstream*>( &(ostringstream() << j) )->str();
            read_ofts_bin(*W.get_ptr_first_coef(i,j), (filename+"["+ss1+"]["+ss2+"].bin"));
        }
    }
}


//########################################################################################
// End of implementation of the matrix template class
//########################################################################################
