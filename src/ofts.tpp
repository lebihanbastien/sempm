//########################################################################################
// Implementation of the Ofts template class
//########################################################################################

/**
 * \file ofts.tpp
 * \brief Fourier-Taylor series template class. See ofts.h for details.
 * \author BLB
 */

//----------------------------------------------------------------------------------------
//Create
//----------------------------------------------------------------------------------------
/**
 *  \brief Default constructor of the class Ofts<T>.
 */
template<typename T> Ofts<T>::Ofts()
{
    int i, index;
    m_ofts_nvar     = REDUCED_NV;
    m_ofts_order  = OFTS_ORDER;
    m_ofs_nvar    = Csts::OFS_NV;
    corder = OFS_ORDER;

    //New array of coefficients
    m_ofts_coefs = (Ofs<complex double>*) calloc(binomial(m_ofts_nvar+m_ofts_order,m_ofts_nvar), sizeof(Ofs<complex double>)); //new T[binomial(m_ofts_nvar+m_ofts_order, m_ofts_nvar)]();
    for(unsigned k= 0; k< binomial(m_ofts_nvar + m_ofts_order, m_ofts_nvar); k++) m_ofts_coefs[k] = Ofs<complex double>(corder);

    //Allocation of the homogeneous polynomials
    m_ofts_term = new Oftsh<T>*[m_ofts_order+1];
    //m_ofts_term = (Oftsh<T>**) calloc(binomial(m_ofts_nvar+m_ofts_order,m_ofts_nvar), sizeof(Oftsh<T>*));//new Oftsh<T>*[m_ofts_order+1];

    //Allocation of the coefficients
    index = 0;
    for(i=0; i<=m_ofts_order; i++)
    {
        //Allocation of each hp
        m_ofts_term[i] = new Oftsh<T>(m_ofts_nvar, i);
        //Link h to coefs at each level of the tree
        m_ofts_term[i]->link_coefs(m_ofts_coefs+index);
        index+=  Manip::nmon(m_ofts_nvar,i);
    }
}

/**
 *  \brief Constructor with given order/number of variables both for the Fourier-Taylor series and the coefficients.
 */
template<typename T> Ofts<T>::Ofts(int newNv, int newOrder, int newCnv, int newCorder)
{
    int i, index;

    m_ofts_nvar = newNv;
    m_ofts_order = newOrder;
    m_ofs_nvar = newCnv;
    corder = newCorder;

    //New array of coefficients
    m_ofts_coefs = (T*) calloc(binomial(m_ofts_nvar+m_ofts_order,m_ofts_nvar), sizeof(T)); //new T[binomial(m_ofts_nvar+m_ofts_order, m_ofts_nvar)]();

    //Allocation of the homogeneous polynomials
    m_ofts_term = new Oftsh<T>*[m_ofts_order+1];

    //Allocation of the coefficients
    index = 0;
    for(i=0; i<=m_ofts_order; i++)
    {
        //Allocation of each hp
        m_ofts_term[i] = new Oftsh<T>(m_ofts_nvar, i);//allocate_homog(m_ofts_nvar, i);
        //Link h to coefs at each level of the tree
        m_ofts_term[i]->link_coefs(m_ofts_coefs+index);
        //m_ofts_term[i]->link_coefs(new T[Manip::nmon(m_ofts_nvar,i)]);
        index+=  Manip::nmon(m_ofts_nvar,i);
    }
}

/**
 *  \brief Constructor with given order/number of variables both for the Fourier-Taylor series and the coefficients, Ofsc case.
 */
template<> inline Ofts< Ofs<complex double> >::Ofts(int newNv, int newOrder, int newCnv, int newCorder)
{
    int i, index;

    m_ofts_nvar = newNv;
    m_ofts_order = newOrder;
    m_ofs_nvar = newCnv;
    corder = newCorder;

    //New array of coefficients
    m_ofts_coefs = (Ofs<complex double>*) calloc(binomial(m_ofts_nvar+m_ofts_order,m_ofts_nvar), sizeof(Ofs<complex double>)); //new T[binomial(m_ofts_nvar+m_ofts_order, m_ofts_nvar)]();
    for(unsigned k= 0; k< binomial(m_ofts_nvar + m_ofts_order, m_ofts_nvar); k++) m_ofts_coefs[k]  = Ofs<complex double>(corder);

    //Allocation of the homogeneous polynomials
    m_ofts_term = new Oftsh< Ofs<complex double> >*[m_ofts_order+1];

    //Allocation of the coefficients
    index = 0;
    for(i=0; i<=m_ofts_order; i++)
    {
        //Allocation of each hp
        m_ofts_term[i] = new Oftsh< Ofs<complex double> >(m_ofts_nvar, i);//allocate_homog(m_ofts_nvar, i);
        //Link h to coefs at each level of the tree
        m_ofts_term[i]->link_coefs(m_ofts_coefs+index);
        //m_ofts_term[i]->link_coefs(new T[Manip::nmon(m_ofts_nvar,i)]);
        index+=  Manip::nmon(m_ofts_nvar,i);
    }
}

/**
 *  \brief Constructor with given order/number of variables both for the Fourier-Taylor series and the coefficients. Ofts< Ofsd > case.
 */
template<> inline Ofts< Ofsd >::Ofts(int newNv, int newOrder, int newCnv, int newCorder)
{
    int i, index;

    m_ofts_nvar = newNv;
    m_ofts_order = newOrder;
    m_ofs_nvar = newCnv;
    corder = newCorder;

    //New array of coefficients
    //-------------------------------------------------------
    //Replacing T *m_ofts_coefs = new T[binomial(m_ofts_nvar+m_ofts_order, m_ofts_nvar)]();
    m_ofts_coefs = (Ofsd*) calloc(binomial(m_ofts_nvar+m_ofts_order,m_ofts_nvar), sizeof(Ofsd));
    for(unsigned k= 0; k< binomial(m_ofts_nvar + m_ofts_order, m_ofts_nvar); k++) m_ofts_coefs[k]  = Ofsd(corder);
    //-------------------------------------------------------

    //Allocation of the homogeneous polynomials
    m_ofts_term = new Oftsh< Ofsd >*[m_ofts_order+1];

    //Allocation of the coefficients
    index = 0;
    for(i=0; i<=m_ofts_order; i++)
    {
        //Allocation of each hp
        m_ofts_term[i] = new Oftsh< Ofsd >(m_ofts_nvar, i);//allocate_homog(m_ofts_nvar, i);
        //Link h to coefs at each level of the tree
        m_ofts_term[i]->link_coefs(m_ofts_coefs+index);
        index+=  Manip::nmon(m_ofts_nvar,i);
    }
}

/**
 *  \brief Constructor from a given Ofts object (without any link).
 */
template<typename T> Ofts<T>::Ofts(Ofts<T> const& b)
{
    int nrc, index;

    //Same m_ofts_nvar/m_ofts_order
    m_ofts_nvar = b.m_ofts_nvar;
    m_ofts_order = b.m_ofts_order;
    m_ofs_nvar = b.m_ofs_nvar;
    corder = b.corder;

    //Copy of all the coefficients at every order in new array
    //m_ofts_coefs = new T[binomial(m_ofts_nvar+b.m_ofts_order, b.m_ofts_nvar)]();
    m_ofts_coefs = (Ofs<complex double>*) calloc(binomial(m_ofts_nvar+m_ofts_order,m_ofts_nvar), sizeof(Ofs<complex double>)); //new T[binomial(m_ofts_nvar+m_ofts_order, m_ofts_nvar)]();
    for(unsigned k= 0; k< binomial(m_ofts_nvar + m_ofts_order, m_ofts_nvar); k++) m_ofts_coefs[k]  = Ofs<complex double>(corder);

    index = 0;
    for(int nrc=0; nrc<= m_ofts_order; nrc++)
    {
        for (int i=0; i< Manip::nmon(m_ofts_nvar, nrc); i++)  m_ofts_coefs[index+i] = b.m_ofts_term[nrc]->get_coef(i);
        index+=  Manip::nmon(m_ofts_nvar,nrc);
    }

    //Allocation of the homogeneous polynomials
    m_ofts_term = new Oftsh<T>*[m_ofts_order+1];

    //Allocation of the coefficients
    index = 0;
    for(nrc=0; nrc<=m_ofts_order; nrc++)
    {
        //Allocation of each hp
        m_ofts_term[nrc] = new Oftsh<T>(m_ofts_nvar, nrc);//allocate_homog(m_ofts_nvar, i);
        //Link h to coefs at each level of the tree
        m_ofts_term[nrc]->link_coefs(m_ofts_coefs+index);
        index+=  Manip::nmon(m_ofts_nvar,nrc);
    }
}

/**
 *  \brief  An operator. Constructor from a given Ofts object (only the coefficients).
 */
template<typename T> Ofts<T>& Ofts<T>::operator = (Ofts<T> const& b)
{
    if(this != &b)
    {
        int nrc, index;

        //Same m_ofts_nvar/m_ofts_order
        m_ofts_nvar = b.m_ofts_nvar;
        m_ofts_order = b.m_ofts_order;
        m_ofs_nvar = b.m_ofs_nvar;
        corder = b.corder;

        //Copy of all the coefficients at every order in new array
        T *coef0 = new T[binomial(m_ofts_nvar+m_ofts_order, m_ofts_nvar)]();
        index = 0;
        for(int nrc=0; nrc<= m_ofts_order; nrc++)
        {
            for (int i=0; i< Manip::nmon(m_ofts_nvar, nrc); i++)  coef0[index+i] = b.m_ofts_term[nrc]->get_coef(i);
            index+=  Manip::nmon(m_ofts_nvar,nrc);
        }

        delete m_ofts_term;

        //Allocation of the homogeneous polynomials
        m_ofts_term = new Oftsh<T>*[m_ofts_order+1];

        //Allocation of the coefficients
        index = 0;
        for(nrc=0; nrc<=m_ofts_order; nrc++)
        {
            //Allocation of each hp
            m_ofts_term[nrc] = new Oftsh<T>(m_ofts_nvar, nrc);//allocate_homog(m_ofts_nvar, i);
            //Link h to coefs at each level of the tree
            m_ofts_term[nrc]->link_coefs(coef0+index);
            index+=  Manip::nmon(m_ofts_nvar,nrc);
        }

    }

    return *this;
}


//----------------------------------------------------------------------------------------
//Delete
//----------------------------------------------------------------------------------------
/**
 *  \brief Default destructor of the class Ofts<T>.
 *         WARNING: memory leak here, through the terms of type Oftsh.
 */
template<typename T> Ofts<T>::~Ofts<T>()
{
    //if(m_ofts_coefs != NULL) delete m_ofts_coefs;
    if(m_ofts_term != NULL)
    {
        //Certainly a problem at this point: since the delete routine of Oftsh is empty,
        // only the first leaf of the Oftsh tree is deleted...
        for(int i =0; i<= m_ofts_order ; i++) delete m_ofts_term[i];
    }
}

//----------------------------------------------------------------------------------------
//Copy
//----------------------------------------------------------------------------------------
/**
 *  \brief  Linked copy from a given Ofs object (exact same object is obtained).
 */
template<typename T> Ofts<T>& Ofts<T>::lcopy (Ofts<T> const& b)
{
    m_ofts_order = b.m_ofts_order;
    m_ofts_nvar = b.m_ofts_nvar;
    m_ofts_term = b.m_ofts_term;
    m_ofs_nvar = b.m_ofs_nvar;
    corder = b.corder;
    return *this;
}

/**
 *  \brief  Copy from a given Ofts object (only the coefficients).
 *
 *  Restricted to same order, same number of variables
 */
template<typename T> Ofts<T>& Ofts<T>::ccopy (Ofts<T> const& b)
{
    if(m_ofts_order != b.m_ofts_order || m_ofts_nvar != b.m_ofts_nvar || corder != b.corder || m_ofs_nvar != b.m_ofs_nvar)
    {
        cout << "Error using ccopy : the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        //Copy of all the coefficients at every order
        for(int nrc=0; nrc<= m_ofts_order; nrc++) for (int i=0; i< Manip::nmon(m_ofts_nvar, nrc); i++)  m_ofts_term[nrc]->set_coef(b.m_ofts_term[nrc]->get_coef(i), i);
        return *this;
    }
}

/**
 *  \brief  Copy from a given Ofts object (only the coefficients) at order nrc.
 *
 *  Restricted to same order, same number of variables
 */
template<typename T> Ofts<T>& Ofts<T>::ccopy (Ofts<T> const& b, int const& nrc)
{
    if(nrc > min(m_ofts_order, b.m_ofts_order) || m_ofts_nvar != b.m_ofts_nvar || corder != b.corder || m_ofs_nvar != b.m_ofs_nvar)
    {
        cout << "Error using ccopy : the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        //Copy of all the coefficients at every order
        for (int i=0; i< Manip::nmon(m_ofts_nvar, nrc); i++)  m_ofts_term[nrc]->set_coef(b.m_ofts_term[nrc]->get_coef(i), i);
        return *this;
    }
}


//----------------------------------------------------------------------------------------
//Setters
//----------------------------------------------------------------------------------------
/**
 *  \brief Sets a coefficient at a given order \c ord and a given position \c i at this order in the series.
 *
 * Set of given coefficient at term ord and position i in this term
 *   - ord gives the order of the homogeneous polynomial (hp)
 *   - i gives the positions wihtin this hp
 * Ex: set_coef(m, 2, 0) sets the x^2 coef
 * Ex: set_coef(m, 2, 1) sets the x*y coef
 */
template<typename T> void Ofts<T>::set_coef(T const& m, int ord, int i)
{
    m_ofts_term[ord]->set_coef(m, i);
}

/**
 *  \brief Adds a coefficient at a given order \c ord and a given position \c i at this order in the series.
 *
 * Set of given coefficient at term ord and position i in this term
 *   - ord gives the order of the homogeneous polynomial (hp)
 *   - i gives the positions wihtin this hp
 * Ex: add_coef(m, 2, 0) adds to the x^2 coef
 * Ex: add_coef(m, 2, 1) adds to the x*y coef
 */
template<typename T> void Ofts<T>::add_coef(T const& m, int ord, int i)
{
    m_ofts_term[ord]->add_coef(m, i);
}

/**
 *  \brief Sets of a given double/complex (typename U) subcoefficient everywhere (at each order and in each coefficient).
 */
template<typename T> template < typename U > void Ofts<T>::set_all_coefs(U const& m)
{
    for(int nrc=0; nrc<= m_ofts_order; nrc++) for (int i=0; i< Manip::nmon(m_ofts_nvar, nrc); i++)  m_ofts_term[nrc]->set_sub_coef(m, i);
}

/**
 *  \brief Sets random coefficients to all positions in the series.
 */
template<typename T> void Ofts<T>::set_random_coefs()
{
    for(int nrc=0; nrc<= m_ofts_order/2; nrc++)
    {
        m_ofts_term[nrc]->set_random_coefs();
    }
}

/**
 *  \brief Sets of a given U subcoefficient at order zero of the coefficient at position \c i of term of order \c n.
 */
template<typename T> template < typename U > void Ofts<T>::set_coef0(U const& m, int const& ord, int const& i)
{
    //m_ofts_term[ord]->setT0Coef(m, i);
    m_ofts_term[ord]->set_sub_coef(m, i, 0);
}


//----------------------------------------------------------------------------------------
//Getters
//----------------------------------------------------------------------------------------
/**
 *  \brief  Gets the order of the series.
 */
template<typename T> int Ofts<T>::get_order() const
{
    return m_ofts_order;
}

/**
 *  \brief  Gets the order of the coefficients.
 */
template<typename T> int Ofts<T>::get_coef_order() const
{
    return corder;
}

/**
 *  \brief  Gets the number of variables of series.
 */
template<typename T> int Ofts<T>::get_nvar() const
{
    return m_ofts_nvar;
}

/**
 *  \brief  Gets the number of variables of the coefficients.
 */
template<typename T> int Ofts<T>::get_coef_nvar() const
{
    return m_ofs_nvar;
}

/**
 *  \brief  Gets the adress of the coefficient at order \c ord and position \c pos
 */
template<typename T>  T* Ofts<T>::get_coef(int const& ord, int const& pos) const
{
    if(ord > m_ofts_order || pos >= Manip::nmon(m_ofts_nvar, ord))
    {
        cout << "Error in get_coef: out of range. First m_ofts_term is returned" << endl;
        cout << "Requested m_ofts_order: " << ord << ", Maximum allowed: " <<  m_ofts_order << endl;
        cout << "Requested pos: " << pos << ", Maximum allowed: " <<  Manip::nmon(m_ofts_nvar, ord) << endl;
        return this->m_ofts_term[0]->get_ptr_first_coef();
    }
    else return this->m_ofts_term[ord]->get_ptr_first_coef()+pos;
}

/**
 *  \brief  Gets the adress of the term at order \c ord
 */
template<typename T> Oftsh<T>* Ofts<T>::get_term(int const& ord) const
{
    if(ord > m_ofts_order)
    {
        cout << "Error in get_term: out of range. First m_ofts_term is returned" << endl;
        cout << "Requested m_ofts_order: " << ord << ", Maximum allowed: " <<  m_ofts_order << endl;
        return this->m_ofts_term[0];
    }
    else return this->m_ofts_term[ord];
}

/**
 *  \brief  Gets the adress of the Ofts object
 */
template <typename T> Ofts<T>* Ofts<T>::get_ptr() const
{
    return (Ofts<T>*) this;
}



//----------------------------------------------------------------------------------------
//Zeroing
//----------------------------------------------------------------------------------------
/**
 *  \brief  Sets all coefficients to zero.
 */
template<typename T> void Ofts<T>::zero()
{
    for(int nrc=0; nrc<= m_ofts_order; nrc++) m_ofts_term[nrc]->zero();
}


//--------------------------------------------------------------------------------
//Operations
//--------------------------------------------------------------------------------
//------------------
// Conjugate
//------------------
/**
 *  \brief Conjugates  all terms (Oftsh object), and only them! To be used with evaluate_conjugate to have the true conjugate.
 */
template<typename T> Ofts<T>& Ofts<T>::conjugate()
{
    for(int nrc=0; nrc<= m_ofts_order; nrc++) m_ofts_term[nrc]->conjugate();
    return *this;
}

/**
 *  \brief Conjugates the order \c nrc. To be used with evaluate_conjugate to have the true conjugate.
 */
template<typename T> Ofts<T>& Ofts<T>::conjugate(int const& nrc)
{
    m_ofts_term[nrc]->conjugate();
    return *this;
}

//------------------
// smult
//------------------
/**
 *  \brief  An operation. Adds the product: \c this \f$  += m a \f$ with m a coefficient.
 */
template<typename T> Ofts<T>& Ofts<T>::ofts_smult_t(Ofts<T> const& a, T const& m)
{
    if(m_ofts_order != a.m_ofts_order || m_ofts_nvar != a.m_ofts_nvar)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        for(int k=0; k <= m_ofts_order; k++) m_ofts_term[k]->oftsh_smult_t(*a.m_ofts_term[k], m);  //ps[k] += m*a[k]
        return *this;
    }
}

/**
 *  \brief  An operation. Adds the product: \c this \f$  += m a \f$ with m a coefficient, at order n
 *  \param  a: a reference to an Oftsh object
 *  \param  m: a reference to a coefficient
 *  \param  n: a reference to the order to update
 */
template<typename T> Ofts<T>& Ofts<T>::ofts_smult_t(Ofts<T> const& a, T const& m, int const& n)
{
    if(m_ofts_order != a.m_ofts_order || m_ofts_nvar != a.m_ofts_nvar)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        m_ofts_term[n]->oftsh_smult_t(*a.m_ofts_term[n], m);  //ps[n] += m*a[n]
        return *this;
    }
}


/**
 *  \brief  An operation. Sets the product: \c this \f$  = m a \f$ with m a coefficient, at order n
 *  \param  a: a reference to an Oftsh object
 *  \param  m: a reference to a coefficient
 *  \param  n: a reference to the order to update
 */
template<typename T> Ofts<T>& Ofts<T>::ofts_mult_t(Ofts<T> const& a, T const& m, int const& n)
{
    if(m_ofts_order != a.m_ofts_order || m_ofts_nvar != a.m_ofts_nvar)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        m_ofts_term[n]->oftsh_mult_t(*a.m_ofts_term[n], m);  //ps[n] += m*a[n]
        return *this;
    }
}

//---------------------------------------

/**
 *  \brief  An operation. Adds the product: \c this \f$  += c a \f$ with c a subcoefficient.
 *  \param  a: a reference to an Oftsh object
 *  \param  c: a reference to a subcoefficient
 */
template<typename T> template<typename U> Ofts< Ofs<U> >& Ofts<T>::ofts_smult_u(Ofts< Ofs<U> > const& a, U const& c)
{
    if(m_ofts_order != a.m_ofts_order || m_ofts_nvar != a.m_ofts_nvar)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        for(int k=0; k <= m_ofts_order; k++) m_ofts_term[k]->oftsh_smult_u(*a.m_ofts_term[k], c);  //ps[k] += m*a[k]
        return *this;
    }
}

/**
 *  \brief  An operation. Adds the product: \c this \f$  += c a \f$ with c a subcoefficient at order k.
 *  \param  a: a reference to an Oftsh object
 *  \param  c: a reference to a subcoefficient
 */
template<typename T> template<typename U> Ofts< Ofs<U> >& Ofts<T>::ofts_smult_u(Ofts< Ofs<U> > const& a, U const& c, int const& k)
{
    if(k > m_ofts_order || m_ofts_nvar != a.m_ofts_nvar)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        m_ofts_term[k]->oftsh_smult_u(*a.m_ofts_term[k], c);  //ps[k] += m*a[k]
        return *this;
    }
}

//---------------------------------------

/**
 *  \brief  An operation. Adds the product: \c this \f$  += c m a \f$ with m a coefficient and c a subcoefficient.
 *  \param  a: a reference to an Oftsh object
 *  \param  c: a reference to a subcoefficient
 */
template<typename T> template<typename U> Ofts< Ofs<U> >& Ofts<T>::ofts_smult_tu(Ofts< Ofs<U> > const& a, T const& m, U const& c)
{
    if(m_ofts_order != a.m_ofts_order || m_ofts_nvar != a.m_ofts_nvar)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {

        int k;
        for(k=0; k <= m_ofts_order; k++)
        {
            m_ofts_term[k]->oftsh_smult_tu(*a.m_ofts_term[k], m, c);  //ps[k] += m*a[k]
        }

        return *this;
    }
}

/**
 *  \brief  An operation. Adds the product: \c this \f$  += c m a \f$ with m a coefficient and c a subcoefficient, at m_ofts_order k.
 *  \param  a: a reference to an Oftsh object
 *  \param  c: a reference to a subcoefficient
 */
template<typename T> template<typename U> Ofts< Ofs<U> >& Ofts<T>::ofts_smult_tu(Ofts< Ofs<U> > const& a, T const& m, U const& c, int const& k)
{
    if(k > m_ofts_order  || m_ofts_nvar != a.m_ofts_nvar)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        m_ofts_term[k]->oftsh_smult_tu(*a.m_ofts_term[k], m, c);  //ps[k] += m*a[k]
        return *this;
    }
}


//------------------
// mult
//------------------
/**
 *  \brief  An operation. Sets the product: \c this \f$  += m a \f$ with m a coefficient.
 *  \param  a: a reference to an Oftsh object
 *  \param  m: a reference to a coefficient
 */
template<typename T> Ofts<T>& Ofts<T>::ofts_mult_t(Ofts<T> const& a, T const& m)
{
    if(m_ofts_order != a.m_ofts_order || m_ofts_nvar != a.m_ofts_nvar)
    {
        cout << "Error using mult: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        for(int k=0; k <= m_ofts_order; k++)
        {
            m_ofts_term[k]->mult(*a.m_ofts_term[k], m);  //ps[k] = m*a[k], contains the zeroing
        }

        return *this;
    }
}

/**
 *  \brief  An operation. Sets the product: \c this \f$  = c a \f$ with c a subcoefficient at order k.
 *  \param  a: a reference to an Oftsh object
 *  \param  c: a reference to a subcoefficient
 */
template<typename T> template<typename U> Ofts< Ofs<U> >& Ofts<T>::ofts_mult_u(Ofts< Ofs<U> > const& a, U const& c, int const& k)
{
    if(k > m_ofts_order || m_ofts_nvar != a.m_ofts_nvar)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        m_ofts_term[k]->oftsh_mult_u(*a.m_ofts_term[k], c);  //ps[k] += m*a[k]
        return *this;
    }
}

/**
 *  \brief  An operation. Sets the product: \c this \f$  = c a \f$ with c a subcoefficient.
 *  \param  a: a reference to an Oftsh object
 *  \param  c: a reference to a subcoefficient
 */
template<typename T> template<typename U> Ofts< Ofs<U> >& Ofts<T>::ofts_mult_u(Ofts< Ofs<U> > const& a, U const& c)
{
    if(m_ofts_order != a.m_ofts_order || m_ofts_nvar != a.m_ofts_nvar)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {

        int k;
        for(k=0; k <= m_ofts_order; k++)
        {
            m_ofts_term[k]->oftsh_mult_u(*a.m_ofts_term[k], c);  //ps[k] += m*a[k]
        }

        return *this;
    }
}

/**
 *  \brief  An operation. Sets the product: \c this \f$  = c m a \f$ with m a coefficient and c a subcoefficient.
 *  \param  a: a reference to an Oftsh object
 *  \param  c: a reference to a subcoefficient
 */
template<typename T> template<typename U> Ofts< Ofs<U> >& Ofts<T>::ofts_mult_tu(Ofts< Ofs<U> > const& a, T const& m, U const& c)
{
    if(m_ofts_order != a.m_ofts_order || m_ofts_nvar != a.m_ofts_nvar)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {

        int k;
        for(k=0; k <= m_ofts_order; k++)
        {
            m_ofts_term[k]->mult(*a.m_ofts_term[k], m, c);  //ps[k] += m*a[k]
        }

        return *this;
    }
}


//------------------
// sfsum
//------------------
/**
 *  \brief  An operation. Adds the product: \c this \f$  = ma*a + mb*b \f$ with ma and mb coefficients.
 *  \param  a: a reference to an Ofts object
 *  \param  ma: a reference to a coefficient
 *  \param  b: a reference to an Ofts object
 *  \param  mn: a reference to a coefficient
 */
template<typename T> Ofts<T>& Ofts<T>::ofts_sfsum_t(Ofts<T> const& a, T const& ma, Ofts<T> const& b, T const& mb)
{
    if(m_ofts_order != a.m_ofts_order || m_ofts_nvar != a.m_ofts_nvar || m_ofts_order != b.m_ofts_order || m_ofts_nvar != b.m_ofts_nvar)
    {
        cout << "Error using sfsum: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        for(int k=0; k <= m_ofts_order; k++)
        {
            m_ofts_term[k]->oftsh_smult_t(*a.m_ofts_term[k], ma);   //ps[k]  += ma*a[k]
            m_ofts_term[k]->oftsh_smult_t(*b.m_ofts_term[k], mb);   //ps[k]  += mb*b[k]
        }
        return *this;
    }

}

/**
 *  \brief  An operation. Adds the product: \c this \f$  += ma*a + mb*b \f$ with ma and mb coefficients at order n.
 *  \param  a: a reference to an Ofts object
 *  \param  ma: a reference to a coefficient
 *  \param  b: a reference to an Ofts object
 *  \param  mn: a reference to a coefficient
 */
template<typename T> Ofts<T>& Ofts<T>::ofts_sfsum_t(Ofts<T> const& a, T const& ma, Ofts<T> const& b, T const& mb, int n)
{
    if(m_ofts_order != a.m_ofts_order || m_ofts_nvar != a.m_ofts_nvar || m_ofts_order != b.m_ofts_order || m_ofts_nvar != b.m_ofts_nvar)
    {
        cout << "Error using sfsum: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        m_ofts_term[n]->oftsh_smult_t(*a.m_ofts_term[n], ma);   //ps[k]  += ma*a[k]
        m_ofts_term[n]->oftsh_smult_t(*b.m_ofts_term[n], mb);   //ps[k]  += mb*b[k]
        return *this;
    }

}

/**
 *  \brief  An operation. Adds the product: \c this \f$  += ma*a + mb*b + mc*c \f$ with ma and mb coefficients at order n.
 *  \param  a: a reference to an Ofts object
 *  \param  ma: a reference to a coefficient
 *  \param  b: a reference to an Ofts object
 *  \param  mb: a reference to a coefficient
 *  \param  c: a reference to an Ofts object
 *  \param  mc: a reference to a coefficient
 *  \param  n: a reference to the order to update
 */
template<typename T> Ofts<T>& Ofts<T>::ofts_sfsum_tt(Ofts<T> const& a, T const& ma, Ofts<T> const& b, T const& mb, Ofts<T> const& c, T const& mc, int n)
{
    if(m_ofts_order != a.m_ofts_order || m_ofts_nvar != a.m_ofts_nvar || m_ofts_order != b.m_ofts_order || m_ofts_nvar != b.m_ofts_nvar|| m_ofts_order != c.m_ofts_order || m_ofts_nvar != c.m_ofts_nvar)
    {
        cout << "Error using sfsum: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        m_ofts_term[n]->oftsh_smult_t(*a.m_ofts_term[n], ma);   //ps[k]  += ma*a[k]
        m_ofts_term[n]->oftsh_smult_t(*b.m_ofts_term[n], mb);   //ps[k]  += mb*b[k]
        m_ofts_term[n]->oftsh_smult_t(*c.m_ofts_term[n], mc);   //ps[k]  += mb*b[k]
        return *this;
    }

}

//------------------
// fsum
//------------------

/**
 *  \brief  An operation. Sets the product: \c this \f$  = ma*a + mb*b \f$ with ma and mb coefficients.
 *  \param  a: a reference to an Ofts object
 *  \param  ma: a reference to a coefficient
 *  \param  b: a reference to an Ofts object
 *  \param  mb: a reference to a coefficient
 */
template<typename T> Ofts<T>& Ofts<T>::ofts_fsum_t(Ofts<T> const& a, T const& ma, Ofts<T> const& b, T const& mb)
{
    if(m_ofts_order != a.m_ofts_order || m_ofts_nvar != a.m_ofts_nvar || m_ofts_order != b.m_ofts_order || m_ofts_nvar != b.m_ofts_nvar)
    {
        cout << "Error using fsum: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        int k;
        for(k=0; k <= m_ofts_order; k++)
        {
            m_ofts_term[k]->oftsh_mult_t(*a.m_ofts_term[k], ma);   //ps[k]   = ma*a[k], contains the zeroing
            m_ofts_term[k]->oftsh_smult_t(*b.m_ofts_term[k], mb);   //ps[k]  += mb*b[k]
        }
        return *this;
    }

}

/**
 *  \brief  An operation. Sets the product: \c this \f$  = ma*a + mb*b \f$ with ma and mb coefficients at order n.
 *  \param  a: a reference to an Ofts object
 *  \param  ma: a reference to a coefficient
 *  \param  b: a reference to an Ofts object
 *  \param  mb: a reference to a coefficient
 *  \param  n: a reference to the order to update
 */
template<typename T> Ofts<T>& Ofts<T>::ofts_fsum_t(Ofts<T> const& a, T const& ma, Ofts<T> const& b, T const& mb, int n)
{
    if(m_ofts_order != a.m_ofts_order || m_ofts_nvar != a.m_ofts_nvar || m_ofts_order != b.m_ofts_order || m_ofts_nvar != b.m_ofts_nvar)
    {
        cout << "Error using fsum: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        m_ofts_term[n]->oftsh_mult_t(*a.m_ofts_term[n], ma);    //ps[k]   = ma*a[k], contains the zeroing
        m_ofts_term[n]->oftsh_smult_t(*b.m_ofts_term[n], mb);   //ps[k]  += mb*b[k]
        return *this;
    }

}

/**
 *  \brief  An operation. Sets the product: \c this \f$  = ca*a + cb*b \f$ with ca and cb subcoefficients.
 *  \param  a: a reference to an Ofts object
 *  \param  ca: a reference to a subcoefficient
 *  \param  b: a reference to an Ofts object
 *  \param  cb: a reference to a subcoefficient
 */
template<typename T> template<typename U> Ofts<T>& Ofts<T>::ofts_fsum_u(Ofts<T> const& a, U const& ca, Ofts<T> const& b, U const& cb)
{
    if(m_ofts_order != a.m_ofts_order || m_ofts_nvar != a.m_ofts_nvar || m_ofts_order != b.m_ofts_order || m_ofts_nvar != b.m_ofts_nvar)
    {
        cout << "Error using fsum: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        for(int k=0; k <= m_ofts_order; k++)
        {
            m_ofts_term[k]->oftsh_mult_u(*a.m_ofts_term[k], ca);        //ps[k]   = ca*a[k], contains the zeroing
            m_ofts_term[k]->oftsh_smult_u(*b.m_ofts_term[k], cb);       //ps[k]  += mb*b[k]
        }
        return *this;
    }

}

/**
 *  \brief  An operation. Sets the product: \c this \f$  = ca*a + cb*b \f$ with ca and cb subcoefficients at order m
 *  \param  a: a reference to an Ofts object
 *  \param  ca: a reference to a subcoefficient
 *  \param  b: a reference to an Ofts object
 *  \param  cb: a reference to a subcoefficient
 *  \param  m: a reference to the order to update
 */
template<typename T> template<typename U> Ofts<T>& Ofts<T>::ofts_fsum_u(Ofts<T> const& a, U const& ca, Ofts<T> const& b, U const& cb, int const& m)
{
    if(m_ofts_order != a.m_ofts_order || m_ofts_nvar != a.m_ofts_nvar || m_ofts_order != b.m_ofts_order || m_ofts_nvar != b.m_ofts_nvar)
    {
        cout << "Error using fsum: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        m_ofts_term[m]->oftsh_mult_u(*a.m_ofts_term[m], ca);        //ps[k]   = ca*a[k], contains the zeroing
        m_ofts_term[m]->oftsh_smult_u(*b.m_ofts_term[m], cb);       //ps[k]  += mb*b[k]
        return *this;
    }

}


//------------------
// sprod
//------------------
/**
 *  \brief  An operation. Adds the product: \c this \f$  += a*b \f$.
 *
 *  Handle the case for which n >= max(a.m_ofts_order, b.m_ofts_order)
 */
template<typename T> Ofts<T>& Ofts<T>::ofts_sprod(Ofts<T> const& a, Ofts<T> const& b)
{
    int k, i, i0, i1;
    //Product
    for(k=0; k<=m_ofts_order; k++)
    {
        i0 = min(b.m_ofts_order, k);
        i1 = min(a.m_ofts_order, k);
        for(i= k-i0; i<=i1; i++) m_ofts_term[k]->oftsh_sprod(*a.m_ofts_term[i], *b.m_ofts_term[k-i]);
    }
    return *this;
}

/**
 *  \brief  An operation. Adds the product: \c this \f$  += a*b \f$ at order n.
 *
 *  Handle the case for which n >= max(a.m_ofts_order, b.m_ofts_order)
 */
template<typename T> Ofts<T>& Ofts<T>::ofts_sprod(Ofts<T> const& a, Ofts<T> const& b, int const& n)
{
    int i, i0, i1;
    i0 = min(b.m_ofts_order, n);
    i1 = min(a.m_ofts_order, n);
    //Product
    for(i= n-i0; i<=i1; i++) m_ofts_term[n]->oftsh_sprod(*a.m_ofts_term[i], *b.m_ofts_term[n-i]);
    return *this;
}

//------------------
// smprod
//------------------
/**
 *  \brief  An operation. Adds the product: \c this \f$  += m*a*b \f$ with m a coefficient.
 *  \param  a: a reference to an Ofts object
 *  \param  b: a reference to an Ofts object
 *  \param  n: a reference to a coefficient
 */
template<typename T> Ofts<T>& Ofts<T>::ofts_smprod_t(Ofts<T> const& a, Ofts<T> const& b, T const& m, T& temp)
{
    //------------------------------------------------------------------------------------
    // This routine does not use a temporary variable in m_ofts_order to avoid
    // additional Ofts object
    //------------------------------------------------------------------------------------
    int k, i, i0, i1;
    //Product
    for(k=0; k<=m_ofts_order; k++)
    {
        i0 = min(b.m_ofts_order, k);
        i1 = min(a.m_ofts_order, k);
        for(i= k-i0; i<=i1; i++) m_ofts_term[k]->oftsh_smprod_t(*a.m_ofts_term[i], *b.m_ofts_term[k-i], m, temp);
    }
    return *this;
}

/**
 *  \brief  An operation. Adds the product: \c this \f$  += m*a*b \f$ with m a coefficient at order n
 *  \param  a: a reference to an Ofts object
 *  \param  b: a reference to an Ofts object
 *  \param  n: a reference to a coefficient
 */
template<typename T> Ofts<T>& Ofts<T>::ofts_smprod_t(Ofts<T> const& a, Ofts<T> const& b, T const& m, int const& n, T& temp)
{
    //------------------------------------------------------------------------------------
    // This routine does not use a temporary variable in m_ofts_order to avoid
    // additional Ofts object
    //------------------------------------------------------------------------------------
    int i, i0, i1;

    i0 = min(b.m_ofts_order, n);
    i1 = min(a.m_ofts_order, n);
    for(i= n-i0; i<=i1; i++) m_ofts_term[n]->oftsh_smprod_t(*a.m_ofts_term[i], *b.m_ofts_term[n-i], m, temp);

    return *this;
}

/**
 *  \brief  An operation. Adds the product: \c this \f$  += c*a*b \f$ with c a subcoefficient.
 *  \param  a: a reference to an Ofts object
 *  \param  b: a reference to an Ofts object
 *  \param  c: a reference to a subcoefficient
 *
 *  NOTE: ofts_smprod_u currently uses a temporary Ofts<T> object, which is not recommended.
 *  This should be changed to match ofts_smprod_t standards (no temporary object to Ofts level).
 *  However, the routine is currently only used in tests routines (ofts_test.cpp), hence this
 *  matter is not critical.
 */
template<typename T> template<typename U> Ofts<T>& Ofts<T>::ofts_smprod_u(Ofts<T> const& a, Ofts<T> const& b, U const& c)
{
    //temp = 0 sizeOf(this)
    Ofts<T> temp(this->m_ofts_nvar, this->m_ofts_order, this->m_ofs_nvar, this->corder);
    //temp = a*b
    temp.ofts_sprod(a,b);
    //this = c*temp
    this->ofts_smult_u(temp,c);
    return *this;
}

/**
 *  \brief  An operation. Adds the product: \c this \f$  += m*a*b \f$ with m a coefficient at order n.
 *  \param  a: a reference to an Ofts object
 *  \param  b: a reference to an Ofts object
 *  \param  m: a reference to a coefficient
 *  \param  n: a reference to the order to update
 *
 *  NOTE: ofts_smprod_u currently uses a temporary Ofts<T> object, which is not recommended.
 *  This should be changed to match ofts_smprod_t standards (no temporary object to Ofts level).
 *  However, the routine is currently only used in tests routines (ofts_test.cpp), hence this
 *  matter is not critical.
 */
template<typename T> template<typename U> Ofts<T>& Ofts<T>::ofts_smprod_u(Ofts<T> const& a, Ofts<T> const& b, U const& c, int const& n)
{
    //temp = 0 sizeOf(this)
    Ofts<T> temp(this->m_ofts_nvar, this->m_ofts_order, this->m_ofs_nvar, this->corder);
    //temp = a*b at order n
    temp.ofts_sprod(a,b,n);
    //this = c*temp at order n
    this->ofts_smult_u(temp,c, n);
    return *this;
}


//------------------
// prod
//------------------
/**
 *  \brief  An operation. Sets the product: \c this \f$  = m*a*b \f$ with m a coefficient.
 *  \param  a: a reference to an Ofts object
 *  \param  b: a reference to an Ofts object
 *  \param  m: a reference to a coefficient
 */
template<typename T> Ofts<T>& Ofts<T>::mprod(Ofts<T> const& a, Ofts<T> const& b, T const& m)
{
    this->zero();
    this->smprod(a,b,m);
    return *this;
}

/**
 *  \brief  An operation. Sets the product: \c this \f$  = m*a*b \f$ with c a subcoefficient.
 *  \param  a: a reference to an Ofts object
 *  \param  b: a reference to an Ofts object
 *  \param  c: a reference to a subcoefficient
 */
template<typename T> template<typename U> Ofts<T>& Ofts<T>::mprod(Ofts<T> const& a, Ofts<T> const& b, U const& c)
{
    this->zero();
    this->smprod(a,b,c);
    return *this;
}

/**
 *  \brief  An operation. Sets the product: \c this \f$  = a*b \f$.
 *  \param  a: a reference to an Ofts object
 *  \param  b: a reference to an Ofts object
 */
template<typename T> Ofts<T>& Ofts<T>::prod(Ofts<T> const& a, Ofts<T> const& b)
{
    this->zero();
    this->sprod(a,b);
    return *this;
}


//------------------
// pows
//------------------
/**
 *   \brief Power function: p = a^alpha at order n. This generic routine is not directly used, only specialized versions are (in particular for the Ofts<Ofs> form).
 **/
template<typename T> template<typename U> Ofts<T>& Ofts<T>::ofts_pows(Ofts<T> const& a,  U const& alpha, int const& n)
{
    T x0;
    x0 = coef0s(a);
    if(n==0)
    {
        this->acoef0s(cpow(x0, alpha));
    }
    else
    {
        //Sets every coefficients to zero for order n
        this->m_ofts_term[n]->zero();
        for(int j=0; j<= n-1; j++) this->m_ofts_term[n]->smprod(*a.m_ofts_term[n-j], *this->m_ofts_term[j], alpha*(n-j)-j);// smprodh(ps->m_ofts_term[k], s->m_ofts_term[k-j], ps->m_ofts_term[j], alpha*(k-j)-j);
        this->m_ofts_term[n]->mult(1.0/(n*x0)); //multh(ps->m_ofts_term[k], 1.0/(x0*k));
    }

    return *this;
}

/**
 *   \brief Power function: p = a^alpha at order n. Ofsc case: ONLY FOR OFS_ORDER = 0
 **/
template<> template<typename U> Ofts< Ofsc >& Ofts< Ofsc >::ofts_pows(Ofts< Ofsc > const& a,  U const& alpha, int const& n)
{
    if(n==0)
    {
        Ofsc temp(corder);
        Ofsc temp2(corder);
        //Initialization of order zero
        temp.ofs_pows(*coef0s(a), alpha);
        this->acoef0s(temp);
    }
    else
    {
        Ofsc temp(corder);
        Ofsc temp2(corder);
        Ofsc temp3(corder);
        //temp2 = 1/x0
        temp2.ofs_pows(*coef0s(a), -1.0+0.0*I);
        //Sets every coefficients to zero for order n
        this->m_ofts_term[n]->zero();
        //Recurrence scheme @order n
        int i0 = min(a.m_ofts_order, n);
        for(int j= n-i0; j< n; j++)
        {
            temp.ofs_mult(temp2, (alpha*(n-j)-j)/n);
            this->m_ofts_term[n]->oftsh_smprod_t(*a.m_ofts_term[n-j], *this->m_ofts_term[j], temp, temp3);
        }
    }

    return *this;
}

//------------------
// Order 0 routines
//------------------
/**
 * \brief Returns the address of the first coefficient of order 0 of the taylor series s
 */
template<typename T> T* Ofts<T>::coef0s(Ofts<T> const& a)
{
    return a.m_ofts_term[0][0].get_ptr_first_coef();
}

/**
 * \brief Sets the coefficient of order 0 of the taylor series s equal to x0
 */
template<typename T> void Ofts<T>::acoef0s(T const& x0)
{
    this->m_ofts_term[0][0].set_coef(x0,0);
}

//----------------------------------------------------------------------------------------
//Operations with TFS coefficients - pure operations
//----------------------------------------------------------------------------------------
//----------------
// Frequency to Time domain
//---------------
/**
 *  \brief Puts all Ofs coefficients in a into Tfs coefficients and puts it into this.
 **/
template<typename T> Ofts<T>& Ofts<T>::tfs_from_ofs(Ofts<T> const& a)
{
    Ofs< cdouble > tfs(corder);
    for(int nrc=0; nrc<= m_ofts_order; nrc++)
    {
        //Current homogeneous polynomial
        for (int i=0; i< Manip::nmon(m_ofts_nvar, nrc); i++)
        {
            tfs.tfs_from_ofs(a.m_ofts_term[nrc]->get_coef(i));
            this->m_ofts_term[nrc]->set_coef(tfs, i);
        }
    }
    return *this;
}

/**
 *  \brief Puts all Ofs coefficients in a into Tfs coefficients and puts it into this. Inline.
 **/
template<typename T> Ofts<T>& Ofts<T>::tfs_from_ofs_inline()
{
    Ofs< cdouble > tfs(corder);
    for(int nrc=0; nrc<= m_ofts_order; nrc++)
    {
        //Current homogeneous polynomial
        for (int i=0; i< Manip::nmon(m_ofts_nvar, nrc); i++)
        {
            (this->m_ofts_term[nrc]->get_ptr_first_coef()+i)->tfs_from_ofs_inline(tfs);
        }
    }
    return *this;
}

/**
 *  \brief Puts all Ofs coefficients in a into Tfs coefficients and puts it into this. Inline.
 **/
template<typename T> Ofts<T>& Ofts<T>::tfs_from_ofs_inline(int nrc)
{
    Ofs< cdouble > tfs(corder);
    //Current homogeneous polynomial
    for (int i=0; i< Manip::nmon(m_ofts_nvar, nrc); i++)
    {
        (this->m_ofts_term[nrc]->get_ptr_first_coef()+i)->tfs_from_ofs_inline(tfs);
    }
    return *this;
}

/**
 *  \brief Puts all Tfs coefficients in a into ofs coefficients and puts it into this.
 *         Makes use of FFT routines from GSL on each Fourier coefficient.
 **/
template<typename T> Ofts<T>& Ofts<T>::tfs_to_ofs(Ofts<T> const& a)
{
    Ofs< cdouble > ofs(corder);
    for(int nrc=0; nrc<= m_ofts_order; nrc++)
    {
        //Current homogeneous polynomial
        for (int i=0; i< Manip::nmon(m_ofts_nvar, nrc); i++)
        {
            ofs.tfs_to_ofs(a.m_ofts_term[nrc]->get_coef(i));
            this->m_ofts_term[nrc]->set_coef(ofs, i);
        }
    }
    return *this;
}

/**
 *  \brief Puts all Tfs coefficients in a into ofs coefficients and puts it into this. Inline
 *         Makes use of FFT routines from GSL on each Fourier coefficient.
 **/
template<typename T> Ofts<T>& Ofts<T>::tfs_to_ofs_inline()
{
    for(int nrc=0; nrc<= m_ofts_order; nrc++)
    {
        //Current homogeneous polynomial
        for (int i=0; i< Manip::nmon(m_ofts_nvar, nrc); i++)
        {
            (this->m_ofts_term[nrc]->get_ptr_first_coef()+i)->tfs_to_ofs_inline();
        }
    }
    return *this;
}

/**
 *  \brief Puts all Tfs coefficients in a into ofs coefficients and puts it into this. Inline
 *         Makes use of FFT routines from GSL on each Fourier coefficient.
 **/
template<typename T> Ofts<T>& Ofts<T>::tfs_to_ofs_inline(int nrc)
{
    //Current homogeneous polynomial
    for (int i=0; i< Manip::nmon(m_ofts_nvar, nrc); i++)
    {
        (this->m_ofts_term[nrc]->get_ptr_first_coef()+i)->tfs_to_ofs_inline();
    }
    return *this;
}


//----------------
// Pows
//----------------
/**
 *   \brief Power function: p = a^alpha, with Tfs coefficients.
 *          WARNING: for alpha < 0, a good precision is achieved ONLY if
 *          the (Fourier) order 0 of the (Taylor) order 0 is >> wrt to the rest of the (Fourier) coefficients of the (Taylor) order 0.
 **/
template<typename T> template<typename U> Ofts<T>& Ofts<T>::tfts_pows(Ofts<T> const& a,  U const& alpha)
{
    cout << "tfts_pows is an empty shell if T != Ofsc" << endl;
    return *this;
}

/**
 *   \brief Power function: p = a^alpha, with Tfs coefficients. Double complex version.
 *          WARNING: for alpha < 0, a good precision is achieved ONLY if
 *          the (Fourier) order 0 of the (Taylor) order 0 is >> wrt to the rest of the (Fourier) coefficients of the (Taylor) order 0.
 **/
template<> template<typename U> Ofts< Ofsc >& Ofts< Ofsc >::tfts_pows(Ofts< Ofsc > const& a,   U const& alpha)
{
    Ofsc temp(corder);

    //Initialization of order zero
    coef0s(*this)->tfs_pows(*coef0s(a), alpha);

    //temp2 = 1/x0
    temp.tfs_pows(*coef0s(a), -1.0+0.0*I);

    //Recurrence scheme
    for(int k=1; k <= m_ofts_order; k++)
    {
        //Sets every coefficients to zero for m_ofts_order k
        this->m_ofts_term[k]->zero();
        //Loop on all previously computed homogeneous terms
        for(int j=0; j<= k-1; j++)
        {
            this->m_ofts_term[k]->tftsh_smprod_tu(*a.m_ofts_term[k-j], *this->m_ofts_term[j], temp, (alpha*(k-j)-j)/k);
        }
    }
    return *this;
}

/**
 *   \brief Power function: p = a^alpha, with Tfs coefficients at order n
 *          WARNING: for alpha < 0, a good precision is achieved ONLY if
 *          the (Fourier) order 0 of the (Taylor) order 0 is >> wrt to the rest of the (Fourier) coefficients of the (Taylor) order 0.
 **/
template<typename T> template<typename U> Ofts<T>& Ofts<T>::tfts_pows(Ofts<T> const& a,  U const& alpha, int const& n)
{
    cout << "tfts_pows is an empty shell if T != Ofsc" << endl;
    return *this;
}

/**
 *   \brief Power function: p = a^alpha, with Tfs coefficients at order n. Double complex version.
 *          WARNING: for alpha < 0, a good precision is achieved ONLY if
 *          the (Fourier) order 0 of the (Taylor) order 0 is >> wrt to the rest of the (Fourier) coefficients of the (Taylor) order 0.
 **/
template<> template<typename U> Ofts< Ofsc >& Ofts< Ofsc >::tfts_pows(Ofts< Ofsc > const& a,   U const& alpha, int const& n)
{
    if(n == 0)
    {
        //Initialization of order zero
        coef0s(*this)->tfs_pows(*coef0s(a), alpha);
    }
    else
    {
        Ofsc temp(corder);
        //temp2 = 1/x0
        temp.tfs_pows(*coef0s(a), -1.0+0.0*I);
        //Sets every coefficients to zero for order n
        this->m_ofts_term[n]->zero();
        //Loop on previously computed coefficient
        for(int j=0; j<= n-1; j++)
        {
            this->m_ofts_term[n]->tftsh_smprod_tu(*a.m_ofts_term[n-j], *this->m_ofts_term[j], temp, (alpha*(n-j)-j)/n);
        }

    }

    return *this;
}

/**
 *   \brief Compute the order zero of the power function: p = a^alpha, with Tfs coefficients at order n
 *          WARNING: for alpha < 0, a good precision is achieved ONLY if
 *          the (Fourier) order 0 of the (Taylor) order 0 is >> wrt to the rest of the (Fourier) coefficients of the (Taylor) order 0.
 **/
template<typename T> template<typename U> Ofts<T>& Ofts<T>::tfts_spows_zero(Ofts<T> const& a,  U const& alpha, int const& n)
{
    cout << "tfts_pows_zero is an empty shell if T != Ofsc" << endl;
    return *this;
}

/**
 *   \brief Compute the order zero of the power function: p = a^alpha, with Tfs coefficients at order n. Double complex version.
 *          WARNING: for alpha < 0, a good precision is achieved ONLY if
 *          the (Fourier) order 0 of the (Taylor) order 0 is >> wrt to the rest of the (Fourier) coefficients of the (Taylor) order 0.
 **/
template<> template<typename U> Ofts< Ofsc >& Ofts< Ofsc >::tfts_spows_zero(Ofts< Ofsc > const& a,   U const& alpha, int const& n)
{
    if(n == 0)
    {
        //Initialization of order zero
        coef0s(*this)->tfs_pows(*coef0s(a), alpha);
    }
    else
    {
        Ofsc temp(corder);
        //temp2 = 1/x0
        temp.tfs_pows(*coef0s(a), -1.0+0.0*I);
        //Only order 0
        int j = 0;
        this->m_ofts_term[n]->tftsh_smprod_tu(*a.m_ofts_term[n-j], *this->m_ofts_term[j], temp, (alpha*(n-j)-j)/n);
    }

    return *this;
}



//----------------
// smult
//----------------
/**
 *  \brief  An operation. Adds the product: \c this \f$  += m a \f$ with m a coefficient, at order n
 *  \param  a: a reference to an Oftsh object
 *  \param  m: a reference to a coefficient
 *  \param  n: a reference to the order to update
 */
template<typename T> Ofts<T>& Ofts<T>::tfts_smult_t(Ofts<T> const& a, T const& m, int const& n)
{
    if(m_ofts_order != a.m_ofts_order || m_ofts_nvar != a.m_ofts_nvar)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        m_ofts_term[n]->tftsh_smult_t(*a.m_ofts_term[n], m);  //ps[n] += m*a[n]
        return *this;
    }
}

/**
 *  \brief  An operation. Adds the product: \c this \f$  += c a \f$ with c a subcoefficient at order n.
 *  \param  a: a reference to an Oftsh object
 *  \param  c: a reference to a subcoefficient
 */
template<typename T> template<typename U> Ofts< Ofs<U> >& Ofts<T>::tfts_smult_u(Ofts< Ofs<U> > const& a, U const& c, int const& n)
{
    if(n > m_ofts_order || m_ofts_nvar != a.m_ofts_nvar)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        m_ofts_term[n]->tftsh_smult_u(*a.m_ofts_term[n], c);  //ps[k] += m*a[k]
        return *this;
    }
}

//----------------
// sprod
//----------------
/**
 *  \brief  An operation. Adds the product: \c this \f$  += a*b \f$ at order n.
 *
 *  Handle the case for which n >= max(a.m_ofts_order, b.m_ofts_order)
 */
template<typename T> Ofts<T>& Ofts<T>::tfts_sprod(Ofts<T> const& a, Ofts<T> const& b, int const& n)
{
    int i, i0, i1;
    i0 = min(b.m_ofts_order, n);
    i1 = min(a.m_ofts_order, n);
    //Product
    for(i= n-i0; i<=i1; i++) m_ofts_term[n]->tftsh_sprod(*a.m_ofts_term[i], *b.m_ofts_term[n-i]);
    return *this;
}

/**
 *  \brief  An operation. Adds the order zero of the product: \c this \f$  += a*b \f$ at order n.
 *
 *  Handle the case for which n >= max(a.m_ofts_order, b.m_ofts_order)
 */
template<typename T> Ofts<T>& Ofts<T>::tfts_sprod_zero(Ofts<T> const& a, Ofts<T> const& b, int const& n)
{
    int i0 = min(b.m_ofts_order, n);
    int i1 = min(a.m_ofts_order, n);
    //Product
    m_ofts_term[n]->tftsh_sprod(*a.m_ofts_term[n-i0], *b.m_ofts_term[i0]);
    m_ofts_term[n]->tftsh_sprod(*a.m_ofts_term[i1], *b.m_ofts_term[n-i1]);
    return *this;
}

/**
 *  \brief  An operation. Adds the product: \c this \f$  += m*a*b \f$ with m a coefficient at order n
 *  \param  a: a reference to an Ofts object
 *  \param  b: a reference to an Ofts object
 *  \param  n: a reference to a coefficient
 */
template<typename T> template<typename U> Ofts<T>& Ofts<T>::tfts_smprod_u(Ofts<T> const& a, Ofts<T> const& b, U const& c, int const& n)
{
    //-----------------------------------------------------------------------
    //Other possibility: use a temporary variable in Ofs, and not Ofts
    int i, i0, i1;
    i0 = min(b.m_ofts_order, n);
    i1 = min(a.m_ofts_order, n);
    for(i= n-i0; i<=i1; i++) m_ofts_term[n]->tftsh_smprod_u(*a.m_ofts_term[i], *b.m_ofts_term[n-i], c);

    return *this;
    //-----------------------------------------------------------------------
}


/**
 *  \brief  An operation. Adds the order zero of the product: \c this \f$  += m*a*b \f$ with m a coefficient at order n
 *  \param  a: a reference to an Ofts object
 *  \param  b: a reference to an Ofts object
 *  \param  n: a reference to a coefficient
 */
template<typename T> template<typename U> Ofts<T>& Ofts<T>::tfts_smprod_u_zero(Ofts<T> const& a, Ofts<T> const& b, U const& c, int const& n)
{
    //-----------------------------------------------------------------------
    //Other possibility: use a temporary variable in Ofs, and not Ofts
    int i0 = min(b.m_ofts_order, n);
    int i1 = min(a.m_ofts_order, n);
    m_ofts_term[n]->tftsh_smprod_u(*a.m_ofts_term[n-i0], *b.m_ofts_term[i0], c);
    m_ofts_term[n]->tftsh_smprod_u(*a.m_ofts_term[i1], *b.m_ofts_term[n-i1], c);
    return *this;
    //-----------------------------------------------------------------------
}

//----------------
// sfsum
//----------------
/**
 *  \brief  An operation. Adds the product: \c this \f$  += ma*a + mb*b \f$ with ma and mb coefficients at order n.
 *  \param  a: a reference to an Ofts object
 *  \param  ma: a reference to a coefficient
 *  \param  b: a reference to an Ofts object
 *  \param  mn: a reference to a coefficient
 */
template<typename T> Ofts<T>& Ofts<T>::tfts_sfsum_t(Ofts<T> const& a, T const& ma, Ofts<T> const& b, T const& mb, int n)
{
    if(m_ofts_order != a.m_ofts_order || m_ofts_nvar != a.m_ofts_nvar || m_ofts_order != b.m_ofts_order || m_ofts_nvar != b.m_ofts_nvar)
    {
        cout << "Error using sfsum: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        m_ofts_term[n]->tftsh_smult_t(*a.m_ofts_term[n], ma);   //ps[k]  += ma*a[k]
        m_ofts_term[n]->tftsh_smult_t(*b.m_ofts_term[n], mb);   //ps[k]  += mb*b[k]
        return *this;
    }

};

/**
 *  \brief  An operation. Adds the product: \c this \f$  += ma*a + mb*b + mc*c \f$ with ma and mb coefficients at order n.
 *  \param  a: a reference to an Ofts object
 *  \param  ma: a reference to a coefficient
 *  \param  b: a reference to an Ofts object
 *  \param  mb: a reference to a coefficient
 *  \param  c: a reference to an Ofts object
 *  \param  mc: a reference to a coefficient
 *  \param  n: a reference to the order to update
 */
template<typename T> Ofts<T>& Ofts<T>::tfts_sfsum_tt(Ofts<T> const& a, T const& ma, Ofts<T> const& b, T const& mb, Ofts<T> const& c, T const& mc, int n)
{
    if(m_ofts_order != a.m_ofts_order || m_ofts_nvar != a.m_ofts_nvar || m_ofts_order != b.m_ofts_order || m_ofts_nvar != b.m_ofts_nvar|| m_ofts_order != c.m_ofts_order || m_ofts_nvar != c.m_ofts_nvar)
    {
        cout << "Error using sfsum: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        m_ofts_term[n]->tftsh_smult_t(*a.m_ofts_term[n], ma);   //ps[k]  += ma*a[k]
        m_ofts_term[n]->tftsh_smult_t(*b.m_ofts_term[n], mb);   //ps[k]  += mb*b[k]
        m_ofts_term[n]->tftsh_smult_t(*c.m_ofts_term[n], mc);   //ps[k]  += mb*b[k]
        return *this;
    }

}

/**
 *  \brief  An operation. Sets the product: \c this \f$  = ca*a + cb*b \f$ with ca and cb subcoefficients at order m
 *  \param  a: a reference to an Ofts object
 *  \param  ca: a reference to a subcoefficient
 *  \param  b: a reference to an Ofts object
 *  \param  cb: a reference to a subcoefficient
 *  \param  m: a reference to the order to update
 */
template<typename T> template<typename U> Ofts<T>& Ofts<T>::tfts_fsum_u(Ofts<T> const& a, U const& ca, Ofts<T> const& b, U const& cb, int const& m)
{
    if(m_ofts_order != a.m_ofts_order || m_ofts_nvar != a.m_ofts_nvar || m_ofts_order != b.m_ofts_order || m_ofts_nvar != b.m_ofts_nvar)
    {
        cout << "Error using fsum: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
        return *this;
    }
    else
    {
        m_ofts_term[m]->tftsh_mult_u(*a.m_ofts_term[m], ca);        //ps[k]   = ca*a[k], contains the zeroing
        m_ofts_term[m]->tftsh_smult_u(*b.m_ofts_term[m], cb);       //ps[k]  += mb*b[k]
        return *this;
    }

}

//----------------
// der
//----------------
/**
 *   \brief Partial derivative at order m: works for m_ofts_order = a.m_ofts_order. TFS format.
 *          WARNING: order m is derived and set in order m-1 of this!
 *          If m==0, nothing is done.
 **/
template<typename T> Ofts<T>& Ofts<T>::tfts_der(Ofts< T > const& a, int ni, int m)
{
    if(m==0)
    {
        cout << "Warning in tfts_der @order m (OFTS object). " << endl;
        cout << "Partial derivatives @order m=0 has been required by user. " << endl;
        cout << "Nothing is done. " << endl;

    }
    else m_ofts_term[m-1]->tfts_derh(*a.m_ofts_term[m], ni); //WARNING: order m is derived and set in order m-1 of this!

    return *this;
}

/**
 *   \brief Same as der but the result is added to the current Ofts instance.
 **/
template<typename T> Ofts<T>& Ofts<T>::tfts_sder(Ofts< T > const& a, int ni, int m)
{
    if(m==0)
    {
        cout << "Warning in tfts_sder @order m (OFTS object). " << endl;
        cout << "Partial derivatives @order m=0 has been required by user. " << endl;
        cout << "Nothing is done. " << endl;

    }
    else m_ofts_term[m-1]->tfts_sderh(*a.m_ofts_term[m], ni); //WARNING: order m is derived and set in order m-1 of this!

    return *this;
}


//----------------------------------------------------------------------------------------
//Operations with specific conditions on the Fourier-Taylor series
//----------------------------------------------------------------------------------------
/**
 *  \brief Power function: p = a^alpha at order n, for Oftsd objects with Fourier
 *         coefficients in Frequency domain.
 *         Moreover, the order zero of the FT series must be unitary.
 **/
template<> template<typename U> Ofts< Ofsd >& Ofts< Ofsd >::ofts_pows(Ofts< Ofsd > const& a,   U const& alpha, int const& n)
{
    if(n==0)
    {
        //Initialization of order zero
        this->acoef0s(*coef0s(a));  //a[0]^alpha = 1.0^alpha = 1;
        return *this;
    }
    else
    {
        //Sets every coefficients to zero for m_ofts_order k
        this->m_ofts_term[n]->zero();
        int i0 = min(a.m_ofts_order, n);
        for(int i= n-i0; i< n; i++) this->m_ofts_term[n]->oftsh_smprod_u(*a.m_ofts_term[n-i], *this->m_ofts_term[i], (double) (alpha*(n-i)-i)/n);
        return *this;
    }
}


/**
 *  \brief Power function: p = a^alpha at order n, for Oftsd objects with Fourier
 *         coefficients in Frequency domain.
 *         Moreover, the order zero of the FT series a[0] must satisfy: a[0] >> a[i], for all i > 0
 *         In this routine, the coefficient a0inv = 1/(a[0]) and a0palpha = a[0]^alpha
 **/
template<> template<typename U> Ofts< Ofsd >& Ofts< Ofsd >::pows(Ofts< Ofsd > const& a,  //intitial Ofts
                                                                 Ofsd a0inv,             //inverse of order 0
                                                                 Ofsd a0palpha,          //order 0 ^alpha
                                                                 U const& alpha)         //power coef
{
    //Initialization of order zero
    this->acoef0s(a0palpha);  //a[0]^alpha = a0palpha provided
    //Recurrence scheme
    for(int k=1; k <= m_ofts_order; k++)
    {
        //Sets every coefficients to zero for m_ofts_order k
        this->m_ofts_term[k]->zero();
        for(int j=0; j<= k-1; j++)
        {
            this->m_ofts_term[k]->oftsh_smprod_u(*a.m_ofts_term[k-j], *this->m_ofts_term[j], (double) (alpha*(k-j)-j)/k);
            this->m_ofts_term[k]->oftsh_mult_t(*this->m_ofts_term[k], a0inv);
        }

    }
    return *this;
}


//----------------------------------------------------------------------------------------
//Stream
//----------------------------------------------------------------------------------------
/**
 *   \brief Stream operator << for Ofts objects.
 **/
template<typename T> std::ostream& operator << (std::ostream& stream, Ofts<T> const& ofts)
{
    int i,j, nrc;
    int k[ofts.m_ofts_nvar];

    for(nrc=0; nrc<= ofts.m_ofts_order; nrc++)
    {
        //Update k
        k[0] = nrc;
        for(i=1; i<ofts.m_ofts_nvar; i++) k[i] = 0;
        //Current homogeneous polynomial
        for (i=0; i< Manip::nmon(ofts.m_ofts_nvar, nrc); i++)
        {
            for(j=0; j<ofts.m_ofts_nvar; j++) stream <<   setiosflags(ios::right) <<  k[j] << " ";
            stream << endl;
            stream << std::showpos << setiosflags(ios::scientific)  << setprecision(16)  <<  ofts.m_ofts_term[nrc]->get_coef(i) << std::noshowpos << endl;

            if(i< Manip::nmon(ofts.m_ofts_nvar, nrc)-1)  Manip::prxkt(k, ofts.m_ofts_nvar);
        }
    }
    return stream;
}

/**
 *   \brief Stream operator for Ofts objects. Print only the order zero of each coefficients.
 **/
template<typename T> void Ofts<T>::fprint_0(ofstream& stream)
{
    int i,j, nrc;
    int k[m_ofts_nvar];
    for(nrc=0; nrc<= m_ofts_order; nrc++)
    {
        //Update k
        k[0] = nrc;
        for(i=1; i<m_ofts_nvar; i++) k[i] = 0;
        //Current homogeneous polynomial
        for (i=0; i< Manip::nmon(m_ofts_nvar, nrc); i++)
        {
            //for(j=0; j<m_ofts_nvar; j++) stream <<  setw(2) << setiosflags(ios::right) <<  k[j] << " ";
            for(j=0; j<m_ofts_nvar; j++) stream <<   setiosflags(ios::right) <<  k[j] << " ";
            stream << endl;
            stream << std::showpos << setiosflags(ios::scientific)  << setprecision(15);
            m_ofts_term[nrc]->get_coef(i).fprint_0(stream);
            stream << std::noshowpos << endl;
            if(i< Manip::nmon(m_ofts_nvar, nrc)-1)  Manip::prxkt(k, m_ofts_nvar);
        }
    }
}



//----------------------------------------------------------------------------------------
//Evaluate
//----------------------------------------------------------------------------------------
/**
 *  \brief Generic routine for the evaluation of an Ofts object at the state X: z = this(X).
 **/
template<typename T> template<typename U> void Ofts<T>::evaluate(U X[], T& z)
{
    z.zero();
    for(int k = m_ofts_order; k >= 0 ; k--)
    {
        m_ofts_term[k]->sevaluate(X, z);
    }
}

/**
 *  \brief Generic routine for the conjugate evaluation of an Ofts object at the state X: z = conj(this(X)).
 **/
template<typename T> template<typename U> void Ofts<T>::evaluate_conjugate(U X[], T& z)
{
    z.zero();
    for(int k = m_ofts_order; k >= 0 ; k--)
    {
        m_ofts_term[k]->sevaluate_conjugate(X, z);
    }
}

/**
 *  \brief Toutine for the evaluation of an Ofts object at the state X, at order m: z = [this(X)]_m.
 **/
template<typename T> template<typename U> void Ofts<T>::evaluate(U X[], T& z, int const& m, int const& ofs_order)
{
    z.zero();
    for(int k = m; k >= 0 ; k--)
    {
        m_ofts_term[k]->sevaluate(X, z, ofs_order);
    }
}

/**
 *  \brief Routine for the evaluation of an Ofts object at the state X, at order m, and at time t: z = [this(X, t)]_(<=m).
 **/
template<typename T> template<typename U> cdouble Ofts<T>::fevaluate(U X[], double const& t, int const& m, int const& ofs_order)
{
    //Initialize the state
    cdouble z = 0.0+0.0*I;

    //Initialize the cosinus/sinus arrays
    double cR[ofs_order];
    double sR[ofs_order];

    cR[0] = cos(t);
    sR[0] = sin(t);
    for(int i = 1; i< ofs_order; i++)
    {
        cR[i] =  cR[i-1]*cR[0] - sR[i-1]*sR[0];
        sR[i] =  sR[i-1]*cR[0] + cR[i-1]*sR[0];
    }

    //Initialize the exponents array
    int kv[m_ofts_nvar];

    //Loop on all desired orders
    for(int k = m; k >= 0 ; k--)
    {
        m_ofts_term[k]->fevaluate(X, z, kv, cR, sR, ofs_order);
    }

    return z;
}

/**
 *  \brief Routine for the evaluation of an Ofts object at the state X, at order m, and at time t: z = [this(X, t)]_(<=m).
 *
 *         Contrary to the routine fevaluate(U X[], double const& t, int const& m, int const& ofs_order), the cosinus/sinus arrays are given as inputs:
 *              - cR[] = [cos(t), ..., cos(ofs_order*t)]
 *              - sR[] = [sin(t), ..., sin(ofs_order*t)]
 **/
template<typename T> template<typename U> cdouble Ofts<T>::fevaluate(U X[], double cR[], double sR[], int const& m, int const& ofs_order)
{
    //Initialize the state
    cdouble z = 0.0+0.0*I;

    //Initialize the exponents array
    int kv[m_ofts_nvar];

    //Loop on all desired orders
    for(int k = m; k >= 0 ; k--)
    {
        m_ofts_term[k]->fevaluate(X, z, kv, cR, sR, ofs_order);
    }

    return z;
}

/**
 *  \brief Contribution of the order m of this to the evaluation: z = [this(X)]_=m).
 **/
template<typename T> template<typename U> void Ofts<T>::contribution(U X[], T& z, int const& m)
{
    z.zero();
    m_ofts_term[m]->sevaluate(X, z);
}


//----------------------------------------------------------------------------------------
//Derivation
//----------------------------------------------------------------------------------------
/**
 *  \brief Partial derivative wrt to the variable z[ni], with ni = 1,...n: \c this \f$ = \frac{\partial a}{\partial z_{n_i}} \f$.
 **/
template<typename T> Ofts<T>& Ofts<T>::der(Ofts<T> const& a, int ni)
{
    for(int k = 0; k < m_ofts_order ; k++) //Careful here: the sum goes up to m_ofts_order-1!
    {
        m_ofts_term[k]->derh(*a.m_ofts_term[k+1], ni);
    }
    return *this;
}

/**
 *  \brief Partial derivative wrt to the variable z[ni], with ni = 1,...n: \c this \f$ += \frac{\partial a}{\partial z_{n_i}} \f$.
 **/
template<typename T> Ofts<T>& Ofts<T>::sder(Ofts< T > const& a, int ni)
{
    for(int k = 0; k < m_ofts_order ; k++) //Careful here: the sum goes up to m_ofts_order-1!
    {
        m_ofts_term[k]->sderh(*a.m_ofts_term[k+1], ni);
    }
    return *this;
}

/**
 *  \brief Partial derivative wrt to the variable z[ni], with ni = 1,...n, at order m of the expansions:
 *          \f$ this_{m-1}  = \left[ \frac{\partial a}{\partial z_{n_i}} \right]_m \f$.
 *
 *          Note: order m is derived and set in order m-1 of this. If m==0, nothing is done, but a warning is sent to the user.
 **/
template<typename T> Ofts<T>& Ofts<T>::der(Ofts< T > const& a, int ni, int m)
{
    if(m==0)
    {
        cout << "Warning in der @order m (OFTS object). " << endl;
        cout << "Partial derivatives @order m=0 has been required by user. " << endl;
        cout << "Nothing is done. " << endl;

    }
    else m_ofts_term[m-1]->derh(*a.m_ofts_term[m], ni); //WARNING: order m is derived and set in order m-1 of this!

    return *this;
}

/**
 *  \brief Partial derivative wrt to the variable z[ni], with ni = 1,...n, at order m of the expansions:
 *          \f$ this_{m-1}  += \left[ \frac{\partial a}{\partial z_{n_i}} \right]_m \f$.
 *
 *          Note: order m is derived and set in order m-1 of this. If m==0, nothing is done, but a warning is sent to the user.
 **/
template<typename T> Ofts<T>& Ofts<T>::sder(Ofts< T > const& a, int ni, int m)
{
    if(m==0)
    {
        cout << "Warning in der @order m (OFTS object). " << endl;
        cout << "Partial derivatives @order m=0 has been required by user. " << endl;
        cout << "Nothing is done. " << endl;

    }
    else m_ofts_term[m-1]->sderh(*a.m_ofts_term[m], ni); //WARNING: order m is derived and set in order m-1 of this!

    return *this;
}

/**
 *  \brief Partial derivative wrt to time:
 *          \f$ this  += \left[ \frac{\partial a}{\partial t} \right] \f$.
 **/
template<typename T> Ofts<T>& Ofts<T>::dot(Ofts<T> const& a, double const&  n)
{
    if(m_ofts_order != a.m_ofts_order || m_ofts_nvar != a.m_ofts_nvar)
    {
        cout << "Error using dot: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
    }
    else
    {
        //Loop
        for(int k = 0; k <= m_ofts_order ; k++) m_ofts_term[k]->dot(*a.m_ofts_term[k], n);
    }
    return *this;
}

/**
 *  \brief Partial derivative wrt to the time, at order k of the expansions:
 *          \f$ this_{k-1}  += \left[ \frac{\partial a}{\partial t} \right]_k \f$.
 **/
template<typename T> Ofts<T>& Ofts<T>::dot(Ofts<T> const& a, double const&  n, int const& k)
{
    if(m_ofts_order != a.m_ofts_order || m_ofts_nvar != a.m_ofts_nvar)
    {
        cout << "Error using dot: the order and/or number of variables does not match. Initial Ofts<T> is returned" << endl;
    }
    else
    {
        m_ofts_term[k]->dot(*a.m_ofts_term[k], n);
    }
    return *this;
}

//----------------------------------------------------------------------------------------
//Integral
//----------------------------------------------------------------------------------------
/**
 *  \brief Primitive wrt to the variable z[ni], with ni = 1,...n: \c this \f$ = \int a dz_{n_i} \f$.
 **/
template<typename T> Ofts<T>& Ofts<T>::sprim(Ofts< T > const& a, int ni)
{
    for(int k = 0; k < m_ofts_order ; k++) //Careful here: the sum goes up to m_ofts_order-1!
    {
        m_ofts_term[k+1]->sprimh(*a.m_ofts_term[k], ni);
    }
    return *this;
}

//----------------------------------------------------------------------------------------
//Norms
//----------------------------------------------------------------------------------------
/**
 *  \brief L1 norm of the term of order m: returns \f$ L_1 \left( [this]_m \right) \f$
 **/
template<typename T> double Ofts<T>::l1norm(int const& m)
{
    return m_ofts_term[m]->l1norm();
}

/**
 *  \brief Infinity norm of the term of order m: returns \f$ L_\infty \left( [this]_m \right) \f$
 **/
template<typename T> double Ofts<T>::linfnorm(int const& m)
{
    return m_ofts_term[m]->linfnorm();
}

/**
 *  \brief  Number of small divisors under a certain value sdmax in the term of order m
 */
template<typename T> int Ofts<T>::nsd(int const& m, int odmax, double sdmax)
{
    return m_ofts_term[m]->nsd(odmax, sdmax);
}


//----------------------------------------------------------------------------------------
//          Reading & writing (I/O)
//----------------------------------------------------------------------------------------
//----------------------------------------------
// Text format, write
//----------------------------------------------
/**
 * \brief Writes a given object W of type \c Ofts<Ofsc >  in a txt files of the form "filename".
 **/
inline void  writeOFTS_txt(Ofts<Ofsc > &W, string filename)
{
    ofstream myfile;
    myfile.open ((filename).c_str(), ios::out);
    myfile << W << endl;
    myfile.close();
}


/**
 * \brief Writes a given vector W of type \c Ofts<Ofsc >  in a txt files of the form "filename+i.txt", with i = 0, length(W)-1.
 **/
inline void  write_vofts_txt(vector<Ofts<Ofsc > > &W, string filename)
{
    ofstream myfile;
    string ss1;
    //Loop on all coefficients
    for(unsigned int i = 0; i < W.size(); i++)
    {
        ss1 = static_cast<ostringstream*>( &(ostringstream() << i) )->str();
        myfile.open ((filename+"["+ss1+"].txt").c_str(), ios::out);
        myfile << W[i] << endl;
        myfile.close();
    }
}

//----------------------------------------------
// Text format, read
//----------------------------------------------
/**
 * \brief Reads a given \c Ofsc  object within a \c Ofts<Ofsc >, in txt format.
 **/
inline void read_ofs_txt(Ofsc &xFFT, ifstream &readStream, int fftN)
{
    //Init
    double ct, cr, ci;
    //Reading
    for(int i = -fftN; i<=fftN; i++)
    {
        readStream >> ct;  //current m_ofts_order
        readStream >> cr;  //real part
        readStream >> ci;  //imag part
        xFFT.set_coef(cr+I*ci, i);
    }
}

/**
 * \brief Reads a given \c Ofts<Ofsc >  object, in txt format.
 **/
inline int read_ofts_txt(Ofts<Ofsc > &x, string filename, int fftN)
{
    //Init
    ifstream readStream;
    string ct;
    //Reading
    readStream.open((filename).c_str(), ios::in);

    //Check that the opening went well
    if (!readStream.is_open())
    {
        cout << "read_ofts_txt. Cannot open file " << filename << endl;
        cout << "Check the text data exist." << endl;
        return GSL_FAILURE;
    }
    else
    {
        for(int k = 0 ; k <= x.get_order() ; k++)
        {
            for(int p = 0; p < Manip::nmon(x.get_nvar(), k); p++)
            {
                //Current kv
                getline(readStream, ct);
                //Reading the coefficient
                read_ofs_txt(*x.get_coef(k,p), readStream, fftN);
                getline(readStream, ct);
                getline(readStream, ct);
            }
        }
        readStream.close();
        return GSL_SUCCESS;
    }
}

/**
 * \brief Reads a given vector W of type \c Ofts<Ofsc >  in a txt files of the form "filename+i.txt", with i = 0, length(W)-1.
 **/
inline void read_vofts_txt(vector<Ofts<Ofsc > >  &W, string filename, int fftN)
{
    string ss1;
    //Loop on all coefficients
    for(unsigned int i = 0; i < W.size(); i++)
    {
        ss1 = static_cast<ostringstream*>( &(ostringstream() << i) )->str();
        read_ofts_txt(W[i], (filename+"["+ss1+"].txt"), fftN);
    }
}


//----------------------------------------------
// Binary format, write
//----------------------------------------------
/**
 * \brief Writes a given \c Ofsc  object within a \c Ofts<Ofsc >, in bin format.
 **/
inline void  write_ofs_bin(Ofsc  &xFFT, fstream &myfile)
{
    //Init
    int fftN = xFFT.get_order();
    double res;

    //Writing
    for(int i = -fftN; i<=fftN; i++)
    {
        //Real part
        res = creal(xFFT.ofs_get_coef(i));
        myfile.write((char*) &res, sizeof(double));

        //Imag part
        res = cimag(xFFT.ofs_get_coef(i));
        myfile.write((char*) &res, sizeof(double));
    }
}

/**
 * \brief Writes a given \c Ofts<Ofsc >  object, in bin format.
 **/
inline void  write_ofts_bin(Ofts<Ofsc > &W, string filename)
{
    fstream myfile;
    myfile.open (filename.c_str(), ios::binary | ios::out);
    //Loop on order
    for(int nrc=0; nrc<= W.get_order(); nrc++)
    {
        //Current homogeneous polynomial
        for (int i=0; i< Manip::nmon(W.get_nvar(), nrc); i++)
        {
            //Write each Ofs coefficient
            write_ofs_bin(*W.get_coef(nrc,i), myfile);
        }
    }
    myfile.close();
}

/**
 * \brief Writes a given vector W of type \c Ofts<Ofsc >  in a binary files of the form "filename+i.bin", with i = 0, length(W)-1.
 **/
inline void  write_vofts_bin(vector<Ofts<Ofsc > > &W, string filename)
{
    string ss1;
    //Loop on all coefficients
    for(unsigned int i = 0; i < W.size(); i++)
    {
        ss1 = static_cast<ostringstream*>( &(ostringstream() << i) )->str();
        write_ofts_bin(W[i], (filename+"["+ss1+"].bin"));
    }
}

//----------------------------------------------
// Binary format, read
//----------------------------------------------
/**
 * \brief Reads a given \c Ofsc  object within a \c Ofts<Ofsc >, in bin format.
 **/
inline void read_ofs_bin(Ofsc  &xFFT, fstream &myfile)
{
    int fftN = xFFT.get_order();
    double cr, ci;
    //Writing
    for(int i = -fftN; i<=fftN; i++)
    {
        //Real part
        myfile.read((char*)&cr, sizeof(double));
        //Imag part
        myfile.read((char*)&ci, sizeof(double));
        //Put in current position
        xFFT.set_coef(cr+I*ci, i);
    }
}

/**
 * \brief Reads a given \c Ofts<Ofsc >  object, in bin format.
 **/
inline int read_ofts_bin(Ofts<Ofsc > &W, string filename)
{
    //Init
    fstream myfile;

    //Open the stream
    myfile.open((filename).c_str(), ios::binary | ios::in);

    //Check that the opening went well
    if (!myfile.is_open())
    {
        cout << "read_ofts_bin. Cannot open file " << filename << endl;
        cout << "read_ofts_bin. Check the binary data exist." << endl;
        return GSL_FAILURE;
    }
    else
    {
        //Loop on order
        for(int nrc=0; nrc<= W.get_order(); nrc++)
        {
            //Current homogeneous polynomial
            for (int i=0; i< Manip::nmon(W.get_nvar(), nrc); i++)
            {
                //Read each Ofs coefficient
                read_ofs_bin(*W.get_coef(nrc,i), myfile);
            }
        }
        myfile.close();
        return GSL_SUCCESS;
    }
}

/**
 * \brief Reads a given vector W of type \c Ofts<Ofsc >  in binary files of the form "filename+i.bin", with i = 0, length(W)-1.
 **/
inline void read_vofts_bin(vector<Ofts<Ofsc > >  &W, string filename, int fftN)
{
    string ss1;
    int status, global_status;
    global_status = 0;

    //Loop on all coefficients
    for(unsigned int i = 0; i < W.size(); i++)
    {
        ss1 = static_cast<ostringstream*>( &(ostringstream() << i) )->str();

        //Read binary format
        status = read_ofts_bin(W[i], (filename+"["+ss1+"].bin"));

        //Try txt format if failure
        if(status != GSL_SUCCESS)
        {
            cout << "read_vofts_bin. Last reading went wrong. Trying to find data in txt format..." << endl;
            status = read_ofts_txt(W[i], (filename+"["+ss1+"].txt"), fftN);
            if(status != GSL_SUCCESS)
            {
                cout << "read_vofts_bin. Txt format also went wrong. Check data manually." << endl;
                global_status--;
            }
            else
            {
                cout << "read_vofts_bin. Success with txt format." << endl;
                global_status++;
            }
        }
    }

    //If global success, we can save in binary format
    if(global_status == (int) W.size())
    {
        int ch;
        cout << "read_vofts_bin. Success with txt format for all components." << endl;
        cout << "Please enter 1 if you want to save the vector in binary files:" << endl;
        scanf("%d",&ch);
        if(ch == 1) write_vofts_bin(W, filename);
    }

}

//----------------------------------------------------------------------------------------
// Binary format, copy in lesser dimension series
//----------------------------------------------------------------------------------------
/**
 * \brief Reads a given \c Oftsc  object, in bin format, and transform it into another Oftsc object, of lesser dimensions
 *        More precisely:
 *        We suppose that W is a very sparse Fourier-Taylor series, with only non-null coefficients along the dimension dim.
 *        Hence, it is possible to entirely store W into a simpler series W1, with only one dimension.
 **/
inline int fromOFTStoOFTS_bin(Oftsc &W, Oftsc &W1, int dim)
{
    //====================================================================================
    //Init & check
    //====================================================================================
    int m_ofts_nvar  = W.get_nvar();  ///number of variables
    //Check that the desired dimension dim is consistent with m_ofts_nvar
    if(dim > m_ofts_nvar-1)
    {
        cout << "fromOFTStoOFTS_bin. dim is greater than the number of dimensions." << endl;
        return -1;
    }

    //====================================================================================
    //Open the stream & check
    //====================================================================================
        //Parameters
        int *kv = (int*) calloc(m_ofts_nvar, sizeof(int)); //exponent

        //================================================================================
        //Loop on order
        //================================================================================
        for(int nrc=0; nrc<= W.get_order(); nrc++)
        {
            //kv = (k 0 0 0 ...)
            kv[0] = nrc;
            for(int i=1; i<m_ofts_nvar; i++) kv[i] = 0;


            //============================================================================
            //Loop on the monomials at order nrc
            //============================================================================
            for (int i=0; i< Manip::nmon(m_ofts_nvar, nrc); i++)
            {
                //If the exponents are non null only on dim, we save the value in W1
                W1.get_coef(nrc, 0)->ccopy(*W.get_coef(nrc,i));

                //Update the exponents
                if(i< Manip::nmon(m_ofts_nvar, nrc)-1)  Manip::prxkt(kv, m_ofts_nvar);
            }

    }

    return 0;
}

/**
 * \brief Reads some Oftsc  objects, in bin format,
 *        and transform it into other Oftsc objects, of lesser dimensions.
 *
 *        Indeed, we know that, if the center-stable or center-unstable manifolds are
 *        computed using the parameterization method, the series in the dimensions
 *        0, and 3 of Wh (parameterization in TFC coordinates) are of the form:
 *
 *          Wh[0] = sum_(k=0)^n c_k s_5^k
 *
 *        I.e. they are FT series of only one variable, namely the variable s5,
 *        last variable of the reduced variables (s1, s2, s3, s4, s5).
 *
 *        To limit the amount of memory needed to stored these series, we can use
 *        one-dimensional FT series. This is done via the present routine.
 *
 *        This routine makes use of fromOFTStoOFTS_bin to read and store into the less dimensions objects.
 *        Then, the results are stored in binary files.
 **/
inline void from_vofts_to_vofts_bin(vector<Oftsc> &W, Oftsc &W1, Oftsc &DW1, string fileout)
{
    //====================================================================================
    //Init
    //====================================================================================
    string ss1;
    int ind = 0;

    //====================================================================================
    //Loop on all coefficients
    //====================================================================================
    for(unsigned int i = 0; i < 6; i++)
    {
        if(i == 0 || i == 3) //only along certain dimensions
        {
            //Store in 1-dim series
            fromOFTStoOFTS_bin(W[i], W1, 4);

            //ss1 = numToString(ind)
            ss1 = static_cast<ostringstream*>( &(ostringstream() << ind) )->str();

            //Store the result
            write_ofts_bin(W1, (fileout+"["+ss1+"].bin"));
            //writeOFTS_txt(W1, (fileout+"["+ss1+"].txt"));

            //Jacobian
            for(int m = 1; m <= W1.get_order(); m++) DW1.der(W1, 1, m);

            //Store the Jacobian
            write_ofts_bin(DW1, (fileout+"d["+ss1+"].bin"));
            //writeOFTS_txt(DW1, (fileout+"d["+ss1+"].txt"));

            //Advance ind
            ind++;
        }
    }
}

//----------------------------------------------------------------------------------------
// Text 2 binary
//----------------------------------------------------------------------------------------
/**
 * \brief Reads a given vector W of type \c Ofts<Ofsc >  in a txt files of the form "filename+i.txt", with i = 0, length(W)-1.
 *        Writes it again in binary form, in binary files of the form "filename+i.bin", with i = 0, length(W)-1.
 **/
inline void txt2bin_vofts(vector<Ofts<Ofsc > > &W, string filename, int fftN)
{
    //Read from txt files
    read_vofts_txt(W, filename, fftN);
    //Write in binary format
    write_vofts_bin(W, filename);
}


/**
 *  \brief Transform a one-variable Ofts< Ofs<U> > fts_z into an Ofs<U> object fs object such that fts_z(epsilon) = fs.
 **/
template <typename U>  void fts2fs(Ofs<U> *fs, Ofts< Ofs<U> > const& fts_z, double epsilon)
{

    int k, l;
    double tempc;
    int n_order_fourier = fts_z.get_coef_order();
    Ofs<U> tempfs(n_order_fourier);
    //Cleaning of fs
    fs->zero();
    //Orders >= 0
    for(k = -n_order_fourier; k <= n_order_fourier; k++)
    {
        tempc = 0;
        for(l = 0 ; l <= fts_z.get_order(); l++)
        {
            //Trigonometric coefficient of m_ofts_order l
            tempfs = fts_z.get_term(l)->get_coef(0);
            //tempc += ulj*epsilon^l
            tempc += tempfs.ofs_get_coef(k)*pow(epsilon, (double) l);
        }
        fs->set_coef(tempc, k);
    }
}

/**
 *  \brief Transform a one-variable Ofts< Ofs<U> > fts_z into a complex number fts_z(epsilon, t).
 **/
template <typename U>  cdouble fts2scalar(Ofts< Ofs<U> > const& fts_z, double epsilon, double t)
{
    Ofs<U> fs(fts_z.get_coef_order());
    fts2fs(&fs, fts_z, epsilon);
    return fs.evaluate(t);
}
