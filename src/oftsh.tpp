//########################################################################################
// Implementation of the Oftsh template class
//########################################################################################
/**
 * \file oftsh.tpp
 * \brief Homogeneous Fourier-Taylor series template class. See oftsh.h for details.
 * \author BLB
 */


//----------------------------------------------------------------------------------------
//Create
//----------------------------------------------------------------------------------------
/**
 *  \brief Default constructor of the class Oftsh<T>.
 */
template<typename T> Oftsh<T>::Oftsh()
{
    m_oftsh_order = 0;
    m_oftsh_nvar  = 0;
    m_oftsh_coef  = new T(0);
}

/**
 *  \brief Constructor with given order and number of variables.
 *
 * Allocates memory  without any link to a coefficient array (requires the inmediate used of link_coefs afterwards).
 */
template<typename T> Oftsh<T>::Oftsh(const int t_oftsh_nvar, const int t_oftsh_order)
{
    m_oftsh_order = t_oftsh_order;
    m_oftsh_nvar  = t_oftsh_nvar;
    m_oftsh_coef  = 0;

    //Tree
    if(m_oftsh_nvar==1)  //1-variable polynomial
    {
        m_oftsh_term = (Oftsh<T>*) calloc(m_oftsh_order+1, sizeof(Oftsh<T>)); //calloc is necessary if we are initializing a fourier-taylor series, with taylor-series coefficients
        if (m_oftsh_term == NULL)
        {
            puts("Oftsh<T>::Oftsh<T>: out of memory (2)");
            exit(1);
        }
    }
    else
    {
        int i;
        //Again, calloc is necessary if we are initializing a fourier-taylor series, with taylor-series coefficients
        //At this step, all sons of the current Oftsh<T> have m_oftsh_nvar = 0, m_oftsh_order = 0 and m_oftsh_coef points @ one T.
        m_oftsh_term = (Oftsh<T>*) calloc(m_oftsh_order+1, sizeof(Oftsh<T>));
        if (m_oftsh_term == NULL)
        {
            puts("Oftsh<T>::Oftsh<T>: out of memory (2)");
            exit(1);
        }

        for(i=0; i<=m_oftsh_order; i++)
        {
            //m_oftsh_term[i].lcopy(Oftsh<T>(t_oftsh_nvar-1, t_oftsh_order-i));     //WHY is it not working here but in the other allocation routines ??
            m_oftsh_term[i].lcopy(*(new Oftsh<T>(t_oftsh_nvar-1, t_oftsh_order-i)));
            //At this step, all sons of the current Oftsh<T> have m_oftsh_nvar = m_oftsh_nvar-1, m_oftsh_order = m_oftsh_order-i and m_oftsh_coef points @ one T.
            //As a consequence, a proper linkage is required to ensured that m_oftsh_coef points at a complete array of coefficients
        }
    }
}

/**
 *  \brief Constructor from a given Oftsh object (without any link).
 */
template<typename T> Oftsh<T>::Oftsh(Oftsh<T> const& b)
{
//-----------------------------------------------------------
// The code between the comment lines is replacing the
// code : this->lcopy(Oftsh<T>(b.m_oftsh_nvar, b.m_oftsh_order));
// Goal: avoid the creation of a temporary Oftsh
//-----------------------------------------------------------
    m_oftsh_order = b.m_oftsh_order;
    m_oftsh_nvar  = b.m_oftsh_nvar;
    m_oftsh_coef  = new T(0);//[ Manip::nmon(m_oftsh_nvar, m_oftsh_order)]();
    //Tree
    if(m_oftsh_nvar==1)  //1-variable polynomial
    {
        m_oftsh_term = (Oftsh<T>*) calloc(m_oftsh_order+1, sizeof(Oftsh<T>));
        if (m_oftsh_term == NULL)
        {
            puts("Oftsh<T>::Oftsh<T>: out of memory (2)");
            exit(1);
        }
    }
    else
    {
        int i;
        m_oftsh_term = (Oftsh<T>*) calloc(m_oftsh_order+1, sizeof(Oftsh<T>));  //At this step, all sons of the current Oftsh<T> have m_oftsh_nvar = 0, m_oftsh_order = 0 and m_oftsh_coef points @ one T.
        if (m_oftsh_term == NULL)
        {
            puts("Oftsh<T>::Oftsh<T>: out of memory (2)");
            exit(1);
        }

        for(i=0; i<=m_oftsh_order; i++)
        {
            m_oftsh_term[i].lcopy(Oftsh<T>(m_oftsh_nvar-1, m_oftsh_order-i));
            //At this step, all sons of the current Oftsh<T> have m_oftsh_nvar = m_oftsh_nvar-1, m_oftsh_order = m_oftsh_order-i and m_oftsh_coef points @ one T.
            //As a consequence, a proper linkage is required to ensured that m_oftsh_coef points at a complete array of coefficients
        }
    }
//-----------------------------------------------------------
// End of replacement code
//-----------------------------------------------------------

    //Copy in linking of the coefficients
    T *coef0 = new T[ Manip::nmon(b.m_oftsh_nvar, b.m_oftsh_order)]();
    for(int i = 0 ; i <  Manip::nmon(b.m_oftsh_nvar, b.m_oftsh_order) ; i++) coef0[i] = b.m_oftsh_coef[i];
    this->link_coefs(coef0);
}

/**
 *  \brief  An operator. Constructor from a given Oftsh object (only the coefficients).
 */
template<typename T> Oftsh<T>& Oftsh<T>::operator = (Oftsh<T> const& b)
{
    if(this != &b)
    {
        delete m_oftsh_term;
        delete m_oftsh_coef;
//-----------------------------------------------------------
// The code between the comment lines is replacing the
// code : this->lcopy(Oftsh<T>(b.m_oftsh_nvar, b.m_oftsh_order));
// Goal: avoid the creation of a temporary Oftsh
//-----------------------------------------------------------
        m_oftsh_order = b.m_oftsh_order;
        m_oftsh_nvar = b.m_oftsh_nvar;
        m_oftsh_coef = new T(0);//[ Manip::nmon(m_oftsh_nvar, m_oftsh_order)]();
        //Tree
        if(m_oftsh_nvar==1)  //1-variable polynomial
        {
            m_oftsh_term = (Oftsh<T>*) calloc(m_oftsh_order+1, sizeof(Oftsh<T>));
            if (m_oftsh_term == NULL)
            {
                puts("Oftsh<T>::Oftsh<T>: out of memory (2)");
                exit(1);
            }
        }
        else
        {
            int i;
            m_oftsh_term = (Oftsh<T>*) calloc(m_oftsh_order+1, sizeof(Oftsh<T>));  //At this step, all sons of the current Oftsh<T> have m_oftsh_nvar = 0, m_oftsh_order = 0 and m_oftsh_coef points @ one T.
            if (m_oftsh_term == NULL)
            {
                puts("Oftsh<T>::Oftsh<T>: out of memory (2)");
                exit(1);
            }

            for(i=0; i<=m_oftsh_order; i++)
            {
                m_oftsh_term[i].lcopy(Oftsh<T>(m_oftsh_nvar-1, m_oftsh_order-i));
                //At this step, all sons of the current Oftsh<T> have m_oftsh_nvar = m_oftsh_nvar-1, m_oftsh_order = m_oftsh_order-i and m_oftsh_coef points @ one T.
                //As a consequence, a proper linkage is required to ensured that m_oftsh_coef points at a complete array of coefficients
            }
        }
//-----------------------------------------------------------
// End of replacement code
//-----------------------------------------------------------
        //Copy in linking of the coefficients
        T *coef0 = new T[ Manip::nmon(b.m_oftsh_nvar, b.m_oftsh_order)]();
        for(int i = 0 ; i <  Manip::nmon(b.m_oftsh_nvar, b.m_oftsh_order) ; i++) coef0[i] = b.m_oftsh_coef[i];
        this->link_coefs(coef0);
    }
    return *this; //same object if returned
}


//----------------------------------------------------------------------------------------
//Delete
//----------------------------------------------------------------------------------------
/**
 *  \brief Default destructor of the class Oftsh<T>.
 *         TO BE DETERMINED: how properly delete with recursivity?
 *         Seems to work fine like this, but memory leak?
 */
template<typename T> Oftsh<T>::~Oftsh<T>()
{
    delete m_oftsh_term;
// COMMENTED FOR NOW
//    if(m_oftsh_coef != 0 && m_oftsh_coef != NULL)
//    {
//        delete m_oftsh_coef;
//        m_oftsh_coef = 0;
//    }
}

//----------------------------------------------------------------------------------------
//Copy
//----------------------------------------------------------------------------------------
/**
 *  \brief  Linked copy from a given Oftsh object (exact same object is obtained).
 *  \param  b: a reference to the Oftsh object to copy
 *  \return a reference to the current object
 *
 *  Note: restricted to same order, same number of variables.
 */
template<typename T> Oftsh<T>& Oftsh<T>::lcopy (Oftsh<T> const& b)
{
    m_oftsh_order = b.m_oftsh_order;
    m_oftsh_nvar = b.m_oftsh_nvar;
    m_oftsh_term = b.m_oftsh_term;
    m_oftsh_coef = b.m_oftsh_coef;
    return *this;
}

/**
 *  \brief  Copy from a given Oftsh object (only the coefficients).
 *  \param  b: a reference to the Oftsh object to copy
 *  \return a reference to the current object
 */
template<typename T> Oftsh<T>& Oftsh<T>::ccopy (Oftsh<T> const& b)
{
    if(m_oftsh_order != b.m_oftsh_order || m_oftsh_nvar != b.m_oftsh_nvar)
    {
        cout << "Error using ccopy : the order and/or number of variables does not match. Initial Oftsh<T> is returned" << endl;
        return *this;
    }
    else
    {
        for(int i = 0 ; i< Manip::nmon(m_oftsh_nvar, m_oftsh_order) ; i++) m_oftsh_coef[i] = b.m_oftsh_coef[i];
        return *this;
    }
}


//----------------------------------------------------------------------------------------
//Linking
//----------------------------------------------------------------------------------------
/**
 *  \brief  Performs linking between the Oftsh<T> and an array of coefficients.
 */
template<typename T> void Oftsh<T>::link_coefs(T *coef0)
{
    if(m_oftsh_nvar==1) //1-variable homogeneous polynomial
    {
        m_oftsh_coef = coef0;
    }
    else
    {
        T *coefc;
        int i;
        //allocate the position of the first coefficient
        m_oftsh_coef = coef0;
        //perform recursice allocation
        for(i=0, coefc=coef0; i<=m_oftsh_order; coefc+= Manip::nmon(m_oftsh_nvar-1,m_oftsh_order-i), i++)
        {
            m_oftsh_term[i].link_coefs(coefc);
        }
    }
}


//----------------------------------------------------------------------------------------
//Setters
//----------------------------------------------------------------------------------------
/**
 *  \brief Sets a coefficient at a given position in the polynomial.
 */
template<typename T> void Oftsh<T>::set_coef(T const& value, int pos)
{
    if(pos >=0 && pos <  Manip::nmon(m_oftsh_nvar, m_oftsh_order)) m_oftsh_coef[pos] = value;
    else cout << "Error in set_coef: position is out of scope\n No coefficient is set." << endl;
}

/**
 *  \brief Adds a coefficient at a given position in the polynomial.
 */
template<typename T> void Oftsh<T>::add_coef(T const& value, int pos)
{
    if(pos >=0 && pos <  Manip::nmon(m_oftsh_nvar, m_oftsh_order)) m_oftsh_coef[pos] += value;
    else cout << "Error in set_coef: position is out of scope\n No coefficient is set." << endl;
}

/**
 *  \brief Sets all subcoefficient of the coefficient \c pos to \c value.
 */
template<typename T> template<typename U> void Oftsh<T>::set_sub_coef(U value, int pos)
{
    if(pos >=0 && pos <  Manip::nmon(m_oftsh_nvar, m_oftsh_order)) m_oftsh_coef[pos].set_all_coefs(value);
    else cout << "Error in set_coef: position is out of scope\n No coefficient is set." << endl;
}

/**
 *  \brief Sets subcoefficient \c i of the coefficient \c pos to \c value.
 */
template<typename T> template<typename U> void Oftsh<T>::set_sub_coef(U value, int pos, int i)
{
    if(pos >=0 && pos <  Manip::nmon(m_oftsh_nvar, m_oftsh_order)) m_oftsh_coef[pos].set_coef(value, i);
    else cout << "Error in set_coef: position is out of scope\n No coefficient is set." << endl;
}

/**
 *  \brief Sets random coefficients to all positions in the polynomial.
 */
template<typename T> void Oftsh<T>::set_random_coefs()
{
    for(int pos = 0; pos< Manip::nmon(m_oftsh_nvar, m_oftsh_order); pos++)
    {
            m_oftsh_coef[pos].set_random_coefs();
            //m_oftsh_coef[pos]/=(double)(m_oftsh_order+1);
    }
}

//----------------------------------------------------------------------------------------
//Getters
//----------------------------------------------------------------------------------------
/**
 *  \brief  Gets the first child in the tree.
 */
template<typename T> Oftsh<T> Oftsh<T>::get_term() const
{
    return m_oftsh_term[0];
}

/**
 *  \brief  Gets the child \c i in the tree.
 *
 *  If position is out of scope, the first m_oftsh_term is returned, as in Oftsh<T>::get_term()
 */
template<typename T> Oftsh<T> Oftsh<T>::get_term(int i) const
{
    if(i >=0 && i <= m_oftsh_order) return m_oftsh_term[i];
    else
    {
        cout << "Error in Oftsh<T>::get_term(int i): position is out of scope." << endl;
        cout << "First m_oftsh_term is returned" << endl;
        return m_oftsh_term[0];
    }
}

/**
 *  \brief  Gets the address of the first coefficient
 */
template<typename T> T* Oftsh<T>::get_ptr_first_coef() const
{
    return m_oftsh_coef;
}

/**
 *  \brief  Gets the coefficient at a given position.
 *
 *  If position is out of scope, the first coefficient is returned.
 */
template<typename T>  T Oftsh<T>::get_coef(int i) const
{
    if(i >=0 && i <  Manip::nmon(m_oftsh_nvar, m_oftsh_order)) return m_oftsh_coef[i];
    else
    {
        cout << "Error in Oftsh<T>::get_coef(int i): position is out of scope." << endl;
        cout << "First coefficient is returned" << endl;
        return m_oftsh_coef[0];
    }
}

/**
 *  \brief  Gets the order of the polynomial.
 */
template<typename T> int Oftsh<T>::get_order() const
{
    return m_oftsh_order;
}

/**
 *  \brief  Gets the number of variables of the polynomial.
 */
template<typename T> int Oftsh<T>::get_nvar() const
{
    return m_oftsh_nvar;
}


//----------------------------------------------------------------------------------------
//Zeroing
//----------------------------------------------------------------------------------------
/**
 *  \brief  Sets all coefficients to zero.
 */
template<typename T> void Oftsh<T>::zero()
{
    for(int i = 0 ; i<  Manip::nmon(m_oftsh_nvar, m_oftsh_order); i++) this->set_coef(0.0, i);
}

/**
 *  \brief  Sets all coefficients to zero. Ofsd case
 */
template<> inline void Oftsh< Ofsd >::zero()
{
    for(int i = 0 ; i<  Manip::nmon(m_oftsh_nvar, m_oftsh_order); i++) this->m_oftsh_coef[i].zero();
}

/**
 *  \brief  Sets all coefficients to zero. Ofsc case
 */
template<> inline void Oftsh< Ofsc >::zero()
{
    for(int i = 0 ; i<  Manip::nmon(m_oftsh_nvar, m_oftsh_order); i++) this->m_oftsh_coef[i].zero();
}


//----------------------------------------------------------------------------------------
//Operations
//----------------------------------------------------------------------------------------
//------------------
// Conjugate
//------------------
/**
 *  \brief Conjugates the coefficients the Oftsh object (and only them!). To be used with
 *         evaluate_conjugate to have the true conjugate.
 *   The conjugate of the series \f$ T_n = \sum \limits_{|r| = n} c_r x^r \f$ is
 *  \f$ \bar{T}_n = \sum \sum \limits_{|r| = n} \bar{c}_{r} x^r\f$
 */
template<typename T> Oftsh<T>& Oftsh<T>::conjugate()
{
    for(int pos=0; pos< Manip::nmon(m_oftsh_nvar, m_oftsh_order); pos++) m_oftsh_coef[pos].conjugate();
    return *this;
}

//------------------
// Smult
//------------------
/**
 *  \brief  An operation. Adds the product: \c this \f$  += c a \f$ with c a coefficient.
 *          Warning: set for Ofs coefficients by default.
 */
template<typename T> Oftsh<T>& Oftsh<T>::oftsh_smult_t(Oftsh<T> const& a, T const& c)
{
    if(m_oftsh_order != a.m_oftsh_order || m_oftsh_nvar != a.m_oftsh_nvar)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Oftsh<T> is returned" << endl;
        return *this;
    }
    else
    {
        //Sum
        for (int i=0 ; i <  Manip::nmon(m_oftsh_nvar, m_oftsh_order); i++) m_oftsh_coef[i].ofs_sprod(a.m_oftsh_coef[i], c); //m_oftsh_coef[i] += c*a.m_oftsh_coef[i]; //
        return *this;
    }
}

/**
 *  \brief  An operation. Adds the product: \c this \f$  += c1*c2 a \f$
 *          with c1 and c2 coefficients. Warning: set for Ofs coefficient by default.
 */
template<typename T> Oftsh<T>& Oftsh<T>::oftsh_smult_tt(Oftsh<T> const& a, T const& c1, T const& c2, T& temp)
{
    if(m_oftsh_order != a.m_oftsh_order || m_oftsh_nvar != a.m_oftsh_nvar)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Oftsh<T> is returned" << endl;
        return *this;
    }
    else
    {
        //Sum
        for (int i=0 ; i <  Manip::nmon(m_oftsh_nvar, m_oftsh_order); i++) m_oftsh_coef[i].ofs_smprod_t(a.m_oftsh_coef[i], c1, c2, temp); //m_oftsh_coef[i] += c*a.m_oftsh_coef[i]; //
        return *this;
    }
}

//------------------

/**
 *  \brief  An operation. Adds the product: \c this \f$  += m a \f$ with m a
 *          subcoefficient. Oftsh< Ofs<U> > case.
 */
template<typename T> template<typename U> Oftsh< Ofs<U> >& Oftsh<T>::oftsh_smult_u(Oftsh< Ofs<U> > const& a, U const& c)
{
    if(m_oftsh_order != a.m_oftsh_order || m_oftsh_nvar != a.m_oftsh_nvar)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Oftsh< Ofs<double> > is returned" << endl;
        return *this;
    }
    else
    {
        //Sum
        for (int i=0 ; i <  Manip::nmon(m_oftsh_nvar, m_oftsh_order); i++) m_oftsh_coef[i].ofs_smult(a.m_oftsh_coef[i], c);
        return *this;
    }
}

/**
 *  \brief  An operation. Adds the product: \c this \f$  += m * r_a * a \f$ with m
 *          a subcoefficient and ra a coefficient. Oftsh< Ofs<U> > case.
 */
template<typename T> template<typename U> Oftsh< Ofs<U> >& Oftsh<T>::oftsh_smult_tu(Oftsh< Ofs<U> > const& a, Ofs<U> const& ra, U const& c)
{
    if(m_oftsh_order != a.m_oftsh_order || m_oftsh_nvar != a.m_oftsh_nvar)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Oftsh< Ofs<double> > is returned" << endl;
        return *this;
    }
    else
    {
        //Sum
        for (int i=0 ; i <  Manip::nmon(m_oftsh_nvar, m_oftsh_order); i++) m_oftsh_coef[i].ofs_smprod(a.m_oftsh_coef[i], ra, c);  //sprod(a.m_oftsh_coef[i], c*ra)
        return *this;
    }
}




//------------------
// mult
//------------------
/**
 *  \brief  An operation. Sets the product: \c this \f$  = c a \f$ with c a coefficient.
 *          Warning: set for Ofs coefficient by default.
 */
template<typename T> Oftsh<T>& Oftsh<T>::oftsh_mult_t(Oftsh<T> const& a, T const& c)
{
    if(m_oftsh_order != a.m_oftsh_order || m_oftsh_nvar != a.m_oftsh_nvar)
    {
        cout << "Error using mult: the order and/or number of variables does not match. Initial Oftsh<T> is returned" << endl;
        return *this;
    }
    else
    {
        //Sum
        for (int i=0 ; i <  Manip::nmon(m_oftsh_nvar, m_oftsh_order); i++) m_oftsh_coef[i].ofs_prod(a.m_oftsh_coef[i], c);
        return *this;
    }
}

//------------------

/**
 *  \brief  An operation. Sets the product: \c this \f$  = m a \f$ with m a
 *          subcoefficient. Oftsh< Ofs<U> > case.
 */
template<typename T> template<typename U> Oftsh< Ofs<U> >& Oftsh<T>::oftsh_mult_u(Oftsh< Ofs<U> > const& a, U const& c)
{
    if(m_oftsh_order != a.m_oftsh_order || m_oftsh_nvar != a.m_oftsh_nvar)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Oftsh< Ofs<double> > is returned" << endl;
        return *this;
    }
    else
    {
        //Sum
        for (int i=0 ; i <  Manip::nmon(m_oftsh_nvar, m_oftsh_order); i++) m_oftsh_coef[i].ofs_mult(a.m_oftsh_coef[i], c);
        return *this;
    }
}


//------------------
// sprod
//------------------
/**
 *  \brief  An operation. Adds the product: \c this \f$ += a \times b \f$ with
 *          a and b Oftsh objects
 */
template<typename T> Oftsh<T>& Oftsh<T>::oftsh_sprod(Oftsh<T> const& a, Oftsh<T> const& b)
{
    if (m_oftsh_nvar > 2)  //3+-variate polynomial
    {
        if (a.m_oftsh_order)
        {
            if (b.m_oftsh_order)
            {
                Oftsh<T> *aa , *bb , *pp , *pp0 , *af , *bf ;
                af = a.m_oftsh_term+a.m_oftsh_order;
                bf = b.m_oftsh_term+b.m_oftsh_order;
                for ( aa= a.m_oftsh_term , pp0= this->m_oftsh_term ; aa< af ; aa++ , pp0++)
                {
                    for (bb= b.m_oftsh_term , pp= pp0 ; bb< bf ; bb++ , pp++)  pp->oftsh_sprod(*aa , *bb);
                    pp->oftsh_smult_t(*aa, *(bb->m_oftsh_coef));
                }
                for ( bb= b.m_oftsh_term , pp= pp0 ; bb< bf ; pp->oftsh_smult_t(*bb, *(aa->m_oftsh_coef)), bb++ , pp++);
                *(pp->m_oftsh_coef)+= *(aa->m_oftsh_coef) * *(bb->m_oftsh_coef);
            }
            else this->oftsh_smult_t(a, *(b.m_oftsh_coef)); //b is scalar
        }
        else this->oftsh_smult_t(b, *(a.m_oftsh_coef));   //a is scalar
    }
    else if (m_oftsh_nvar == 2)  //2-variate homogeneous polynomial
    {
        T *aa , *bb , *pp , *pp0 , *af , *bf;
        af = a.m_oftsh_coef +a.m_oftsh_order;
        bf = b.m_oftsh_coef +b.m_oftsh_order;
        for ( aa= a.m_oftsh_coef , pp0= this->m_oftsh_coef ; aa<= af ; aa++ , pp0++)
            for ( bb= b.m_oftsh_coef , pp= pp0 ; bb<= bf ; *pp+= *aa * *bb , bb++ , pp++ ) ;
    }
    else
    {
        *(this->m_oftsh_coef)+= *(a.m_oftsh_coef ) * *(b.m_oftsh_coef); //1-variate homogeneous polynomial
    }
    return *this;
}

//------------------

/**
 *  \brief  An operation. Adds the product: \c this \f$ += a \times b \f$ with
 *          a and b Oftsh objects. Oftsh< Ofsd > double case.
 * Specialization of the routine:
 * template<typename T> Oftsh<T>& Oftsh<T>::oftsh_sprod(Oftsh<T> const& a, Oftsh<T> const& b)
 */
template<> inline Oftsh< Ofsd >& Oftsh< Ofsd >::oftsh_sprod(Oftsh< Ofsd > const& a, Oftsh< Ofsd > const& b)
{
    if (m_oftsh_nvar > 2)  //3+-variate polynomial
    {
        if (a.m_oftsh_order)
        {
            if (b.m_oftsh_order)
            {
                for(int i = 0; i < a.m_oftsh_order ; i++)
                {
                    for(int j = 0; j < b.m_oftsh_order; j++) m_oftsh_term[i+j].oftsh_sprod(a.m_oftsh_term[i], b.m_oftsh_term[j]);
                    m_oftsh_term[i+b.m_oftsh_order].oftsh_smult_t(a.m_oftsh_term[i], b.m_oftsh_term[b.m_oftsh_order].m_oftsh_coef[0]);
                }
                for(int j = 0; j < b.m_oftsh_order; j++) m_oftsh_term[a.m_oftsh_order+j].oftsh_smult_t(b.m_oftsh_term[j], a.m_oftsh_term[a.m_oftsh_order].m_oftsh_coef[0]);
                m_oftsh_term[a.m_oftsh_order+b.m_oftsh_order].m_oftsh_coef[0].ofs_sprod(a.m_oftsh_term[a.m_oftsh_order].m_oftsh_coef[0], b.m_oftsh_term[b.m_oftsh_order].m_oftsh_coef[0]);
            }
            else this->oftsh_smult_t(a, *(b.m_oftsh_coef)); //b is scalar
        }
        else this->oftsh_smult_t(b, *(a.m_oftsh_coef));   //a is scalar
    }
    else if (m_oftsh_nvar == 2)  //2-variate homogeneous polynomial
    {
        for(int i = 0; i<= a.m_oftsh_order; i++) for(int j = 0; j <= b.m_oftsh_order; j++) m_oftsh_coef[i+j].ofs_sprod(a.m_oftsh_coef[i], b.m_oftsh_coef[j]);
    }
    else
    {
        //1-variate homogeneous polynomial
        m_oftsh_coef[0].ofs_sprod(a.m_oftsh_coef[0], b.m_oftsh_coef[0]);
    }
    return *this;
}

/**
 *  \brief  An operation. Adds the product: \c this \f$ += a \times b \f$
 *          with a and b Oftsh objects. Oftsh< Ofsc > cdouble case.
 * Specialization of the routine:
 * template<typename T> Oftsh<T>& Oftsh<T>::oftsh_sprod(Oftsh<T> const& a, Oftsh<T> const& b)
 */
template<> inline Oftsh< Ofsc >& Oftsh< Ofsc >::oftsh_sprod(Oftsh< Ofsc > const& a, Oftsh< Ofsc > const& b)
{
    if (m_oftsh_nvar > 2)  //3+-variate polynomial
    {
        if (a.m_oftsh_order)
        {
            if (b.m_oftsh_order)
            {
                for(int i = 0; i < a.m_oftsh_order ; i++)
                {
                    for(int j = 0; j < b.m_oftsh_order; j++) m_oftsh_term[i+j].oftsh_sprod(a.m_oftsh_term[i], b.m_oftsh_term[j]);
                    m_oftsh_term[i+b.m_oftsh_order].oftsh_smult_t(a.m_oftsh_term[i], b.m_oftsh_term[b.m_oftsh_order].m_oftsh_coef[0]);
                }
                for(int j = 0; j < b.m_oftsh_order; j++) m_oftsh_term[a.m_oftsh_order+j].oftsh_smult_t(b.m_oftsh_term[j], a.m_oftsh_term[a.m_oftsh_order].m_oftsh_coef[0]);
                m_oftsh_term[a.m_oftsh_order+b.m_oftsh_order].m_oftsh_coef[0].ofs_sprod(a.m_oftsh_term[a.m_oftsh_order].m_oftsh_coef[0], b.m_oftsh_term[b.m_oftsh_order].m_oftsh_coef[0]);
            }
            else this->oftsh_smult_t(a, *(b.m_oftsh_coef)); //b is scalar
        }
        else this->oftsh_smult_t(b, *(a.m_oftsh_coef));   //a is scalar
    }
    else if (m_oftsh_nvar == 2)  //2-variate homogeneous polynomial
    {
        for(int i = 0; i<= a.m_oftsh_order; i++) for(int j = 0; j <= b.m_oftsh_order; j++) m_oftsh_coef[i+j].ofs_sprod(a.m_oftsh_coef[i], b.m_oftsh_coef[j]);
    }
    else
    {
        //1-variate homogeneous polynomial
        m_oftsh_coef[0].ofs_sprod(a.m_oftsh_coef[0], b.m_oftsh_coef[0]);
    }
    return *this;
}


//------------------
// smprod
//------------------
/**
 *  \brief  An operation. Adds the product: \c this \f$ += c a \times b \f$
 *          with a and b Oftsh objects and c a coefficient.
 */
template<typename T> Oftsh<T>& Oftsh<T>::oftsh_smprod_t(Oftsh< T > const& a, Oftsh< T> const& b, T const& c, T& temp)
{
    if (m_oftsh_nvar > 2)  //3+-variate polynomial
    {
        if (a.m_oftsh_order)
        {
            if (b.m_oftsh_order)
            {
                Oftsh<T> *aa , *bb , *pp , *pp0 , *af , *bf ;
                af = a.m_oftsh_term+a.m_oftsh_order;
                bf = b.m_oftsh_term+b.m_oftsh_order;
                for ( aa= a.m_oftsh_term , pp0= this->m_oftsh_term ; aa< af ; aa++ , pp0++)
                {
                    for (bb= b.m_oftsh_term , pp= pp0 ; bb< bf ; bb++ , pp++)  pp->oftsh_smprod_t(*aa , *bb, c, temp);
                    pp->oftsh_smult_tt(*aa, *(bb->m_oftsh_coef), c, temp);
                }
                for ( bb= b.m_oftsh_term , pp= pp0 ; bb< bf ; pp->oftsh_smult_tt(*bb, *(aa->m_oftsh_coef), c, temp), bb++ , pp++);
                pp->m_oftsh_coef->ofs_smprod_t(*(aa->m_oftsh_coef), *(bb->m_oftsh_coef), c, temp);
            }
            else this->oftsh_smult_tt(a, *(b.m_oftsh_coef), c, temp); //b is scalar
        }
        else this->oftsh_smult_tt(b, *(a.m_oftsh_coef), c, temp);   //a is scalar
    }
    else if (m_oftsh_nvar == 2)  //2-variate homogeneous polynomial
    {
        T  *aa , *bb , *pp , *pp0 , *af , *bf;
        af = a.m_oftsh_coef +a.m_oftsh_order;
        bf = b.m_oftsh_coef +b.m_oftsh_order;
        for ( aa= a.m_oftsh_coef , pp0= this->m_oftsh_coef ; aa<= af ; aa++ , pp0++)
            for ( bb= b.m_oftsh_coef , pp= pp0 ; bb<= bf ; bb++ , pp++ )  pp->ofs_smprod_t(*aa, *bb, c, temp);//*pp+= *aa * *bb;
    }
    else
    {
        this->m_oftsh_coef->ofs_smprod_t(*(a.m_oftsh_coef ), *(b.m_oftsh_coef), c, temp); //1-variate homogeneous polynomial
    }
    return *this;

//Other implementation less pointers
/*    if (m_oftsh_nvar > 2)  //3+-variate polynomial
    {
        if (a.m_oftsh_order)
        {
            if (b.m_oftsh_order)
            {

                for(int i = 0; i < a.m_oftsh_order ; i++)
                {
                    for(int j = 0; j < b.m_oftsh_order; j++) m_oftsh_term[i+j].oftsh_smprod_t(a.m_oftsh_term[i], b.m_oftsh_term[j], c, temp);
                    m_oftsh_term[i+b.m_oftsh_order].oftsh_smult_tt(a.m_oftsh_term[i], b.m_oftsh_term[b.m_oftsh_order].m_oftsh_coef[0], c, temp);
                }

                for(int j = 0; j < b.m_oftsh_order; j++) m_oftsh_term[a.m_oftsh_order+j].oftsh_smult_tt(b.m_oftsh_term[j], a.m_oftsh_term[a.m_oftsh_order].m_oftsh_coef[0], c, temp);
                m_oftsh_term[a.m_oftsh_order+b.m_oftsh_order].m_oftsh_coef[0].ofs_smprod_t(a.m_oftsh_term[a.m_oftsh_order].m_oftsh_coef[0], b.m_oftsh_term[b.m_oftsh_order].m_oftsh_coef[0], c, temp);
            }
            else this->oftsh_smult_tt(a, *(b.m_oftsh_coef), c, temp); //b is scalar
        }
        else this->oftsh_smult_tt(b, *(a.m_oftsh_coef), c, temp);   //a is scalar
    }
    else if (m_oftsh_nvar == 2)  //2-variate homogeneous polynomial
    {
          for(int i = 0; i<= a.m_oftsh_order; i++) for(int j = 0; j <= b.m_oftsh_order; j++) m_oftsh_coef[i+j].ofs_smprod_t(a.m_oftsh_coef[i], b.m_oftsh_coef[j], c, temp);
    }
    else
    {
        m_oftsh_coef[0].ofs_sprod(a.m_oftsh_coef[0], b.m_oftsh_coef[0]);
    }
    return *this;
*/
}

/**
 *  \brief  An operation. Adds the product: \c this \f$ += m a \times b \f$
 *          with a and b Oftsh objects. Oftsh< Ofs<U> > case.
 */
template<typename T> template<typename U> Oftsh< Ofs<U> >& Oftsh<T>::oftsh_smprod_u(Oftsh< Ofs<U> > const& a, Oftsh< Ofs<U> > const& b, U const& m)
{
    if (m_oftsh_nvar > 2)  //3+-variate polynomial
    {
        if (a.m_oftsh_order)
        {
            if (b.m_oftsh_order)
            {
                Oftsh< Ofs<U> > *aa , *bb , *pp , *pp0 , *af , *bf ;
                af = a.m_oftsh_term+a.m_oftsh_order;
                bf = b.m_oftsh_term+b.m_oftsh_order;
                for ( aa= a.m_oftsh_term , pp0= this->m_oftsh_term ; aa< af ; aa++ , pp0++)
                {
                    for (bb= b.m_oftsh_term , pp= pp0 ; bb< bf ; bb++ , pp++)  pp->oftsh_smprod_u(*aa , *bb, m);
                    pp->oftsh_smult_tu(*aa, *(bb->m_oftsh_coef), m);
                }
                for ( bb= b.m_oftsh_term , pp= pp0 ; bb< bf ; pp->oftsh_smult_tu(*bb, *(aa->m_oftsh_coef), m), bb++ , pp++);
                pp->m_oftsh_coef->ofs_smprod(*(aa->m_oftsh_coef), *(bb->m_oftsh_coef), m);
            }
            else this->oftsh_smult_tu(a, *(b.m_oftsh_coef),m); //b is scalar
        }
        else this->oftsh_smult_tu(b, *(a.m_oftsh_coef), m);   //a is scalar
    }
    else if (m_oftsh_nvar == 2)  //2-variate homogeneous polynomial
    {
        Ofs<U>  *aa , *bb , *pp , *pp0 , *af , *bf;
        af = a.m_oftsh_coef +a.m_oftsh_order;
        bf = b.m_oftsh_coef +b.m_oftsh_order;
        for ( aa= a.m_oftsh_coef , pp0= this->m_oftsh_coef ; aa<= af ; aa++ , pp0++)
            for ( bb= b.m_oftsh_coef , pp= pp0 ; bb<= bf ; bb++ , pp++ )  pp->ofs_smprod(*aa, *bb, m);//*pp+= *aa * *bb;
    }
    else
    {
        this->m_oftsh_coef->ofs_smprod(*(a.m_oftsh_coef ), *(b.m_oftsh_coef), m); //1-variate homogeneous polynomial
    }
    return *this;
}



//----------------------------------------------------------------------------------------
// TFS operations
//----------------------------------------------------------------------------------------
//------------------
// sprod
//------------------
/**
 *  \brief  An operation. Adds the product: \c this \f$ += a \times b \f$
 *          with a and b Oftsh objects
 */
template<typename T> Oftsh<T>& Oftsh<T>::tftsh_sprod(Oftsh<T> const& a, Oftsh<T> const& b)
{
    if (m_oftsh_nvar > 2)  //3+-variate polynomial
    {
        if (a.m_oftsh_order)
        {
            if (b.m_oftsh_order)
            {
                Oftsh<T> *aa , *bb , *pp , *pp0 , *af , *bf ;
                af = a.m_oftsh_term+a.m_oftsh_order;
                bf = b.m_oftsh_term+b.m_oftsh_order;
                for ( aa= a.m_oftsh_term , pp0= this->m_oftsh_term ; aa< af ; aa++ , pp0++)
                {
                    for (bb= b.m_oftsh_term , pp= pp0 ; bb< bf ; bb++ , pp++)  pp->tftsh_sprod(*aa , *bb);
                    pp->tftsh_smult_t(*aa, *(bb->m_oftsh_coef));
                }
                for ( bb= b.m_oftsh_term , pp= pp0 ; bb< bf ; pp->tftsh_smult_t(*bb, *(aa->m_oftsh_coef)), bb++ , pp++);
                pp->m_oftsh_coef->tfs_sprod(*(aa->m_oftsh_coef), *(bb->m_oftsh_coef));
            }
            else this->tftsh_smult_t(a, *(b.m_oftsh_coef)); //b is scalar
        }
        else this->tftsh_smult_t(b, *(a.m_oftsh_coef));   //a is scalar
    }
    else if (m_oftsh_nvar == 2)  //2-variate homogeneous polynomial
    {
        T  *aa , *bb , *pp , *pp0 , *af , *bf;
        af = a.m_oftsh_coef +a.m_oftsh_order;
        bf = b.m_oftsh_coef +b.m_oftsh_order;
        for ( aa= a.m_oftsh_coef , pp0= this->m_oftsh_coef ; aa<= af ; aa++ , pp0++)
            for ( bb= b.m_oftsh_coef , pp= pp0 ; bb<= bf ; bb++ , pp++ )  pp->tfs_sprod(*aa, *bb);
    }
    else
    {
        this->m_oftsh_coef->tfs_sprod(*(a.m_oftsh_coef ), *(b.m_oftsh_coef));//1-variate homogeneous polynomial
    }
    return *this;
}

//------------------
// smprod
//------------------
/**
 *  \brief  An operation. Adds the product: \c this \f$ += c a \times b \f$
 *          with a and b Oftsh objects and c a coefficient.
 */
template<typename T> Oftsh<T>& Oftsh<T>::tftsh_smprod_t(Oftsh< T > const& a, Oftsh< T> const& b, T const& c)
{
    if (m_oftsh_nvar > 2)  //3+-variate polynomial
    {
        if (a.m_oftsh_order)
        {
            if (b.m_oftsh_order)
            {
                Oftsh<T> *aa , *bb , *pp , *pp0 , *af , *bf ;
                af = a.m_oftsh_term+a.m_oftsh_order;
                bf = b.m_oftsh_term+b.m_oftsh_order;
                for ( aa= a.m_oftsh_term , pp0= this->m_oftsh_term ; aa< af ; aa++ , pp0++)
                {
                    for (bb= b.m_oftsh_term , pp= pp0 ; bb< bf ; bb++ , pp++)  pp->tftsh_smprod_t(*aa , *bb, c);
                    pp->tftsh_smult_tt(*aa, *(bb->m_oftsh_coef), c);
                }
                for ( bb= b.m_oftsh_term , pp= pp0 ; bb< bf ; pp->tftsh_smult_tt(*bb, *(aa->m_oftsh_coef), c), bb++ , pp++);
                pp->m_oftsh_coef->tfs_smprod_t(*(aa->m_oftsh_coef), *(bb->m_oftsh_coef), c);
            }
            else this->tftsh_smult_tt(a, *(b.m_oftsh_coef), c); //b is scalar
        }
        else this->tftsh_smult_tt(b, *(a.m_oftsh_coef), c);   //a is scalar
    }
    else if (m_oftsh_nvar == 2)  //2-variate homogeneous polynomial
    {
        T  *aa , *bb , *pp , *pp0 , *af , *bf;
        af = a.m_oftsh_coef +a.m_oftsh_order;
        bf = b.m_oftsh_coef +b.m_oftsh_order;
        for ( aa= a.m_oftsh_coef , pp0= this->m_oftsh_coef ; aa<= af ; aa++ , pp0++)
            for ( bb= b.m_oftsh_coef , pp= pp0 ; bb<= bf ; bb++ , pp++ )  pp->tfs_smprod_t(*aa, *bb, c);//*pp+= *aa * *bb;
    }
    else
    {
        this->m_oftsh_coef->tfs_smprod_t(*(a.m_oftsh_coef ), *(b.m_oftsh_coef), c); //1-variate homogeneous polynomial
    }
    return *this;
}

/**
 *  \brief  An operation. Adds the product: \c this \f$ += c a \times b \f$
 *          with a and b Oftsh objects and c a coefficient.
 */
template<typename T> template<typename U> Oftsh<T>& Oftsh<T>::tftsh_smprod_tu(Oftsh< T > const& a, Oftsh< T> const& b, T const& c, U const& m)
{
    if (m_oftsh_nvar > 2)  //3+-variate polynomial
    {
        if (a.m_oftsh_order)
        {
            if (b.m_oftsh_order)
            {
                Oftsh<T> *aa , *bb , *pp , *pp0 , *af , *bf ;
                af = a.m_oftsh_term+a.m_oftsh_order;
                bf = b.m_oftsh_term+b.m_oftsh_order;
                for ( aa= a.m_oftsh_term , pp0= this->m_oftsh_term ; aa< af ; aa++ , pp0++)
                {
                    for (bb= b.m_oftsh_term , pp= pp0 ; bb< bf ; bb++ , pp++)  pp->tftsh_smprod_tu(*aa , *bb, c, m);
                    pp->tftsh_smult_ttu(*aa, *(bb->m_oftsh_coef), c, m);
                }
                for ( bb= b.m_oftsh_term , pp= pp0 ; bb< bf ; pp->tftsh_smult_ttu(*bb, *(aa->m_oftsh_coef), c, m), bb++ , pp++);
                pp->m_oftsh_coef->tfs_smprod_tu(*(aa->m_oftsh_coef), *(bb->m_oftsh_coef), c, m);
            }
            else this->tftsh_smult_ttu(a, *(b.m_oftsh_coef), c, m); //b is scalar
        }
        else this->tftsh_smult_ttu(b, *(a.m_oftsh_coef), c, m);   //a is scalar
    }
    else if (m_oftsh_nvar == 2)  //2-variate homogeneous polynomial
    {
        T  *aa , *bb , *pp , *pp0 , *af , *bf;
        af = a.m_oftsh_coef +a.m_oftsh_order;
        bf = b.m_oftsh_coef +b.m_oftsh_order;
        for ( aa= a.m_oftsh_coef , pp0= this->m_oftsh_coef ; aa<= af ; aa++ , pp0++)
            for ( bb= b.m_oftsh_coef , pp= pp0 ; bb<= bf ; bb++ , pp++ )  pp->tfs_smprod_tu(*aa, *bb, c, m);//*pp+= *aa * *bb;
    }
    else
    {
        this->m_oftsh_coef->tfs_smprod_tu(*(a.m_oftsh_coef ), *(b.m_oftsh_coef), c, m); //1-variate homogeneous polynomial
    }
    return *this;
}

/**
 *  \brief  An operation. Adds the product: \c this \f$ += m a \times b \f$ with a and b
 *          Oftsh objects. Oftsh< Ofs<U> > case.
 */
template<typename T> template<typename U> Oftsh< Ofs<U> >& Oftsh<T>::tftsh_smprod_u(Oftsh< Ofs<U> > const& a, Oftsh< Ofs<U> > const& b, U const& m)
{
    if (m_oftsh_nvar > 2)  //3+-variate polynomial
    {
        if (a.m_oftsh_order)
        {
            if (b.m_oftsh_order)
            {
                Oftsh< Ofs<U> > *aa , *bb , *pp , *pp0 , *af , *bf ;
                af = a.m_oftsh_term+a.m_oftsh_order;
                bf = b.m_oftsh_term+b.m_oftsh_order;
                for ( aa= a.m_oftsh_term , pp0= this->m_oftsh_term ; aa< af ; aa++ , pp0++)
                {
                    for (bb= b.m_oftsh_term , pp= pp0 ; bb< bf ; bb++ , pp++)  pp->tftsh_smprod_u(*aa , *bb, m);
                    pp->tftsh_smult_tu(*aa, *(bb->m_oftsh_coef), m);
                }
                for ( bb= b.m_oftsh_term , pp= pp0 ; bb< bf ; pp->tftsh_smult_tu(*bb, *(aa->m_oftsh_coef), m), bb++ , pp++);
                pp->m_oftsh_coef->tfs_smprod(*(aa->m_oftsh_coef), *(bb->m_oftsh_coef), m);
            }
            else this->tftsh_smult_tu(a, *(b.m_oftsh_coef),m); //b is scalar
        }
        else this->tftsh_smult_tu(b, *(a.m_oftsh_coef), m);   //a is scalar
    }
    else if (m_oftsh_nvar == 2)  //2-variate homogeneous polynomial
    {
        Ofs<U>  *aa , *bb , *pp , *pp0 , *af , *bf;
        af = a.m_oftsh_coef +a.m_oftsh_order;
        bf = b.m_oftsh_coef +b.m_oftsh_order;
        for ( aa= a.m_oftsh_coef , pp0= this->m_oftsh_coef ; aa<= af ; aa++ , pp0++)
            for ( bb= b.m_oftsh_coef , pp= pp0 ; bb<= bf ; bb++ , pp++ )  pp->tfs_smprod(*aa, *bb, m);//*pp+= *aa * *bb;
    }
    else
    {
        this->m_oftsh_coef->tfs_smprod(*(a.m_oftsh_coef ), *(b.m_oftsh_coef), m); //1-variate homogeneous polynomial
    }
    return *this;
}

//------------------
// smult
//------------------
/**
 *  \brief  An operation. Adds the product: \c this \f$  += m a \f$
 *          with m a subcoefficient. Oftsh< Ofs<U> > case.
 */
template<typename T> template<typename U> Oftsh< Ofs<U> >& Oftsh<T>::tftsh_smult_u(Oftsh< Ofs<U> > const& a, U const& c)
{
    if(m_oftsh_order != a.m_oftsh_order || m_oftsh_nvar != a.m_oftsh_nvar)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Oftsh< Ofs<double> > is returned" << endl;
        return *this;
    }
    else
    {
        //Sum: note ofs_smult can be used because it does the same as tfs_smult would do.
        for (int i=0 ; i <  Manip::nmon(m_oftsh_nvar, m_oftsh_order); i++) m_oftsh_coef[i].ofs_smult(a.m_oftsh_coef[i], c);
        return *this;
    }
}

/**
 *  \brief  An operation. Sets the product: \c this \f$  = m a \f$
 *          m a subcoefficient. Oftsh< Ofs<U> > case.
 */
template<typename T> template<typename U> Oftsh< Ofs<U> >& Oftsh<T>::tftsh_mult_u(Oftsh< Ofs<U> > const& a, U const& c)
{
    if(m_oftsh_order != a.m_oftsh_order || m_oftsh_nvar != a.m_oftsh_nvar)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Oftsh< Ofs<double> > is returned" << endl;
        return *this;
    }
    else
    {
        //Sum: note ofs_mult can be used because it does the same as tfs_smult would do.
        for (int i=0 ; i <  Manip::nmon(m_oftsh_nvar, m_oftsh_order); i++) m_oftsh_coef[i].ofs_mult(a.m_oftsh_coef[i], c);
        return *this;
    }
}

/**
 *  \brief  An operation. Adds the product: \c this \f$  += c a \f$
 *          with c a coefficient. Warning: set for Ofs coefficient by default.
 */
template<typename T> Oftsh<T>& Oftsh<T>::tftsh_smult_t(Oftsh<T> const& a, T const& c)
{
    if(m_oftsh_order != a.m_oftsh_order || m_oftsh_nvar != a.m_oftsh_nvar)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Oftsh<T> is returned" << endl;
        return *this;
    }
    else
    {
        //Sum
        for (int i=0 ; i <  Manip::nmon(m_oftsh_nvar, m_oftsh_order); i++) m_oftsh_coef[i].tfs_sprod(a.m_oftsh_coef[i], c); //m_oftsh_coef[i] += c*a.m_oftsh_coef[i]; //
        return *this;
    }
}

/**
 *  \brief  An operation. Adds the product: \c this \f$  += c1*c2 a \f$
 *          with c1 and c2 coefficients. Warning: set for Ofs coefficient by default.
 */
template<typename T> Oftsh<T>& Oftsh<T>::tftsh_smult_tt(Oftsh<T> const& a, T const& c1, T const& c2)
{
    if(m_oftsh_order != a.m_oftsh_order || m_oftsh_nvar != a.m_oftsh_nvar)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Oftsh<T> is returned" << endl;
        return *this;
    }
    else
    {
        //Sum
        for (int i=0 ; i <  Manip::nmon(m_oftsh_nvar, m_oftsh_order); i++) m_oftsh_coef[i].tfs_smprod_t(a.m_oftsh_coef[i], c1, c2); //m_oftsh_coef[i] += c*a.m_oftsh_coef[i]; //
        return *this;
    }
}

/**
 *  \brief  An operation. Adds the product: \c this \f$  += c1*c2 a \f$
 *          with c1 and c2 coefficients. Warning: set for Ofs coefficient by default.
 */
template<typename T> template<typename U> Oftsh<T>& Oftsh<T>::tftsh_smult_ttu(Oftsh<T> const& a, T const& c1, T const& c2, U const& m)
{
    if(m_oftsh_order != a.m_oftsh_order || m_oftsh_nvar != a.m_oftsh_nvar)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Oftsh<T> is returned" << endl;
        return *this;
    }
    else
    {
        //Sum
        for (int i=0 ; i <  Manip::nmon(m_oftsh_nvar, m_oftsh_order); i++) m_oftsh_coef[i].tfs_smprod_tu(a.m_oftsh_coef[i], c1, c2, m); //m_oftsh_coef[i] += c*a.m_oftsh_coef[i]; //
        return *this;
    }
}

/**
 *  \brief  An operation. Adds the product: \c this \f$  += m * r_a * a \f$
 *          with m a subcoefficient and ra a coefficient. Oftsh< Ofs<U> > case.
 */
template<typename T> template<typename U> Oftsh< Ofs<U> >& Oftsh<T>::tftsh_smult_tu(Oftsh< Ofs<U> > const& a, Ofs<U> const& ra, U const& c)
{
    if(m_oftsh_order != a.m_oftsh_order || m_oftsh_nvar != a.m_oftsh_nvar)
    {
        cout << "Error using smult: the order and/or number of variables does not match. Initial Oftsh< Ofs<double> > is returned" << endl;
        return *this;
    }
    else
    {
        //Sum
        for (int i=0 ; i <  Manip::nmon(m_oftsh_nvar, m_oftsh_order); i++) m_oftsh_coef[i].tfs_smprod(a.m_oftsh_coef[i], ra, c);  //sprod(a.m_oftsh_coef[i], c*ra)
        return *this;
    }
}

//------------------
// derh
//------------------
/**
 *  \brief  An operation. Applies the partial derivative with respect
 *          to the variable \c ni: this \f$ = \frac{\partial a}{\partial x_{ni}} \f$
 *
 *  Notes:
 *  1. If a is of order n, this is of order n-1.
 *  2. Need an extension to autonomous case: WORKS ONLY FOR Ofs<double/cdouble> coefficients.
 */
template<typename T> Oftsh<T>& Oftsh<T>::tfts_derh(Oftsh<T> const& a, int ni)
{
    Oftsh<T> *dd, *pp, *pf;
    pf = a.m_oftsh_term+a.m_oftsh_order;
    if(this->m_oftsh_nvar > ni)  //the number of variables is greater than ni
    {

        for(pp = a.m_oftsh_term, dd=this->m_oftsh_term; pp < pf; pp++,dd++)
        {
            dd->tfts_derh(*pp, ni);
        }
    }
    else
    {
        if(ni==1)
        {
            this->set_coef( ((double)a.m_oftsh_order+0.0*I)*a.get_coef(0), 0);
        }
        else
        {
            for(dd=this->m_oftsh_term, pp=a.m_oftsh_term+1; pp<=pf; pp++,dd++)
            {
                dd->tftsh_mult_u(*pp, (double)(pp - a.m_oftsh_term)+0.0*I);
            }
        }
    }
    return *this;
}

/**
 *  \brief  An operation. Adds the partial derivative with respect to
 *          the variable \c ni: this \f$ += \frac{\partial a}{\partial x_{ni}} \f$
 *
 *  Notes:
 *  1. If a is of order n, this is of order n-1.
 *  2. Need an extension to autonomous case: WORKS ONLY FOR Ofs<double/cdouble> coefficients.
 */
template<typename T> Oftsh<T>& Oftsh<T>::tfts_sderh(Oftsh< T > const& a, int ni)
{
    Oftsh<T> *dd, *pp, *pf;
    pf = a.m_oftsh_term+a.m_oftsh_order;
    if(this->m_oftsh_nvar > ni)  //the number of variables is greater than ni
    {

        for(pp = a.m_oftsh_term, dd=this->m_oftsh_term; pp < pf; pp++,dd++)
        {
            dd->tfts_sderh(*pp, ni);
        }
    }
    else
    {
        if(ni==1)
        {
            this->add_coef( (a.m_oftsh_order+0.0*I)*a.get_coef(0), 0);
        }
        else
        {
            for(dd=this->m_oftsh_term, pp=a.m_oftsh_term+1; pp<=pf; pp++,dd++)
            {
                dd->tftsh_smult_u(*pp, (pp - a.m_oftsh_term)+0.0*I);
            }
        }
    }

    return *this;
}


//----------------------------------------------------------------------------------------
// Derivation
//----------------------------------------------------------------------------------------
/**
 *  \brief  An operation. Applies the partial derivative with respect
 *          to the variable \c ni: this \f$ = \frac{\partial a}{\partial x_{ni}} \f$
 *
 *  Notes:
 *  1. If a is of order n, this is of order n-1.
 *  2. Need an extension to autonomous case: WORKS ONLY FOR Ofs<double/cdouble> coefficients.
 */
template<typename T> Oftsh<T>& Oftsh<T>::derh(Oftsh<T> const& a, int ni)
{
    Oftsh<T> *dd, *pp, *pf;
    pf = a.m_oftsh_term+a.m_oftsh_order;
    if(this->m_oftsh_nvar > ni)  //the number of variables is greater than ni
    {

        for(pp = a.m_oftsh_term, dd=this->m_oftsh_term; pp < pf; pp++,dd++)
        {
            dd->derh(*pp, ni);
        }
    }
    else
    {
        if(ni==1)
        {
            this->set_coef( ((double)a.m_oftsh_order+0.0*I)*a.get_coef(0), 0);
        }
        else
        {
            for(dd=this->m_oftsh_term, pp=a.m_oftsh_term+1; pp<=pf; pp++,dd++)
            {
                dd->oftsh_mult_u(*pp, (double)(pp - a.m_oftsh_term)+0.0*I);
            }
        }
    }
    return *this;
}

/**
 *  \brief  An operation. Applies the partial derivative with respect
 *          to the variable \c ni: this \f$ = \frac{\partial a}{\partial x_{ni}} \f$
 *   Specialization of  Oftsh<T>::derh(Oftsh< T > const& a, int ni) to Ofsd coefficients.
 */
template<> inline Oftsh<Ofsd >& Oftsh<Ofsd >::derh(Oftsh<Ofsd  > const& a, int ni)
{
    Oftsh<Ofsd > *dd, *pp, *pf;
    pf = a.m_oftsh_term+a.m_oftsh_order;
    if(this->m_oftsh_nvar > ni)  //the number of variables is greater than ni
    {

        for(pp = a.m_oftsh_term, dd=this->m_oftsh_term; pp < pf; pp++,dd++)
        {
            dd->derh(*pp, ni);
        }
    }
    else
    {
        if(ni==1)
        {
            this->set_coef( ((double)a.m_oftsh_order)*a.get_coef(0), 0);
        }
        else
        {
            for(dd=this->m_oftsh_term, pp=a.m_oftsh_term+1; pp<=pf; pp++,dd++)
            {
                dd->oftsh_mult_u(*pp, (double) (pp - a.m_oftsh_term));
            }
        }
    }

    return *this;
}

/**
 *  \brief  An operation. Adds the partial derivative with respect
 *          to the variable \c ni: this \f$ += \frac{\partial a}{\partial x_{ni}} \f$
 *
 *  Notes:
 *  1. If a is of order n, this is of order n-1.
 *  2. Need an extension to autonomous case: WORKS ONLY FOR Ofs<double/cdouble> coefficients.
 */
template<typename T> Oftsh<T>& Oftsh<T>::sderh(Oftsh< T > const& a, int ni)
{
    Oftsh<T> *dd, *pp, *pf;
    pf = a.m_oftsh_term+a.m_oftsh_order;
    if(this->m_oftsh_nvar > ni)  //the number of variables is greater than ni
    {

        for(pp = a.m_oftsh_term, dd=this->m_oftsh_term; pp < pf; pp++,dd++)
        {
            dd->sderh(*pp, ni);
        }
    }
    else
    {
        if(ni==1)
        {
            this->add_coef( (a.m_oftsh_order+0.0*I)*a.get_coef(0), 0);
        }
        else
        {
            for(dd=this->m_oftsh_term, pp=a.m_oftsh_term+1; pp<=pf; pp++,dd++)
            {
                dd->oftsh_smult_u(*pp, (pp - a.m_oftsh_term)+0.0*I);
            }
        }
    }

    return *this;
}

/**
 *  \brief  An operation. Adds the partial primitive with respect to the variable \c ni
 *
 *  Notes:
 *  1. If a is of order n, this is of order n+1.
 *  2. Need an extension to autonomous case: WORKS ONLY FOR Ofs<double/cdouble> coefficients.
 */
template<typename T> Oftsh<T>& Oftsh<T>::sprimh(Oftsh< T > const& a, int ni)
{
    Oftsh<T> *dd, *pp, *pf;
    pf = a.m_oftsh_term+a.m_oftsh_order;
    if(this->m_oftsh_nvar > ni)  //the number of variables is greater than ni
    {
        for(pp = a.m_oftsh_term, dd=this->m_oftsh_term; pp <= pf; pp++,dd++)
        {
            dd->sprimh(*pp, ni);
        }
    }
    else
    {
        if(ni==1)
        {
            this->add_coef( 1.0/(a.m_oftsh_order+1.0+0.0*I)*a.get_coef(0), 0);
        }
        else
        {
            for(dd=this->m_oftsh_term+1, pp=a.m_oftsh_term; pp<=pf; pp++,dd++)
            {
                dd->oftsh_smult_u(*pp, 1.0/(pp - a.m_oftsh_term+0.0*I+1.0));
            }
        }
    }

    return *this;
}

/**
 *  \brief  An operation. Adds the partial derivative with respect
 *          to the variable \c ni: this \f$ += \frac{\partial a}{\partial x_{ni}} \f$
 *   Specialization of  Oftsh<T>::sderh(Oftsh< T > const& a, int ni) to Ofsd coefficients.
 */
template<> inline Oftsh<Ofsd >& Oftsh<Ofsd >::sderh(Oftsh< Ofsd > const& a, int ni)
{
    Oftsh<Ofsd > *dd, *pp, *pf;
    pf = a.m_oftsh_term+a.m_oftsh_order;
    if(this->m_oftsh_nvar > ni)  //the number of variables is greater than ni
    {

        for(pp = a.m_oftsh_term, dd=this->m_oftsh_term; pp < pf; pp++,dd++)
        {
            dd->sderh(*pp, ni);
        }
    }
    else
    {
        if(ni==1)
        {
            this->add_coef( ((double)a.m_oftsh_order)*a.get_coef(0), 0);
        }
        else
        {
            for(dd=this->m_oftsh_term, pp=a.m_oftsh_term+1; pp<=pf; pp++,dd++)
            {
                dd->oftsh_smult_u(*pp, (double) (pp - a.m_oftsh_term));
            }
        }
    }

    return *this;
}

/**
 *  \brief  An operation. Set the time derivative of object \c a
 *          with pulsation \f$ \omega = n \f$, so that \c this \f$ = \dot{a} \f$.
 */
template<typename T> Oftsh<T>& Oftsh<T>::dot(Oftsh<T> const& a, double const&  n)
{
    if(m_oftsh_order != a.m_oftsh_order || m_oftsh_nvar != a.m_oftsh_nvar)
    {
        cout << "Error using dot: the order and/or number of variables does not match. Initial Oftsh<T> is returned" << endl;
        return *this;
    }
    else
    {
        //Derivation of all the coefficients
        for (int i=0 ; i <  Manip::nmon(m_oftsh_nvar, m_oftsh_order); i++) m_oftsh_coef[i].dot(a.m_oftsh_coef[i], n);
        return *this;
    }
}

//----------------------------------------------------------------------------------------
// Evaluate
//----------------------------------------------------------------------------------------
/**
 *  \brief  Evaluates the Ofs object at coordinates X and set it in \c z: \c z \f$ = T_n(X) \f$. cdouble version
 */
template<typename T> template<typename U> void Oftsh<T>::evaluate(U X[], T& z)
{
    //Zeroing z
    z.zero();
    int *kv = (int*) calloc(m_oftsh_order, sizeof(int));
    kv[0] = m_oftsh_order;
    for(int i=1; i<m_oftsh_nvar; i++) kv[i] = 0;
    U aux;
    U bux;

    for (int i=0; i< Manip::nmon(m_oftsh_nvar, m_oftsh_order); i++)
    {
        //z += X[ii]^kv[ii]*m_oftsh_coef(i)
        bux = 1.0+0.0*I;
        for(int ii = 0; ii < m_oftsh_nvar; ii++)
        {
            aux = 1.0+0.0*I;
            if(kv[ii] != 0.0)
            {
                for(int j = 1; j <= kv[ii]; j++) aux*= X[ii];
            }
            bux *= aux;
        }

        z.ofs_smult(m_oftsh_coef[i], (U) bux);
        if(i< Manip::nmon(m_oftsh_nvar, m_oftsh_order)-1)  Manip::prxkt(kv, m_oftsh_nvar);
    }
    free(kv);
}

/**
 *  \brief  Evaluates the Ofs object at coordinates X and set it in \c z: \c z \f$ = T_n(X) \f$. double version.
 */
template<typename T>  void Oftsh<T>::evaluate(double X[], T& z)
{
    //Zeroing z
    z.zero();
    int *kv = (int*) calloc(m_oftsh_order, sizeof(int));
    kv[0] = m_oftsh_order;
    for(int i=1; i<m_oftsh_nvar; i++) kv[i] = 0;
    double aux;
    double bux;

    for (int i=0; i< Manip::nmon(m_oftsh_nvar, m_oftsh_order); i++)
    {
        //z += X[ii]^kv[ii]*m_oftsh_coef(i)
        bux = 1.0;
        for(int ii = 0; ii < m_oftsh_nvar; ii++)
        {
            aux = 1.0;
            if(kv[ii] != 0.0)
            {
                for(int j = 1; j <= kv[ii]; j++) aux*= X[ii];
            }
            bux *= aux;
        }

        z.ofs_smult(m_oftsh_coef[i], (double) bux);
        if(i< Manip::nmon(m_oftsh_nvar, m_oftsh_order)-1)  Manip::prxkt(kv, m_oftsh_nvar);
    }
    free(kv);
}


/**
 *  \brief  Evaluates the Ofs object at coordinates X and adds it in \c z: \c z \f$ += T_n(X) \f$. cdouble version
 */
template<typename T> template<typename U> void Oftsh<T>::sevaluate(U X[], T& z)
{
    int *kv = (int*) calloc(m_oftsh_order, sizeof(int));
    kv[0] = m_oftsh_order;
    for(int i=1; i<m_oftsh_nvar; i++) kv[i] = 0;
    U aux, bux;

    for (int i=0; i< Manip::nmon(m_oftsh_nvar, m_oftsh_order); i++)
    {
        //z += X[ii]^kv[ii]*m_oftsh_coef(i)
        bux = (U) (1.0+0.0*I);
        for(int ii = 0; ii < m_oftsh_nvar; ii++)
        {
            aux = (U) (1.0+0.0*I);
            if(kv[ii] != 0.0) for(int j = 1; j <= kv[ii]; j++) aux*= X[ii];
            bux *= aux;
        }
        z.ofs_smult(m_oftsh_coef[i], (U) bux);
        if(i< Manip::nmon(m_oftsh_nvar, m_oftsh_order)-1)  Manip::prxkt(kv, m_oftsh_nvar);
    }
    free(kv);
}

/**
 *  \brief  Evaluates the Ofs object at coordinates X and adds it in \c z: \c z \f$ += T_n(X) \f$ at order ofs_order. cdouble version
 */
template<typename T> template<typename U> void Oftsh<T>::sevaluate(U X[], T& z, int const& ofs_order)
{
    int *kv = (int*) calloc(m_oftsh_nvar, sizeof(int));
    kv[0] = m_oftsh_order;
    for(int i=1; i<m_oftsh_nvar; i++) kv[i] = 0;
    U aux, bux;

    for (int i=0; i< Manip::nmon(m_oftsh_nvar, m_oftsh_order); i++)
    {
        //Evaluate one coefficient
        if(!m_oftsh_coef[i].is_null(1))
        {
            //z += X[ii]^kv[ii]*m_oftsh_coef(i)
            bux = (U) (1.0+0.0*I);
            for(int ii = 0; ii < m_oftsh_nvar; ii++)
            {
                aux = (U) (1.0+0.0*I);
                if(kv[ii] != 0.0) for(int j = 1; j <= kv[ii]; j++) aux*= X[ii];
                bux *= aux;
            }
            z.ofs_smult(m_oftsh_coef[i], (U) bux, ofs_order);
       }
        if(i< Manip::nmon(m_oftsh_nvar, m_oftsh_order)-1)  Manip::prxkt(kv, m_oftsh_nvar); //update the exponents
    }
    free(kv);
}


/**
 *  \brief  Evaluates the Ofs object at coordinates X and time t and adds it in \c z: \c z \f$ += T_n(X) \f$. cdouble version
 */
template<typename T> template<typename U> void Oftsh<T>::sevaluatedc(U X[], U& z, double const& t, int const& ofs_order)
{
    int *kv = (int*) calloc(m_oftsh_nvar, sizeof(int));
    kv[0] = m_oftsh_order;
    for(int i=1; i<m_oftsh_nvar; i++) kv[i] = 0;
    U aux, bux;

    for (int i=0; i< Manip::nmon(m_oftsh_nvar, m_oftsh_order); i++)
    {
        //z += X[ii]^kv[ii]*m_oftsh_coef(i)
        bux = (U) (1.0+0.0*I);
        for(int ii = 0; ii < m_oftsh_nvar; ii++)
        {
            aux = (U) (1.0+0.0*I);
            if(kv[ii] != 0.0) for(int j = 1; j <= kv[ii]; j++) aux*= X[ii];
            bux *= aux;
        }
        z += m_oftsh_coef[i].evaluate(t, ofs_order)*bux;
        if(i< Manip::nmon(m_oftsh_nvar, m_oftsh_order)-1)  Manip::prxkt(kv, m_oftsh_nvar);
    }
    free(kv);
}


/**
 *  \brief  Evaluates the Ofs object at coordinates X and time t and adds it in \c z: \c z \f$ += T_n(X) \f$. cdouble version
 */
template<typename T> template<typename U> void Oftsh<T>::fevaluate(U X[], U& z, int kv[], double cR[], double sR[], int const& ofs_order)
{
    //cout << "oftsh_fevaluate. m_oftsh_order = " << m_oftsh_order << endl;
    U aux, bux;
    //Initialize kv
    kv[0] = m_oftsh_order;
    for(int i=1; i<m_oftsh_nvar; i++) kv[i] = 0;
    //Loop on all monomials
    for (int i=0; i< Manip::nmon(m_oftsh_nvar, m_oftsh_order); i++)
    {
        //z += X[ii]^kv[ii]*m_oftsh_coef(i)
        bux = (U) (1.0+0.0*I);
        for(int ii = 0; ii < m_oftsh_nvar; ii++)
        {
            aux = (U) (1.0+0.0*I);
            if(kv[ii] != 0.0) for(int j = 1; j <= kv[ii]; j++) aux*= X[ii];
            bux *= aux;
        }
        z += m_oftsh_coef[i].fevaluate(cR, sR, ofs_order)*bux;
        if(i< Manip::nmon(m_oftsh_nvar, m_oftsh_order)-1)  Manip::prxkt(kv, m_oftsh_nvar);
    }
}

/**
 *  \brief  Evaluates the Ofs object at coordinates X and adds it in \c z: \c z \f$ += T_n(X) \f$. double version
 */
template<typename T> void Oftsh<T>::sevaluate(double X[], T& z)
{
    int *kv = (int*) calloc(m_oftsh_order, sizeof(int));
    kv[0] = m_oftsh_order;
    for(int i=1; i<m_oftsh_nvar; i++) kv[i] = 0;
    double aux, bux;

    for (int i=0; i< Manip::nmon(m_oftsh_nvar, m_oftsh_order); i++)
    {
        //z += X[ii]^kv[ii]*m_oftsh_coef(i)
        bux = 1.0;
        for(int ii = 0; ii < m_oftsh_nvar; ii++)
        {
            aux = 1.0;
            if(kv[ii] != 0.0) for(int j = 1; j <= kv[ii]; j++) aux*= X[ii];
            bux *= aux;
        }

        z.ofs_smult(m_oftsh_coef[i], bux);
        if(i< Manip::nmon(m_oftsh_nvar, m_oftsh_order)-1)  Manip::prxkt(kv, m_oftsh_nvar);
    }
    free(kv);
}


/**
 *  \brief  Evaluates the Ofs object at the conjugate of coordinates X and set it in \c z: \c z \f$ = T_n(\bar{X}) \f$.
 */
template<typename T> template<typename U> void Oftsh<T>::evaluate_conjugate(U X[], T& z)
{
    //Zeroing z
    z.zero();
    int kv[m_oftsh_nvar];
    kv[0] = m_oftsh_order;
    for(int i=1; i<m_oftsh_nvar; i++) kv[i] = 0;
    U aux;
    U bux;

    for (int i=0; i< Manip::nmon(m_oftsh_nvar, m_oftsh_order); i++)
    {
        //z += X[ii]^kv[ii]*m_oftsh_coef(i)
        bux = 1.0+0.0*I;
        for(int ii = 0; ii < m_oftsh_nvar; ii++)
        {
            aux = 1.0+0.0*I;
            if(kv[ii] != 0.0)
            {
                for(int j = 1; j <= kv[ii]; j++) aux*= conj(X[ii]);
            }
            bux *= aux;
        }

        z.ofs_smult(m_oftsh_coef[i], (U) bux);
        if(i< Manip::nmon(m_oftsh_nvar, m_oftsh_order)-1)  Manip::prxkt(kv, m_oftsh_nvar);
    }
}

/**
 *  \brief  Evaluates the Ofs object at the conjugate of coordinates X and set it in \c z: \c z \f$ = T_n(\bar{X}) \f$.
 */
template<typename T> void Oftsh<T>::evaluate_conjugate(double X[], T& z)
{
    //Zeroing z
    z.zero();
    int kv[m_oftsh_nvar];
    kv[0] = m_oftsh_order;
    for(int i=1; i<m_oftsh_nvar; i++) kv[i] = 0;
    double aux;
    double bux;

    for (int i=0; i< Manip::nmon(m_oftsh_nvar, m_oftsh_order); i++)
    {
        //z += X[ii]^kv[ii]*m_oftsh_coef(i)
        bux = 1.0;
        for(int ii = 0; ii < m_oftsh_nvar; ii++)
        {
            aux = 1.0;
            if(kv[ii] != 0.0)
            {
                for(int j = 1; j <= kv[ii]; j++) aux*= X[ii];
            }
            bux *= aux;
        }

        z.ofs_smult(m_oftsh_coef[i], (double) bux);
        if(i< Manip::nmon(m_oftsh_nvar, m_oftsh_order)-1)  Manip::prxkt(kv, m_oftsh_nvar);
    }
}


/**
 *  \brief  Evaluates the Ofs object at the conjugate of coordinates X and set it in \c z: \c z \f$ += T_n(\bar{X}) \f$.
 */
template<typename T> template<typename U> void Oftsh<T>::sevaluate_conjugate(U X[], T& z)
{
    int kv[m_oftsh_nvar];
    kv[0] = m_oftsh_order;
    for(int i=1; i<m_oftsh_nvar; i++) kv[i] = 0;
    U aux;
    U bux;

    for (int i=0; i< Manip::nmon(m_oftsh_nvar, m_oftsh_order); i++)
    {
        //z += X[ii]^kv[ii]*m_oftsh_coef(i)
        bux = 1.0+0.0*I;
        for(int ii = 0; ii < m_oftsh_nvar; ii++)
        {
            aux = 1.0+0.0*I;
            if(kv[ii] != 0.0)
            {
                for(int j = 1; j <= kv[ii]; j++) aux*= conj(X[ii]);
            }
            bux *= aux;
        }

        z.ofs_smult(m_oftsh_coef[i], (U) bux);
        if(i< Manip::nmon(m_oftsh_nvar, m_oftsh_order)-1)  Manip::prxkt(kv, m_oftsh_nvar);
    }
}

/**
 *  \brief  Evaluates the Ofs object at the conjugate of coordinates X and set it in \c z: \c z \f$ += T_n(\bar{X}) \f$.
 */
template<typename T> void Oftsh<T>::sevaluate_conjugate(double X[], T& z)
{
    int kv[m_oftsh_nvar];
    kv[0] = m_oftsh_order;
    for(int i=1; i<m_oftsh_nvar; i++) kv[i] = 0;
    double aux;
    double bux;

    for (int i=0; i< Manip::nmon(m_oftsh_nvar, m_oftsh_order); i++)
    {
        //z += X[ii]^kv[ii]*m_oftsh_coef(i)
        bux = 1.0;
        for(int ii = 0; ii < m_oftsh_nvar; ii++)
        {
            aux = 1.0;
            if(kv[ii] != 0.0)
            {
                for(int j = 1; j <= kv[ii]; j++) aux*= X[ii];
            }
            bux *= aux;
        }

        z.ofs_smult(m_oftsh_coef[i], (double) bux);
        if(i< Manip::nmon(m_oftsh_nvar, m_oftsh_order)-1)  Manip::prxkt(kv, m_oftsh_nvar);
    }
}


/**
 *  \brief  Evaluates the \f$ L_1 \f$ norm of the current Oftsh object.
 */
template<typename T> double Oftsh<T>::l1norm()
{
    double l1n = 0.0;
    if(MODEL_TYPE == Csts::QBCP) for(int k =0; k < Manip::nmon(m_oftsh_nvar, m_oftsh_order); k++) l1n += cabs(this->get_coef(k).l1norm());
    else if(MODEL_TYPE == Csts::CRTBP) for(int k =0; k < Manip::nmon(m_oftsh_nvar, m_oftsh_order); k++) l1n += cabs(this->get_coef(k).ofs_get_coef(0));
    return l1n;
}



/**
 *  \brief  Number of small divisors under a certain value
 */
template<typename T> int Oftsh<T>::nsd(int odmax, double sdmax)
{
    int res = 0;
    for(int k =0; k < Manip::nmon(m_oftsh_nvar, m_oftsh_order); k++) res += this->get_coef(k).nsd(odmax, sdmax);
    return res;
}


/**
 *  \brief  Evaluates the \f$ L_{\infty} \f$ norm of the current Oftsh object.
 */
template<typename T> double Oftsh<T>::linfnorm()
{
    double lin = 0.0;
    if(MODEL_TYPE == Csts::QBCP)
    {
        lin = cabs(this->get_coef(0).l1norm()+0.0*I);
        for(int k =1; k < Manip::nmon(m_oftsh_nvar, m_oftsh_order); k++)
        {
            if(lin < fabs(this->get_coef(k).l1norm())) lin = fabs(this->get_coef(k).l1norm());
        }
    }
    else if(MODEL_TYPE == Csts::CRTBP)
    {
        lin = cabs(this->get_coef(0).ofs_get_coef(0));
        for(int k =1; k < Manip::nmon(m_oftsh_nvar, m_oftsh_order); k++)
        {
           if(lin < cabs(this->get_coef(k).ofs_get_coef(0))) lin = cabs(this->get_coef(k).ofs_get_coef(0));
        }
    }
    return lin;
}

//----------------------------------------------------------------------------------------
// Functions
//----------------------------------------------------------------------------------------
/**
 * \fn template<typename T> Oftsh<T> operator + (Oftsh<T> const& a, Oftsh<T> const& b)
 * \brief An operator. Makes the sum a+b
 */
template<typename T> Oftsh<T> operator + (Oftsh<T> const& a, Oftsh<T> const& b)
{
    Oftsh<T> cop(a);
    cop+=b;
    return cop;
}

/**
 * \fn template<typename T> Oftsh<T> operator - (Oftsh<T> const& a, Oftsh<T> const& b)
 * \brief An operator. Makes the sub a-b
 */
template<typename T> Oftsh<T> operator - (Oftsh<T> const& a, Oftsh<T> const& b)
{
    Oftsh<T> cop(a);
    cop-=b;
    return cop;
}

/**
 *  \brief  Expected error on the product: \c this \f$ += c \times a \f$ at times t.
 *          Works only when a.m_oftsh_order = b.m_oftsh_order which is the default case.
 */
template <typename U> cdouble smult_expct_error(Oftsh<Ofs<U> > const& a, Ofs<U> const& c, U X[], double const& t)
{
    //Initialization
    cdouble result = 0.0+0.0*I;
    int *kv = (int*) calloc(a.get_order(), sizeof(int));
    kv[0] = a.get_order();
    for(int i=1; i< a.get_nvar(); i++) kv[i] = 0;
    U aux, bux;

    //Loop on all the monomials in a
    for (int i=0; i< Manip::nmon(a.get_nvar(), a.get_order()); i++)
    {
        bux = 1.0+0.0*I;
        for(int ii = 0; ii < a.get_nvar(); ii++)
        {
            aux = 1.0+0.0*I;
            if(kv[ii] != 0.0)
            {
                for(int j = 1; j <= kv[ii]; j++) aux*= X[ii];
            }
            bux *= aux;
        }
        //Updating the result
        result += bux*sprod_expct_error(a.get_coef(i), c, t);

        //updating the exponents in kv
        if(i< Manip::nmon(a.get_nvar(), a.get_order())-1)  Manip::prxkt(kv, a.get_nvar());
    }
    free(kv);
    return result;
}

/**
 *  \brief  Expected error on the product: \c this \f$ += c \times a \f$ at times t.
 *          Works only when a.m_oftsh_order = b.m_oftsh_order which is the default case. double version
 */
inline cdouble smult_expct_error(Oftsh<Ofsd > const& a, Ofsd const& c, double X[], double const& t)
{
    //Initialization
    cdouble result = 0.0+0.0*I;
    int *kv = (int*) calloc(a.get_order(), sizeof(int));
    kv[0] = a.get_order();
    for(int i=1; i< a.get_nvar(); i++) kv[i] = 0;
    double aux, bux;

    //Loop on all the monomials in a
    for (int i=0; i< Manip::nmon(a.get_nvar(), a.get_order()); i++)
    {
        bux = 1.0;
        for(int ii = 0; ii < a.get_nvar(); ii++)
        {
            aux = 1.0;
            if(kv[ii] != 0.0)
            {
                for(int j = 1; j <= kv[ii]; j++) aux*= X[ii];
            }
            bux *= aux;
        }
        //Updating the result
        result += bux*sprod_expct_error(a.get_coef(i), c, t);

        //updating the exponents in kv
        if(i< Manip::nmon(a.get_nvar(), a.get_order())-1)  Manip::prxkt(kv, a.get_nvar());
    }
    free(kv);
    return result;
}



//----------------------------------------------------------------------------------------
//Stream
//----------------------------------------------------------------------------------------
/**
 *  \brief  A stream operator
 */
template<typename T> std::ostream& operator << (std::ostream& stream, Oftsh<T> const& oftsh)
{
    int i,j;
    int k[oftsh.m_oftsh_nvar];
    k[0] = oftsh.m_oftsh_order;
    for(i=1; i<oftsh.m_oftsh_nvar; i++) k[i] = 0;

    stream << "#Homogeneous polynomial"    << endl;
    stream << "#Order:     " << oftsh.m_oftsh_order    << endl;
    stream << "#Variables: " << oftsh.m_oftsh_nvar    << endl;
    stream << "--------------------------" << endl;

    for (i=0; i< Manip::nmon(oftsh.m_oftsh_nvar, oftsh.m_oftsh_order); i++)
    {
        for(j=0; j<oftsh.m_oftsh_nvar; j++) stream <<  setw(2) << setiosflags(ios::right) <<  k[j] << " ";
        stream << std::showpos << setiosflags(ios::scientific)  << setprecision(15) << " " <<  oftsh.m_oftsh_coef[i] << std::noshowpos << endl;

        if(i< Manip::nmon(oftsh.m_oftsh_nvar, oftsh.m_oftsh_order)-1)  Manip::prxkt(k, oftsh.m_oftsh_nvar);
    }
    return stream;
}

/**
 *  \brief  A stream operator. Oftsh<complex double>. Specialization of
 *          template<typename T> std::ostream& operator << (std::ostream& stream, Oftsh<T> const& oftsh)
 */
template<> inline std::ostream& operator << (std::ostream& stream, Oftsh<complex double> const& oftsh)
{
    int i,j;
    int k[oftsh.m_oftsh_nvar];
    k[0] = oftsh.m_oftsh_order;
    for(i=1; i<oftsh.m_oftsh_nvar; i++) k[i] = 0;

    stream << "#Homogeneous polynomial"    << endl;
    stream << "#Order:     " << oftsh.m_oftsh_order    << endl;
    stream << "#Variables: " << oftsh.m_oftsh_nvar    << endl;
    stream << "--------------------------" << endl;

    for (i=0; i< Manip::nmon(oftsh.m_oftsh_nvar, oftsh.m_oftsh_order); i++)
    {
        for(j=0; j<oftsh.m_oftsh_nvar; j++) stream <<  setw(2) << setiosflags(ios::right) <<  k[j] << " ";
        stream << std::showpos << setiosflags(ios::scientific)  << setprecision(15) << " " <<  creal(oftsh.m_oftsh_coef[i]) << " " <<  cimag(oftsh.m_oftsh_coef[i]) << std::noshowpos << endl;

        if(i< Manip::nmon(oftsh.m_oftsh_nvar, oftsh.m_oftsh_order)-1)  Manip::prxkt(k, oftsh.m_oftsh_nvar);
    }
    return stream;
}
