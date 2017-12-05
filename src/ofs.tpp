//########################################################################################
// Implementation of the Ofs template class (this file is included in ofs.h)
//########################################################################################

/**
 * \file ofs.tpp
 * \brief Fourier series template class. See ofs.h for details.
 * \author BLB
 */


//----------------------------------------------------------------------------------------
//Create
//----------------------------------------------------------------------------------------
/**
 *  \brief Default constructor of the class Ofs<T>.
 */
template <typename T> Ofs<T>::Ofs()
{
    m_ofs_order = OFS_ORDER;
    m_ofs_coef = new T[2*m_ofs_order+1];
    this->zero(); //all coefficients to zero
}

/**
 *  \brief Constructor with a given order.
 */
template <typename T>
Ofs<T>::Ofs(const int newOrder)
{
    m_ofs_order = newOrder;
    m_ofs_coef = new T[2*m_ofs_order+1];
    this->zero(); //all coefficients to zero
}

/**
 *  \brief Constructor from a given Ofs object.
 */
template <typename T> Ofs<T>::Ofs(Ofs const& b)
{
    m_ofs_order = b.m_ofs_order;
    m_ofs_coef = new T[2*m_ofs_order+1];
    for(int i = -m_ofs_order ; i<= m_ofs_order; i++) m_ofs_coef[i+m_ofs_order] = b.m_ofs_coef[i+m_ofs_order];
}

//----------------------------------------------------------------------------------------
//Delete
//----------------------------------------------------------------------------------------
/**
 *  \brief Default destructor of the class Ofs<T>.
 */
template <typename T> Ofs<T>::~Ofs<T>()
{
    if(m_ofs_coef != NULL) delete m_ofs_coef;
    m_ofs_coef = 0;
}

//----------------------------------------------------------------------------------------
//Copy
//----------------------------------------------------------------------------------------
/**
 *  \brief  Copy from a given Ofs object (only the coefficients).
 */
template <typename T> Ofs<T>& Ofs<T>::ccopy(Ofs<T> const& b)
{
    if(m_ofs_order != b.m_ofs_order)
    {
        cout << "Erreur in Ofs<T>::ccopy: orders do not match. Nothing is done." << endl;
        return *this;
    }
    else
    {
        for(int i = -m_ofs_order ; i<= m_ofs_order; i++) m_ofs_coef[i+m_ofs_order] = b.m_ofs_coef[i+m_ofs_order];//this->set_coef(b.get_coef(i), i);
        return *this;
    }
}

/**
 *  \brief  Linked copy from a given Ofs object (exact same object is obtained).
 */
template <typename T> Ofs<T>& Ofs<T>::lcopy(Ofs<T> const& b)
{
    m_ofs_order = b.m_ofs_order;
    m_ofs_coef = b.m_ofs_coef;
    return *this;
}

//----------------------------------------------------------------------------------------
//Setters
//----------------------------------------------------------------------------------------
/**
 *  \brief Sets a coefficient at a given position in the series.
 */
template <typename T> void Ofs<T>::set_coef(T const& value, int const& pos)
{
    if(fabs(pos) <= m_ofs_order) m_ofs_coef[pos+m_ofs_order] = value;
    else cout << "Error in Ofs<T>::set_coef: position is out of scope. No coefficient is set." << endl;
}

/**
 *  \brief Sets a coefficient at a given position in the series (cdouble version)
 */
template <typename T> template <typename U> void Ofs<T>::set_coef(U const& value, int const& pos)
{
    if(fabs(pos) <= m_ofs_order) m_ofs_coef[pos+m_ofs_order] = value+0.0*I;
    else cout << "Error in Ofs<T>::set_coef: position is out of scope. No coefficient is set." << endl;
}

/**
 *  \brief Adds a coefficient at a given position in the series.
 */
template <typename T> void Ofs<T>::add_coef(T const& value, int const& pos)
{
    if(fabs(pos) <= m_ofs_order) m_ofs_coef[pos+m_ofs_order] += value;
    else cout << "Error in Ofs<T>::add_coef: position is out of scope\n" << endl;
}

/**
 *  \brief Sets a coefficient to all positions in the series.
 */
template <typename T> void Ofs<T>::set_all_coefs(T const& value)
{
    for(int pos = -m_ofs_order; pos <= m_ofs_order; pos++) this->set_coef(value, pos);
}

/**
 *  \brief Sets random coefficients to all positions in the series.
 */
template <typename T> void Ofs<T>::set_random_coefs()
{
    for(int pos = -m_ofs_order; pos <= m_ofs_order; pos++)
    {
        this->set_coef((double)rand()/(10.0*(fabs(pow(pos,7.0))+1)*RAND_MAX), pos);
    }
}

/**
 *  \brief Sets random cdouble coefficients to all positions in the series.
 */
template <> inline void Ofsc::set_random_coefs()
{
    //all coefficients between expect order zero.
    for(int pos = -m_ofs_order; pos <= m_ofs_order; pos++)
    {
        this->set_coef(  (double)rand()/(10.0*(fabs(pow(pos,7.0))+1)*RAND_MAX) +
                      I*(double)rand()/(10.0*(fabs(pow(pos,7.0))+1)*RAND_MAX), pos);
    }
    //Warning: order zero is set real! (important for precision).
    double rd = (double) rand();
    cdouble fac = (cdouble) (rd/(10.0*(fabs(pow(0,7.0))+1)*RAND_MAX) + 0.0I);
    this->set_coef(fac, 0);
}

//----------------------------------------------------------------------------------------
//Getters
//----------------------------------------------------------------------------------------
/**
 *  \brief  Gets the order of the series.
 */
template <typename T> int Ofs<T>::get_order() const
{
    return m_ofs_order;
}

/**
 *  \brief  Gets the pointer address of the Ofs object
 */
template <typename T> Ofs<T>* Ofs<T>::get_ptr() const
{
    return (Ofs<T>*) this;
}

/**
 *  \brief  Gets the coefficient at a given position.
 *
 *   If the position is out of scope, 0.0 is returned by default, and a warning is sent.
 */
template <typename T> T Ofs<T>::tfs_get_coef(int pos) const
{
    if(pos < 2*m_ofs_order+1) return m_ofs_coef[pos];
    else
    {
        cout << "Warning in Ofs<T>::get_coef: position " << pos << " is out of scope. 0.0 is returned by default." << endl;
        return (T) 0.0;
    }
}

/**
 *  \brief  Gets the coefficient at a given position.
 *
 *   If the position is out of scope, 0.0 is returned by default, and a warning is sent.
 */
template <> inline cdouble Ofsc::tfs_get_coef(int pos) const
{
    if(pos < 2*m_ofs_order+1)  return m_ofs_coef[pos];
    else
    {
        cout << "Warning in Ofs<T>::get_coef: position " << pos << " is out of scope. 0.0 is returned by default." << endl;
        return 0.0+0.0*I;
    }
}

/**
 *  \brief  Gets the coefficient at a given position.
 *
 *   If the position is out of scope, 0.0 is returned by default, and a warning is sent.
 */
template <typename T> T Ofs<T>::ofs_get_coef(int pos) const
{
    if(fabs(pos) <= m_ofs_order)  return m_ofs_coef[pos+m_ofs_order];
    else
    {
        cout << "Warning in Ofs<T>::get_coef: position " << pos << " is out of scope. 0.0 is returned by default." << endl;
        return (T) 0.0;
    }
}

/**
 *  \brief  Gets the coefficient at a given position. cdouble case
 */
template <> inline cdouble Ofsc::ofs_get_coef(int pos) const
{
    if(fabs(pos) <= m_ofs_order)  return m_ofs_coef[pos+m_ofs_order];
    else
    {
        cout << "Warning in Ofs<T>::get_coef: position " << pos << " is out of scope. 0.0 is returned by default." << endl;
        return  0.0+0.0*I;
    }
}

/**
 *  \brief  Computes the maximum coefficient in norm.
 */
template <typename T> void Ofs<T>::get_coef_max_norm(double maxAbs[]) const
{
    int n = this->get_order();
    maxAbs[0] = cabs(m_ofs_coef[-n+m_ofs_order]);//cabs(this->get_coef(-n));
    maxAbs[1] = -n;

    //Loop on all the coefficient but -n
    for(int i = -n+1; i <= n; i++)
    {
        if(cabs(m_ofs_coef[i+m_ofs_order]) > maxAbs[0])
        {
            maxAbs[0] = cabs(m_ofs_coef[i+m_ofs_order]);
            maxAbs[1] = i;
        }
    }
}

/**
 *  \brief  Computes the maximum coefficient in norm.
 */
template <typename T> double Ofs<T>::get_coef_max_norm() const
{
    int n = this->get_order();
    double maxAbs = cabs(m_ofs_coef[-n+m_ofs_order]);//cabs(this->get_coef(-n));

    //Loop on all the coefficient but -n
    for(int i = -n+1; i <= n; i++)
    {
        if(cabs(m_ofs_coef[i+m_ofs_order]) > maxAbs)
        {
            maxAbs =cabs(m_ofs_coef[i+m_ofs_order]);
        }
    }

    return maxAbs;
}

/**
 *  \brief  Is the Ofs object equal to zero, at order ofs_order?
 */
template <typename T> bool Ofs<T>::is_null(const int ofs_order) const
{
    for(int i = -min(ofs_order, m_ofs_order); i <= min(ofs_order, m_ofs_order); i++)
    {
        if(cabs(ofs_get_coef(i)) != 0.0) return false;
    }
    return true;
}

//----------------------------------------------------------------------------------------
//Zeroing
//----------------------------------------------------------------------------------------
/**
 *  \brief  Sets all coefficients to zero.
 */
template <typename T> void Ofs<T>::zero()
{
    for(int i = -m_ofs_order ; i<= m_ofs_order; i++) m_ofs_coef[i+m_ofs_order] = 0.0;
}

/**
 *  \brief  Sets all coefficients to zero. cdouble case
 */
template <> inline void Ofsc::zero()
{
    for(int i = -m_ofs_order ; i<= m_ofs_order; i++) m_ofs_coef[i+m_ofs_order] = 0.0+0.0*I;
}


//----------------------------------------------------------------------------------------
// Frequency domain <--> Time domain
//----------------------------------------------------------------------------------------
/**
 *  \brief  From Frequency domain to time domain.
 */
template <typename T> void Ofs<T>::tfs_from_ofs(Ofs<T> const& a)
{
    //---------------------
    //if wrong orders, stop
    //---------------------
    if(m_ofs_order < a.get_order())
    {
        cout << "tfs_from_ofs. Wrong orders. Stop." << endl;
        return;
    }

    int N = 2*m_ofs_order+1;
    //---------------------
    //Copy the coefficients in m_ofs_coef
    //---------------------
    for(int i = 0; i < N; i++) m_ofs_coef[i] = ((Ofs<T>)a).evaluate(i*2*M_PI/((double)N));
}

/**
 *  \brief  Inline from Frequency domain to time domain.
 */
template <typename T> void Ofs<T>::tfs_from_ofs_inline(Ofs<T>& temp)
{
    //---------------------
    //if wrong orders, stop
    //---------------------
    if(m_ofs_order != temp.get_order())
    {
        cout << "tfs_from_ofs. Wrong orders. Stop." << endl;
        return;
    }

    //---------------------
    //Copy in temp
    //---------------------
    temp.ccopy(*this);

    int N = 2*m_ofs_order+1;
    //---------------------
    //Copy the coefficients in m_ofs_coef
    //---------------------
    for(int i = 0; i < N; i++) m_ofs_coef[i] = temp.evaluate(i*2*M_PI/((double)N));
}

/**
 *  \brief  From Time domain to Frequency domain.
 */
template <typename T> void Ofs<T>::tfs_to_ofs(Ofs<T> const& a)
{
    //---------------------
    // If wrong orders, stop
    //---------------------
    if(m_ofs_order > a.get_order())
    {
        cout << "tfs_to_ofs. Wrong orders. Stop." << endl;
        return;
    }

    //---------------------
    // FFT structures
    //---------------------
    int N = 2*a.get_order()+1;
    gsl_vector_complex *data_fft          = gsl_vector_complex_calloc(N);
    gsl_fft_complex_wavetable *wavetable  = gsl_fft_complex_wavetable_alloc (N);
    gsl_fft_complex_workspace *workspace  = gsl_fft_complex_workspace_alloc (N);

    //---------------------
    //Copy a in data_fft
    //---------------------
    for(int i = 0; i < N; i++) gsl_vector_complex_set(data_fft, i, gslc_complex(creal(a.m_ofs_coef[i]), cimag(a.m_ofs_coef[i])));

    //---------------------
    //FFT
    //---------------------
    gsl_fft_complex_forward (data_fft->data, 1, N, wavetable, workspace);

    //---------------------
    //Order 0
    //---------------------
    this->set_coef(+GSL_REAL(gsl_vector_complex_get(data_fft, 0))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_fft, 0))/(double)N,  0);
    //Version without set_coef
    //this->m_ofs_coef[m_ofs_order] = +GSL_REAL(gsl_vector_complex_get(data_fft, 0))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_fft, 0))/(double)N;

    //---------------------
    //Order n
    //---------------------
    for(int i = 1; i<= m_ofs_order; i++)
    {
        //Negative frequecies
        this->set_coef(+GSL_REAL(gsl_vector_complex_get(data_fft, N-i))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_fft, N-i))/(double)N, -i);
        //Positive frequencies
        this->set_coef(+GSL_REAL(gsl_vector_complex_get(data_fft, i))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_fft, i))/(double)N,  i);
        //Version without set_coef
        //this->m_ofs_coef[m_ofs_order-i] = +GSL_REAL(gsl_vector_complex_get(data_fft, N-i))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_fft, N-i))/(double)N;
        //this->m_ofs_coef[m_ofs_order+i] = +GSL_REAL(gsl_vector_complex_get(data_fft, i))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_fft, i))/(double)N;
    }

    //------------------------
    // Free
    //------------------------
    gsl_vector_complex_free(data_fft);
    gsl_fft_complex_wavetable_free (wavetable);
    gsl_fft_complex_workspace_free (workspace);
}

/**
 *  \brief  Inline from Time domain to Frequency domain.
 */
template <typename T> void Ofs<T>::tfs_to_ofs_inline()
{
    //---------------------
    // FFT structures
    //---------------------
    int N = 2*m_ofs_order+1;
    gsl_vector_complex *data_fft          = gsl_vector_complex_calloc(N);
    gsl_fft_complex_wavetable *wavetable  = gsl_fft_complex_wavetable_alloc (N);
    gsl_fft_complex_workspace *workspace  = gsl_fft_complex_workspace_alloc (N);

    //---------------------
    //Copy a in data_fft
    //---------------------
    for(int i = 0; i < N; i++) gsl_vector_complex_set(data_fft, i, gslc_complex(creal(this->m_ofs_coef[i]), cimag(this->m_ofs_coef[i])));

    //---------------------
    //FFT
    //---------------------
    gsl_fft_complex_forward (data_fft->data, 1, N, wavetable, workspace);

    //---------------------
    //Order 0
    //---------------------
    this->set_coef(+GSL_REAL(gsl_vector_complex_get(data_fft, 0))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_fft, 0))/(double)N,  0);
    //Version without set_coef
    //this->m_ofs_coef[m_ofs_order] = +GSL_REAL(gsl_vector_complex_get(data_fft, 0))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_fft, 0))/(double)N;

    //---------------------
    //Order n
    //---------------------
    for(int i = 1; i<= m_ofs_order; i++)
    {
        //Negative frequecies
        this->set_coef(+GSL_REAL(gsl_vector_complex_get(data_fft, N-i))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_fft, N-i))/(double)N, -i);
        //Positive frequencies
        this->set_coef(+GSL_REAL(gsl_vector_complex_get(data_fft, i))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_fft, i))/(double)N,  i);
        //Version without set_coef
        //this->m_ofs_coef[m_ofs_order-i] = +GSL_REAL(gsl_vector_complex_get(data_fft, N-i))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_fft, N-i))/(double)N;
        //this->m_ofs_coef[m_ofs_order+i] = +GSL_REAL(gsl_vector_complex_get(data_fft, i))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_fft, i))/(double)N;
    }

    //------------------------
    // Free
    //------------------------
    gsl_vector_complex_free(data_fft);
    gsl_fft_complex_wavetable_free (wavetable);
    gsl_fft_complex_workspace_free (workspace);
}

/**
 *  \brief  From Time domain to Frequency domain. Alternative version with Fortran code, for real values only
 */
template <typename T> void Ofs<T>::tfs_to_ofs_F(Ofs<T>& a)
{
    //---------------------
    // If wrong orders, stop
    //---------------------
    if(m_ofs_order > a.get_order())
    {
        cout << "tfs_to_ofs. Wrong orders. Stop." << endl;
        return;
    }

    //---------------------
    //FFT
    //---------------------
    int N = 2*m_ofs_order+1;
    int M = m_ofs_order;
    double CSF[M+1], SIF[M+1], F[N];
    for(int i = 0; i < N; i++) F[i]  = creal(a.m_ofs_coef[i]);
    foun_(F, &N, &M, CSF, SIF);


    //---------------------
    //Order 0
    //---------------------
    this->set_coef(CSF[0],  0);

    //---------------------
    //Order n
    //---------------------
    for(int i = 1; i<= m_ofs_order; i++)
    {
        //Negative frequecies
        this->set_coef(0.5*(CSF[i] + I*SIF[i]), -i);
        //Positive frequencies
        this->set_coef(0.5*(CSF[i] - I*SIF[i]),  i);
    }

}

/**
 *  \brief  From Time domain to Frequency domain. Version with externalized GSL structures.
 */
template <typename T> void Ofs<T>::tfs_to_ofs(Ofs<T>& a, gsl_vector_complex *data_fft,
                                              gsl_fft_complex_wavetable *wavetable,
                                              gsl_fft_complex_workspace *workspace)
{
    //---------------------
    // If wrong orders, stop
    //---------------------
    if(m_ofs_order != a.get_order())
    {
        cout << "tfs_from_ofs. Wrong orders. Stop." << endl;
        return;
    }

    //---------------------
    // FFT structures
    //---------------------
    int N = 2*m_ofs_order+1;

    //---------------------
    //Copy a in data_fft
    //---------------------
    for(int i = 0; i < N; i++) gsl_vector_complex_set(data_fft, i, gslc_complex(creal(a.m_ofs_coef[i]), cimag(a.m_ofs_coef[i])));

    //---------------------
    //FFT
    //---------------------
    gsl_fft_complex_forward (data_fft->data, 1, N, wavetable, workspace);

    //---------------------
    //Order 0
    //---------------------
    this->set_coef(+GSL_REAL(gsl_vector_complex_get(data_fft, 0))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_fft, 0))/(double)N,  0);
    //Version without set_coef
    //this->m_ofs_coef[m_ofs_order] = +GSL_REAL(gsl_vector_complex_get(data_fft, 0))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_fft, 0))/(double)N;

    //---------------------
    //Order n
    //---------------------
    for(int i = 1; i<= m_ofs_order; i++)
    {
        //Negative frequecies
        this->set_coef(+GSL_REAL(gsl_vector_complex_get(data_fft, N-i))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_fft, N-i))/(double)N, -i);
        //Positive frequencies
        this->set_coef(+GSL_REAL(gsl_vector_complex_get(data_fft, i))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_fft, i))/(double)N,  i);
        //Version without set_coef
        //this->m_ofs_coef[m_ofs_order-i] = +GSL_REAL(gsl_vector_complex_get(data_fft, N-i))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_fft, N-i))/(double)N;
        //this->m_ofs_coef[m_ofs_order+i] = +GSL_REAL(gsl_vector_complex_get(data_fft, i))/(double)N+I*GSL_IMAG(gsl_vector_complex_get(data_fft, i))/(double)N;
    }
}


//----------------------------------------------------------------------------------------
//Operators (+=, -=, ...)
//----------------------------------------------------------------------------------------
/**
 *  \brief  An operator. Constructor from a given Ofs object (only the coefficients).
 */
template <typename T> Ofs<T>& Ofs<T>::operator = (Ofs<T> const& b)
{
    if(this != &b)
    {
        m_ofs_order = b.m_ofs_order;
        if(m_ofs_coef != NULL) delete m_ofs_coef;
        m_ofs_coef = new T[2*m_ofs_order+1];
        for(int i = -m_ofs_order ; i<= m_ofs_order; i++) m_ofs_coef[i+m_ofs_order] = b.m_ofs_coef[i+m_ofs_order]; //this->set_coef(b.get_coef(i), i);
    }
    return *this; //same object if returned
}

/**
 *  \brief  An operator. Constructor from a given coefficient. The returned object is of order zero.
 */
template <typename T> Ofs<T>& Ofs<T>::operator  = (T const& coef0)
{
    m_ofs_order = 0;
    m_ofs_coef = new T[2*m_ofs_order+1];
    m_ofs_coef[m_ofs_order] = coef0;
    return *this;
}

/**
 * \brief An operator. Compares two Ofs objects
 */
template <typename T> bool Ofs<T>::isEqual(Ofs<T> const& b) const
{
    if(m_ofs_order != b.m_ofs_order) return false;
    else
    {
        bool result = true;
        for(int i = 0 ; i< 2*m_ofs_order + 1; i++) result = result&&(m_ofs_coef[i] == b.m_ofs_coef[i]);
        return result;
    }

}

/**
 *  \brief  An operator. Adds all coefficients term by term  from a given Ofs object.
 *
 *   Allows b.m_ofs_order != m_ofs_order.
 */
template <typename T>  Ofs<T>& Ofs<T>::operator += (Ofs<T> const& b)
{
    if(b.m_ofs_order > m_ofs_order) //if b.m_ofs_order > m_ofs_order, a new array of coefficients must be set
    {
        //Copy m_ofs_coef into temporary array
        T temp[2*m_ofs_order+1];
        for(int i = 0 ; i< 2*m_ofs_order + 1; i++) temp[i] = m_ofs_coef[i];
        //Recreate a good array
        delete m_ofs_coef;
        m_ofs_coef = new T[2*b.m_ofs_order+1];
        //Store the coefficients again
        for(int i = -m_ofs_order ; i<= m_ofs_order; i++) m_ofs_coef[i+b.m_ofs_order] = temp[i+m_ofs_order];
        m_ofs_order = b.m_ofs_order;
    }

    //Adding the coefficients
    for(int i = -b.m_ofs_order ; i <= b.m_ofs_order ; i++) m_ofs_coef[i+m_ofs_order] += b.m_ofs_coef[i+b.m_ofs_order];
    return *this;
}

/**
 *  \brief  An operator. Subtracts all coefficients term by term  from a given Ofs object.
 *
 *   Allows b.m_ofs_order != m_ofs_order.
 */
template <typename T>  Ofs<T>& Ofs<T>::operator -= (Ofs<T> const& b)
{
    if(b.m_ofs_order > m_ofs_order) //if b.m_ofs_order > m_ofs_order, a new array of coefficients must be set
    {
        //Copy m_ofs_coef into temporary array
        T temp[2*m_ofs_order+1];
        for(int i = 0 ; i< 2*m_ofs_order + 1; i++) temp[i] = m_ofs_coef[i];
        //Recreate a good array
        delete m_ofs_coef;
        m_ofs_coef = new T[2*b.m_ofs_order+1];
        //Store the coefficients again
        for(int i = -m_ofs_order ; i<= m_ofs_order; i++) m_ofs_coef[i+b.m_ofs_order] = temp[i+m_ofs_order];
        m_ofs_order = b.m_ofs_order;
    }

    //Adding the coefficients
    for(int i = -b.m_ofs_order ; i <= b.m_ofs_order ; i++) m_ofs_coef[i+m_ofs_order] -= b.m_ofs_coef[i+b.m_ofs_order];
    return *this;
}

/**
 *  \brief  An operator. Multiplies all coefficients by a given \c T coefficient.
 */
template <typename T>  Ofs<T>& Ofs<T>::operator *= (T const& c)
{
    for(int i=0; i<2*m_ofs_order+1; i++) m_ofs_coef[i] *= c;
    return *this;
}

/**
 *  \brief  An operator. Divides all coefficients by a given \c T coefficient.
 */
template <typename T>  Ofs<T>& Ofs<T>::operator /= (T const& c)
{
    for(int i=0; i<2*m_ofs_order+1; i++) m_ofs_coef[i] /= c;
    return *this;
}


//----------------------------------------------------------------------------------------
//Operators (sprod, smult, ...)
//----------------------------------------------------------------------------------------
/**
 *  \brief Conjugates the Ofs object.
 */
template <typename T> void Ofs<T>::conjugate()
{
    Ofs ofs_temp(*this);
    for(int k = -m_ofs_order; k <=m_ofs_order; k++) m_ofs_coef[k+m_ofs_order] = ofs_temp.m_ofs_coef[-k+m_ofs_order];//this->set_coef(ofs_temp.get_coef(-k), k);
}

/**
 *  \brief Conjugates the Ofsc object (the conjugate of each term is taken).
 */
template <> inline void Ofs<cdouble >::conjugate()
{
    Ofs ofs_temp(*this);
    for(int k = -m_ofs_order; k <=m_ofs_order; k++) m_ofs_coef[k+m_ofs_order] = conj(ofs_temp.m_ofs_coef[-k+m_ofs_order]);//this->set_coef(conj(ofs_temp.get_coef(-k)), k);
}

/**
 *  \brief  An operation. Adds the product: \c this \f$ += a \times b \f$.

    Notes:
    1. \c this \f$ += a \times b \f$. with new m_ofs_order = max(m_ofs_order, a.m_ofs_order, b.m_ofs_order): the sum is truncated at max(a.get_order(),b.get_order()).
    2. WARNING: Need improvement: for n <= a.get_order(), b.get_order(), some products are out of scope in: this->add_coef(a.get_coef(p)*b.get_coef(n-p), n). The get_coef function set these coefficients to zero, which guarantees the good result. However, unecessary product are made. psup and pinf must be redefined.
    3. Works fine when a.m_ofs_order = b.m_ofs_order which is the default case.
 */
template <typename T>  void Ofs<T>::ofs_sprod(Ofs<T> const& a, Ofs<T> const& b)
{
    int psup, pinf;
    //Product
    for(int n=-m_ofs_order ; n<= m_ofs_order; n++)
    {
        //psup = n>0? m_ofs_order: n+m_ofs_order;//min(n+m_ofs_order,  m_ofs_order);
        //pinf = n>0? n-m_ofs_order: -m_ofs_order;//max(n-m_ofs_order, -m_ofs_order);
        psup = min(n+m_ofs_order,  m_ofs_order);
        pinf = max(n-m_ofs_order, -m_ofs_order);
        for(int p=pinf; p<= psup; p++) m_ofs_coef[n+m_ofs_order] += a.m_ofs_coef[p+m_ofs_order]*b.m_ofs_coef[n-p+m_ofs_order];
    }
}

/**
 *  \brief  An operation. Adds the product: \c this \f$ += a \times b \times c \f$.
 */
template <typename T>  void Ofs<T>::ofs_smprod_t(Ofs<T> const& a, Ofs<T> const& b, Ofs<T> const& c, Ofs<T>& temp)
{
    //temp = a*b
    temp.ofs_prod(a,b);
    //this = temp*c = a*b*c
    this->ofs_sprod(temp, c);
}

/**
 *  \brief  An operation. Adds the product: \c this \f$  += m a \times b \f$.

    Notes:
    1. \c this \f$ += a \times b \f$. with new m_ofs_order = max(m_ofs_order, a.m_ofs_order, b.m_ofs_order): the sum is truncated at max(a.get_order(),b.get_order()).
    2. WARNING: Need improvement: for n <= a.get_order(), b.get_order(), some products are out of scope in: this->add_coef(a.get_coef(p)*b.get_coef(n-p), n). The get_coef function set these coefficients to zero, which guarantees the good result. However, unecessary product are made. psup and pinf must be redefined.
    3. Works fine when a.m_ofs_order = b.m_ofs_order which is the default case.
 */
template <typename T>  void Ofs<T>::ofs_smprod(Ofs<T> const& a, Ofs<T> const& b, T const& m)
{
    int psup, pinf;
    //Product
    for(int n=-m_ofs_order ; n<= m_ofs_order; n++)
    {
        psup = min(n+m_ofs_order,  m_ofs_order);
        pinf = max(n-m_ofs_order, -m_ofs_order);
        for(int p=pinf; p<= psup; p++)
        {
            //indix
            m_ofs_coef[n+m_ofs_order] += m*a.m_ofs_coef[p+m_ofs_order]*b.m_ofs_coef[n-p+m_ofs_order];
        }
    }
}

/**
 *  \brief  An operation. Sets the product: \c this \f$ = a \times b \f$.
 *
 *  Notes: see smprod.
 */
template <typename T> void Ofs<T>::ofs_prod(Ofs<T> const& a, Ofs<T> const& b)
{
    int psup, pinf;
    //Product

    for(int n=-m_ofs_order ; n<= m_ofs_order; n++)
    {
        //psup = n>0? m_ofs_order: n+m_ofs_order;//min(n+m_ofs_order,  m_ofs_order);
        //pinf = n>0? n-m_ofs_order: -m_ofs_order;//max(n-m_ofs_order, -m_ofs_order);
        psup = min(n+m_ofs_order,  m_ofs_order);
        pinf = max(n-m_ofs_order, -m_ofs_order);
        m_ofs_coef[n+m_ofs_order] = 0.0; //reset

       for(int p=pinf; p<= psup; p++) m_ofs_coef[n+m_ofs_order] += a.m_ofs_coef[p+m_ofs_order]*b.m_ofs_coef[n-p+m_ofs_order];
    }
}

/**
 *  \brief  An operation. Sets the product: \c this \f$ = a \times b \f$. cdouble case
 *
 *  Notes: see smprod.
 */
template <> inline void Ofsc::ofs_prod(Ofsc const& a, Ofsc const& b)
{
    int psup, pinf;
    //Product

    for(int n=-m_ofs_order ; n<= m_ofs_order; n++)
    {
        //psup = n>0? m_ofs_order: n+m_ofs_order;//min(n+m_ofs_order,  m_ofs_order);
        //pinf = n>0? n-m_ofs_order: -m_ofs_order;//max(n-m_ofs_order, -m_ofs_order);
        psup = min(n+m_ofs_order,  m_ofs_order);
        pinf = max(n-m_ofs_order, -m_ofs_order);
        m_ofs_coef[n+m_ofs_order] = 0.0+0.0*I; //reset

       for(int p=pinf; p<= psup; p++) m_ofs_coef[n+m_ofs_order] += a.m_ofs_coef[p+m_ofs_order]*b.m_ofs_coef[n-p+m_ofs_order];
    }
}

/**
 *  \brief  An operation. Sets the product: \c this \f$ = m a \times b \f$.
 *
 *  Notes: see smprod.
 */
template <typename T>  void Ofs<T>::ofs_mprod(Ofs<T> const& a, Ofs<T> const& b, T const& m)
{
    int psup, pinf;
    //Product
    for(int n=-m_ofs_order ; n<= m_ofs_order; n++)
    {
        psup = min(n+m_ofs_order,  m_ofs_order);
        pinf = max(n-m_ofs_order, -m_ofs_order);
        m_ofs_coef[n+m_ofs_order] = 0.0; //reset
        for(int p=pinf; p<= psup; p++)
        {
            //indix
            m_ofs_coef[n+m_ofs_order] += m*a.m_ofs_coef[p+m_ofs_order]*b.m_ofs_coef[n-p+m_ofs_order];
        }
    }
}

/**
 *  \brief  An operation. Sets the product: \c this \f$ = m a \times b \f$. cdouble case
 *
 *  Notes: see smprod.
 */
template <>  inline  void Ofsc::ofs_mprod(Ofsc const& a, Ofsc const& b, cdouble const& m)
{
    int psup, pinf;
    //Product
    for(int n=-m_ofs_order ; n<= m_ofs_order; n++)
    {
        psup = min(n+m_ofs_order,  m_ofs_order);
        pinf = max(n-m_ofs_order, -m_ofs_order);
        m_ofs_coef[n+m_ofs_order] = 0.0+0.0*I; //reset
        for(int p=pinf; p<= psup; p++)
        {
            //indix
            m_ofs_coef[n+m_ofs_order] += m*a.m_ofs_coef[p+m_ofs_order]*b.m_ofs_coef[n-p+m_ofs_order];
        }
    }
}

/**
 *  \brief  An operation. Adds the product: \c this \f$  += m a \f$.
 *
 *  Note: can be used in place.
 */
template<typename T> void Ofs<T>::ofs_smult(Ofs<T> const& a, T const& c)
{
    if(m_ofs_order != a.m_ofs_order)
    {
        cout << "Error using smult: the m_ofs_order of variables does not match. Initial Ofs<T> is returned" << endl;
        return;
    }
    else
    {
        //Sum
        for(int i = -m_ofs_order; i <= m_ofs_order; i++) m_ofs_coef[i+m_ofs_order] += c*a.m_ofs_coef[i+m_ofs_order];//add_coef(c*a.get_coef(i), i);
    }
}

/**
 *  \brief  An operation. Adds the product: \c this \f$  += m a \f$ at a certain order eff_order.
 *
 *  Note: can be used in place.
 */
template<typename T> void Ofs<T>::ofs_smult(Ofs<T> const& a, T const& c, int eff_order)
{
    //Sum
    for(int i = -eff_order; i <= eff_order; i++)
    {
        add_coef(c*a.ofs_get_coef(i), i);
    }
}

/**
 *  \brief  An operation. Sets the product: \c this \f$  = m a \f$.
 *
 *  Note: can be used in place.
 */
template<typename T> void Ofs<T>::ofs_mult(Ofs<T> const& a, T const& c)
{
    if(m_ofs_order != a.m_ofs_order)
    {
        cout << "Error using smult: the m_ofs_order of variables does not match. Initial Ofs<T> is returned" << endl;
        return;
    }
    else
    {
        //Sum
        for(int i=0; i<2*m_ofs_order+1; i++) m_ofs_coef[i] = c*a.m_ofs_coef[i];
    }
}

/**
 *  \brief  An operation. Sets the sum-product: \c this \f$  = m_a a + m_b b \f$.
 *
 *  Note: can be used in place.
 */
template<typename T> void Ofs<T>::ofs_fsum(Ofs<T> const& a, T const& ma, Ofs<T> const& b, T const& mb)
{
    if(m_ofs_order != a.m_ofs_order ||  m_ofs_order != b.m_ofs_order)
    {
        cout << "Error using fsum: the m_ofs_order does not match. Initial Ofs<T> is returned" << endl;
        return;
    }
    else
    {
        for(int i=-m_ofs_order; i<=m_ofs_order; i++)
        {
            m_ofs_coef[i+m_ofs_order] = ma*a.m_ofs_coef[i+m_ofs_order] + mb*b.m_ofs_coef[i+m_ofs_order];//set_coef(ma*a.get_coef(i)+mb*b.get_coef(i), i);
        }
    }
}

/**
 *  \brief  An operation. Performs the expansion this \f$ = a^{\alpha} \f$ when \c a is of the form \c this \f$ a = 1 + d \f$ with \f$ |d|_1 << 1 \f$.
 *
 *  Note: works well with at least a factor \f$ 1e-4 \f$ between the order zero and the order one.
 */
template<typename T> void Ofs<T>::ofs_epspow(Ofs<T> const& a, T const& alpha)
{
    //set to zero
    this->zero();
    //order 0
    this->set_coef(1.0, 0);
    //order 1
    Ofs<T> a0(a);
    a0.set_coef(0.0,0); //a0 = a without the first order
    Ofs<T> ak(a0);
    Ofs<T> akc(a0);
    int n_order_fourier = a.get_order();
    double facinv = 1;
    T m_ofs_coef = alpha;
    //this = this + alpha*a
    this->ofs_smult(ak, facinv*m_ofs_coef);
    for(int k = 2; k<= n_order_fourier; k++)
    {
        //1/fac(k)
        facinv*=1.0/k;
        //alpha*(alpha-1)...*(alpha-k+1)
        m_ofs_coef*= alpha - k + 1;
        //akc = ak;
        akc.ccopy(ak);
        //ak = akc*a0
        ak.ofs_prod(akc, a0);
        //this += facinv*m_ofs_coef*ak
        this->ofs_smult(ak, facinv*m_ofs_coef);
    }
}

/**
 *  \brief An operation. Performs the expansion this \f$ = a^{\alpha} \f$ using inverse FFT, the power function in time domain, and finally direct FFT.
 **/
template<typename T> void Ofs<T>::ofs_pows(Ofs<T> const& a, T const& alpha)
{
    this->tfs_from_ofs(a);
    this->tfs_pows(alpha);
    this->tfs_to_ofs_inline();
}

//----------------------------------------------------------------------------------------
// TFS operations
//----------------------------------------------------------------------------------------
//-----------------
// pows
//-----------------
/**
 *  \brief  An operation. Performs the power this = a^alpha in time domain.
 */
template<typename T> void Ofs<T>::tfs_pows(Ofs<T> const& a, T const& alpha)
{
    for(int i = 0; i < 2*m_ofs_order+1; i++) m_ofs_coef[i] = cpow(a.m_ofs_coef[i], alpha);
}

/**
 *  \brief  An operation. Performs the power this = this^alpha in time domain.
 */
template<typename T> void Ofs<T>::tfs_pows(T const& alpha)
{
    for(int i = 0; i < 2*m_ofs_order+1; i++) m_ofs_coef[i] = cpow(m_ofs_coef[i], alpha);
}

//-----------------
// sprod
//-----------------
/**
 *  \brief  An operation. Adds the product: \c this \f$ += a \times b \f$ in time domain.

    Notes:
    1. \c this \f$ += a \times b \f$. with new m_ofs_order = max(m_ofs_order, a.m_ofs_order, b.m_ofs_order): the sum is truncated at max(a.get_order(),b.get_order()).
    2. WARNING: Need improvement: for n <= a.get_order(), b.get_order(), some products are out of scope in: this->add_coef(a.get_coef(p)*b.get_coef(n-p), n). The get_coef function set these coefficients to zero, which guarantees the good result. However, unecessary product are made. psup and pinf must be redefined.
    3. Works fine when a.m_ofs_order = b.m_ofs_order which is the default case.
 */
template <typename T>  void Ofs<T>::tfs_sprod(Ofs<T> const& a, Ofs<T> const& b)
{
    for(int k=0 ; k< 2*m_ofs_order+1; k++) m_ofs_coef[k] += a.m_ofs_coef[k]*b.m_ofs_coef[k];
}

/**
 *  \brief  An operation. Adds the product: \c this \f$ += a \times b \times c \f$ in time domain.
 */
template <typename T>  void Ofs<T>::tfs_smprod_t(Ofs<T> const& a, Ofs<T> const& b, Ofs<T> const& c)
{
    for(int k=0 ; k< 2*m_ofs_order+1; k++) m_ofs_coef[k] += a.m_ofs_coef[k]*b.m_ofs_coef[k]*c.m_ofs_coef[k];
}

/**
 *  \brief  An operation. Adds the product: \c this \f$ += a \times b \times c \f$ in time domain.
 */
template <typename T> template<typename U> void Ofs<T>::tfs_smprod_tu(Ofs<T> const& a, Ofs<T> const& b, Ofs<T> const& c, U const& m)
{
    for(int k=0 ; k< 2*m_ofs_order+1; k++) m_ofs_coef[k] += m*a.m_ofs_coef[k]*b.m_ofs_coef[k]*c.m_ofs_coef[k];
}

/**
 *  \brief  An operation. Adds the product: \c this \f$ += a \times b \times c \f$ in time domain.
 */
template <typename T> void Ofs<T>::tfs_smprod(Ofs<T> const& a, Ofs<T> const& b,  T const& m)
{
    for(int k=0 ; k< 2*m_ofs_order+1; k++) m_ofs_coef[k] += m*a.m_ofs_coef[k]*b.m_ofs_coef[k];
}

//----------------------------------------------------------------------------------------
//Derivation
//----------------------------------------------------------------------------------------
/**
 *  \brief  An operation. Set the time derivative of object \c a with pulsation \f$ \omega = n \f$, so that \c this \f$ = \dot{a} \f$ in frequency domain.
 */
template<typename T> void Ofs<T>::dot(Ofs<T> const& a, double const& n)
{
    //d(a_k)/dt = k*n*I*a_k
    for(int k=-m_ofs_order; k<= m_ofs_order; k++) m_ofs_coef[k+m_ofs_order] = k*n*I*a.m_ofs_coef[k+m_ofs_order];
}

/**
 *  \brief  An operation. Applies the time derivative with pulsation \f$ \omega = n \f$ directly to this in frequency domain.
 */
template<typename T> void Ofs<T>::dot(double const& n)
{
    //d(a_k)/dt = k*n*I*a_k
    for(int k=-m_ofs_order; k<= m_ofs_order; k++) m_ofs_coef[k+m_ofs_order] = k*n*I*m_ofs_coef[k+m_ofs_order];
}

//----------------------------------------------------------------------------------------
// Functions (+, -,...)
//----------------------------------------------------------------------------------------
/**
 * \fn template<typename T> bool   operator == (Ofs<T> const& a, Ofs<T> const& b)
 * \brief An operator. Compares two Ofs objects
 */
template <typename T>  bool operator == (Ofs<T> const& a, Ofs<T> const& b)
{
    return a.isEqual(b);
}

/**
 * \fn template<typename T> Ofs<T> operator + (Ofs<T> const& a, Ofs<T> const& b)
 * \brief An operator. Makes the sum a+b
 */
template <typename T>  Ofs<T> operator + (Ofs<T> const& a, Ofs<T> const& b)
{
    Ofs<T> cop(a);
    cop+=b;
    return cop;
}

/**
 * \fn template<typename T> Ofs<T> operator - (Ofs<T> const& a, Ofs<T> const& b)
 * \brief An operator. Makes the sum a-b
 */
template <typename T>  Ofs<T> operator - (Ofs<T> const& a, Ofs<T> const& b)
{
    Ofs<T> cop(a);
    cop-=b;
    return cop;
}

/**
 * \fn template<typename T> Ofs<T> operator - (Ofs<T> const& b)
 * \brief Returns -b at order b.m_ofs_order
 */
template <typename T>  Ofs<T> operator - (Ofs<T> const& b)
{
    Ofs<T> ofs(b.get_order());
    for(int i = -b.get_order(); i<= b.get_order() ; i++) ofs.set_coef(-b.ofs_get_coef(i),i);
    return ofs;
}

/**
 * \fn template<typename T> Ofs<T> operator - (Ofs<T> const& b). cdouble case
 * \brief Returns -b at order b.m_ofs_order
 */
template <>  inline  Ofsc operator - (Ofsc const& b)
{
    Ofsc ofs(b.get_order());
    for(int i = -b.get_order(); i<= b.get_order() ; i++) ofs.set_coef(0.0*I-b.ofs_get_coef(i),i);
    return ofs;
}

/**
 * \fn template<typename T> Ofs<T> operator * (Ofs<T> const& a, T const& c)
 * \brief An operator. Makes the product c*a
 */
template <typename T>  Ofs<T> operator * (Ofs<T> const& a, T const& c)
{
    Ofs<T> cop(a);
    cop*=c;
    return cop;
}

/**
 * \fn template<typename T> Ofs<T> operator * (T const& c, Ofs<T> const& a)
 * \brief An operator. Makes the product c*a
 */
template <typename T>  Ofs<T> operator * (T const& c, Ofs<T> const& a)
{
    Ofs<T> cop(a);
    cop*=c;
    return cop;
}

/**
 * \fn template<typename T> Ofs<T> operator / (Ofs<T> const& a, T const& c)
 * \brief An operator. Makes the division a/c
 */
template <typename T>  Ofs<T> operator / (Ofs<T> const& a, T const& c)
{
    Ofs<T> cop(a);
    cop/=c;
    return cop;
}

/**
 * \fn template<typename T> Ofs<T> operator *  (Ofs<T> const& a, Ofs<T> const& b)
 * \brief An operator. Makes the product a*b.
 *
 *  Was supposed to be used instead of sprod/prod. In fact, makes the code a bit unclear
 *  and adds a hidden temporary variable to the implementation.
 *  As a consequence, use sprod/prod in priority.
 */
template <typename T>  Ofs<T>  operator * (Ofs<T> const& a, Ofs<T> const& b)
{
    int n, p, tO, psup, pinf;
    tO = max(a.get_order(), b.get_order());

    //Virgin ofs of same m_ofs_order
    Ofs<T> temp(tO);

    //Product
    for(n=-tO ; n<= tO; n++)
    {
        psup = min(n+a.get_order(),  a.get_order());
        pinf = max(n-b.get_order(), -b.get_order());
        for(p=pinf; p<= psup; p++)
        {
            //indix
            temp.add_coef(a.ofs_get_coef(p)*b.ofs_get_coef(n-p), n);
        }
    }

    return temp;
}


//----------------------------------------------------------------------------------------
// Functions (change of format)
//----------------------------------------------------------------------------------------
/**
 * \fn void inline ofs_double_to_complex(Ofsd& xr, Ofsc& xc)
 * \brief Copy from Ofsd to Ofsc.
 */
void inline ofs_double_to_complex(Ofsd const& xr, Ofsc& xc)
{
    int n_order_fourier = xr.get_order(); //m_ofs_order of the expansion
    if(n_order_fourier != xc.get_order()) //checking that the orders match
    {
        cout << "ofs_double_to_complex: orders do not match." << endl;
    }
    else for(int l = -n_order_fourier; l<=n_order_fourier; l++) xc.set_coef((cdouble) (xr.ofs_get_coef(l)+I*0.0), l); //copy from double to cdouble
}


//----------------------------------------------------------------------------------------
// Functions (Real and Imaginary part)
//----------------------------------------------------------------------------------------
/**
 * \fn void inline ofs_real_part(Ofsc& x, Ofsc& const xr)
 * \brief Takes the real part of an Ofsc object: xr = real(x).
 */
void inline ofs_real_part(Ofsc const& x, Ofsc& xr)
{
    int n_order_fourier = x.get_order();
    //xc = conj(xc)
    Ofsc xc(x);
    xc.conjugate();
    //Storing the real part
    for(int l = -n_order_fourier; l<=n_order_fourier; l++) xr.set_coef(0.5*(x.ofs_get_coef(l) + xc.ofs_get_coef(l)), l);
}

/**
 * \fn void inline ofs_real_part(Ofsc& x, Ofsc& const xr)
 * \brief Takes the real part of an Ofsc object: xi = imag(x).
 */
void inline ofs_imag_part(Ofsc const& x, Ofsc& xi)
{
    int n_order_fourier = x.get_order();
    //xc = conj(xc)
    Ofsc xc(x);
    xc.conjugate();
    //Storing the imag part
    for(int l = -n_order_fourier; l<=n_order_fourier; l++) xi.set_coef(-0.5*I*(x.ofs_get_coef(l) - xc.ofs_get_coef(l)), l);
}


//----------------------------------------------------------------------------------------
// Functions (evaluate)
//----------------------------------------------------------------------------------------
/**
 *  \brief  Evaluates the Ofs object at angle theta (theta = nt) and at a certain order eff_order.
 *
 *  The sum is given in the form: \f$ c_0 + \sum \limits_{k = 1}^{J} \cos(k\theta) (c_k + c_{-k}) + i\sin(k\theta) (c_k - c_{-k}) \f$
 *  in order to avoid any additionnal roundoff errors by using cos and sin functions as little as possible.
 */
template <typename T>  cdouble Ofs<T>::fevaluate(double cR[], double sR[], int eff_order)
{
    cdouble result = 0+0.0*I;
    //Order -k and +k are calculated together to avoid any additionnal roundoff errors
    for(int k=eff_order; k>=1; k--) result+= cR[k-1]*(m_ofs_coef[k+m_ofs_order] + m_ofs_coef[-k+m_ofs_order]) + I*sR[k-1]*(m_ofs_coef[k+m_ofs_order] - m_ofs_coef[-k+m_ofs_order]);
    //Order 0
    result+= m_ofs_coef[0+m_ofs_order];
    //Return result
    return result;
}

/**
 *  \brief  Evaluates the Ofs object at angle theta (theta = nt) and at a certain order eff_order.
 *
 *  The sum is given in the form: \f$ c_0 + \sum \limits_{k = 1}^{J} \cos(k\theta) (c_k + c_{-k}) + i\sin(k\theta) (c_k - c_{-k}) \f$
 *  in order to avoid any additionnal roundoff errors by using cos and sin functions as little as possible.
 */
template <typename T>  cdouble Ofs<T>::evaluate(double const& theta, int eff_order)
{
    cdouble result = 0+0.0*I;
    double cR[eff_order];
    double sR[eff_order];

    //Initialization of the cosinus/sinus arrays
    init_cR_sR(theta, cR, sR, eff_order);

    //Order -k and +k are calculated together to avoid any additionnal roundoff errors
    for(int k=eff_order; k>=1; k--) result+= cR[k-1]*(m_ofs_coef[k+m_ofs_order] + m_ofs_coef[-k+m_ofs_order]) + I*sR[k-1]*(m_ofs_coef[k+m_ofs_order] - m_ofs_coef[-k+m_ofs_order]);
    //Order 0
    result+= m_ofs_coef[0+m_ofs_order];//get_coef(0);


    return result;
    //Obsolete: for(int k=-m_ofs_order; k<=m_ofs_order; k++) result+= get_coef(k)*(cos((double)k*theta) + I*sin((double)k*theta));
}

/**
 *  \brief  Evaluates the Ofs object at angle theta (theta = nt).
 *
 *  The sum is given in the form: \f$ c_0 + \sum \limits_{k = 1}^{J} \cos(k\theta) (c_k + c_{-k}) + i\sin(k\theta) (c_k - c_{-k}) \f$
 *  in order to avoid any additionnal roundoff errors by using cos and sin functions as little as possible.
 */
template <typename T>  cdouble Ofs<T>::evaluate(double const& theta)
{
    cdouble result = 0+0.0*I;
    double cR[m_ofs_order];
    double sR[m_ofs_order];

    //Initialization of the cosinus/sinus arrays
    init_cR_sR(theta, cR, sR, m_ofs_order);

    //Order -k and +k are calculated together to avoid any additionnal roundoff errors
    for(int k=m_ofs_order; k>=1; k--) result+= cR[k-1]*(m_ofs_coef[k+m_ofs_order] + m_ofs_coef[-k+m_ofs_order]) + I*sR[k-1]*(m_ofs_coef[k+m_ofs_order] - m_ofs_coef[-k+m_ofs_order]);
    //Order 0
    result+= m_ofs_coef[0+m_ofs_order];//get_coef(0);

    return result;
    //Obsolete: for(int k=-m_ofs_order; k<=m_ofs_order; k++) result+= get_coef(k)*(cos((double)k*theta) + I*sin((double)k*theta));
}

/**
 *  \brief  Evaluates the derivatives of the Ofs object at angle theta, with pulsation n (theta = nt).
 *
 *  The sum is given in the form: \f$ c_0 + n \sum \limits_{k = 1}^{J} - k \sin(k\theta) (c_k + c_{-k}) + i k \cos(k\theta) (c_k - c_{-k}) \f$
 *  in order to avoid any additionnal roundoff errors by using cos and sin functions as little as possible.
 */
template <typename T>  cdouble Ofs<T>::evaluatedot(double const& theta, double const& n)
{
    cdouble result = 0+0.0*I;
    double cR[m_ofs_order];
    double sR[m_ofs_order];

    //Initialization of the cosinus/sinus arrays
    init_cR_sR(theta, cR, sR, m_ofs_order);

    //Order -k and +k are calculated together to avoid any additionnal roundoff errors
    for(int k=m_ofs_order; k>=1; k--) result+= -n*k*sR[k-1]*(m_ofs_coef[k+m_ofs_order] + m_ofs_coef[-k+m_ofs_order]) + I*n*k*cR[k-1]*(m_ofs_coef[k+m_ofs_order] - m_ofs_coef[-k+m_ofs_order]);

    return result;
    //Obsolete: for(int k=-m_ofs_order; k<=m_ofs_order; k++) result+= get_coef(k)*(cos((double)k*t) + I*sin((double)k*t));
}

/**
 *  \brief  Evaluates the derivatives of the Ofs object at angle theta, with pulsation n (theta = nt), at a certain order eff_order.
 *
 *  The sum is given in the form: \f$ c_0 + n \sum \limits_{k = 1}^{J} - k \sin(k\theta) (c_k + c_{-k}) + i k \cos(k\theta) (c_k - c_{-k}) \f$
 *  in order to avoid any additionnal roundoff errors by using cos and sin functions as little as possible.
 */
template <typename T>  cdouble Ofs<T>::evaluatedot(double const& theta, double const& n, int eff_order)
{
    cdouble result = 0+0.0*I;
    double cR[eff_order];
    double sR[eff_order];

    //Initialization of the cosinus/sinus arrays
    init_cR_sR(theta, cR, sR, eff_order);

    //Order -k and +k are calculated together to avoid any additionnal roundoff errors
    for(int k=eff_order; k>=1; k--) result+= -n*k*sR[k-1]*(m_ofs_coef[k+m_ofs_order] + m_ofs_coef[-k+m_ofs_order]) + I*n*k*cR[k-1]*(m_ofs_coef[k+m_ofs_order] - m_ofs_coef[-k+m_ofs_order]);

    return result;
    //Obsolete: for(int k=-m_ofs_order; k<=m_ofs_order; k++) result+= get_coef(k)*(cos((double)k*t) + I*sin((double)k*t));
}

/**
 *  \brief  Expected error on the product: \c this \f$ += a \times b \f$ at times t. Works only when a.m_ofs_order = b.m_ofs_order which is the default case.
 */
template <typename T>  cdouble sprod_expct_error(Ofs<T> const& a, Ofs<T> const& b, double const& t)
{
    cdouble result = 0.0+0.0*I;
    for(int n= a.get_order()+1 ; n<= 2*a.get_order(); n++)
    {
        for(int p=n-a.get_order(); p<= a.get_order(); p++)
        {
            //indix
            result += a.ofs_get_coef(p)*b.ofs_get_coef(n-p)*(cos(n*t) + I*sin(n*t));
            result += a.ofs_get_coef(-p)*b.ofs_get_coef(-n+p)*(cos(-n*t) + I*sin(-n*t));
        }
    }
    return result;
}

/**
 *  \brief  Expected error on the product: \c this \f$ += c*a \times b \f$ at times t. Works only when a.m_ofs_order = b.m_ofs_order which is the default case.
 */
template <typename T>  cdouble smprod_expct_error(Ofs<T> const& a, Ofs<T> const& b, T const& c, double const& t)
{
    cdouble result = 0.0+0.0*I;
    for(int n= a.get_order()+1 ; n<= 2*a.get_order(); n++)
    {
        for(int p=n-a.get_order(); p<= a.get_order(); p++)
        {
            //indix
            result += c*a.ofs_get_coef(p)*b.ofs_get_coef(n-p)*(cos(n*t) + I*sin(n*t));
            result += c*a.ofs_get_coef(-p)*b.ofs_get_coef(-n+p)*(cos(-n*t) + I*sin(-n*t));
        }
    }
    return result;
}

//Stream
//----------------------------------------------------------------------------------------
/**
 *  \brief  A stream operator
 */
template <typename T>  std::ostream& operator << (std::ostream& stream, Ofs<T> const& ofs)
{
    //stream << "Fourier series" << endl;
    //Order
    //stream << "Order : " << ofs.m_ofs_order << endl;
    //Coefficients
    for(int i = 0 ; i< 2*ofs.m_ofs_order + 1; i++)
    {
        stream << setw(3) << setiosflags(ios::right) << std::showpos << i-ofs.m_ofs_order << "   " <<  setiosflags(ios::scientific) << setprecision(15) << ofs.m_ofs_coef[i] << endl;

    }
    return stream;
}

/**
 *  \brief  A stream operator in the \c double \c complex case.
 */
template <>  inline std::ostream& operator << (std::ostream& stream, Ofs< cdouble > const& ofs)
{
    //Coefficients
    for(int i = 0 ; i< 2*ofs.m_ofs_order + 1; i++)
    {
        stream << setw(3) << setiosflags(ios::right) << std::showpos << i-ofs.m_ofs_order << "   " <<  setiosflags(ios::scientific) << setprecision(16) << creal(ofs.m_ofs_coef[i]) << "  " << cimag(ofs.m_ofs_coef[i]) << endl;

    }
    return stream;
}

/**
 *  \brief  Evaluates an upper bound for the \f$ L_1 \f$ norm of the current Ofs object.
 */
template<typename T> double Ofs<T>::l1norm()
{
//    double theta, temp;
//    double l1n = cabs(evaluate(0.0));
//    int N = 1000;
//    //Loop on a N+1 point grid from theta = 0 to 2pi
//    for(int idt = 1; idt <= N ; idt++)
//    {
//        theta = ((double)idt)/N*2*M_PI;
//        temp = cabs(evaluate(theta));
//        if(temp > l1n) l1n = temp;
//    }
//    return l1n;
    double l1n = 0.0;
    for(int i = 0 ; i< 2*m_ofs_order + 1; i++) l1n += cabs(m_ofs_coef[i]);
    return l1n;
}

/**
 *  \brief  Number of small divisors under a certain value
*/
template<typename T> int Ofs<T>::nsd(int odmax, double sdmax)
{
    int res = 0;
    for(int i = -odmax; i <= odmax; i++) if(cabs(this->ofs_get_coef(i)) < sdmax)
        {
            //cout << i << " sd = " << creal(this->ofs_get_coef(i)) << " " << cimag(this->ofs_get_coef(i)) << endl;
            res++;
        }
    return res;
}

/**
 *  \brief  A stream operator. Print only the autonomous term (term of order zero).
 */
template<typename T> void Ofs<T>::fprint_0(ofstream& stream)
{
    stream << setw(3) << setiosflags(ios::right) << std::showpos << 0 << "   " <<  setiosflags(ios::scientific) << setprecision(15) << creal(this->ofs_get_coef(0)) << "  " << cimag(this->ofs_get_coef(0)) << endl;
}


//----------------------------------------------------------------------------------------
// Reading an OFS from a text file
//----------------------------------------------------------------------------------------
/**
 * \fn void inline read_ofs_txt(Ofsc& xFFT, string filename, int fftN)
 * \brief Reading an Ofsc object from a text file.
 */
void inline read_ofs_txt(Ofsc& xFFT, string filename)
{
    //Init
    ifstream readStream;
    double ct, cr, ci;
    int fftN = xFFT.get_order();

    //Reading
    readStream.open((filename+".txt").c_str());
    for(int i = -fftN; i<=fftN; i++)
    {
        readStream >> ct;  //current m_ofs_order
        readStream >> cr;  //real part
        readStream >> ci;  //imag part
        xFFT.set_coef(cr+I*ci, i);
    }
    readStream.close();
}

//----------------------------------------------------------------------------------------
// Array cos/sin computation
//----------------------------------------------------------------------------------------
/**
 *  \brief Initializes the arrays cR[] and sR[], containing the numbers cos(t), cos(nt), ... cos(ofs_order*n*t) and sin(t), sin(nt), ... sin(ofs_order*n*t), respectively.
 **/
void inline init_cR_sR(double t, double cR[], double sR[], int ofs_order)
{
    cR[0] = cos(t);
    sR[0] = sin(t);
    for(int i = 1; i< ofs_order; i++)
    {
        cR[i] =  cR[i-1]*cR[0] - sR[i-1]*sR[0];
        sR[i] =  sR[i-1]*cR[0] + cR[i-1]*sR[0];
    }
}

//----------------------------------------------------------------------------------------
// Single storage
//----------------------------------------------------------------------------------------
/**
 *  \brief Single storage of a QBTBP Ofs object in txt files + its cosinus/sinus version.
 */
inline void ofs_sst(Ofsc &xOFS, string filename, int flag, string suffix)
{
    int n_order_fourier = xOFS.get_order();
    ofstream curentStream;

    //Storage in txt file
    curentStream.open((filename+suffix+".txt").c_str());
    curentStream <<  xOFS << endl;
    curentStream.close();

    curentStream.open((filename+"c"+suffix+".txt").c_str());
    if(flag) //even case
    {
        //Cosinus expansion version
        curentStream << 0 << " " << creal(xOFS.ofs_get_coef(0)) << endl;
        for(int l = 1; l<=n_order_fourier; l++) curentStream << l << " " << creal(xOFS.ofs_get_coef(-l) + xOFS.ofs_get_coef(l))  <<  endl;

    }
    else //odd case
    {
        //Sinus expansion version
        for(int l = 0; l<=n_order_fourier; l++)
            curentStream << l << " " <<  cimag(xOFS.ofs_get_coef(-l) - xOFS.ofs_get_coef(l)) << endl;
        curentStream.close();
    }
    curentStream.close();
}
