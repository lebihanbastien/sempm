#include "manip.h"

/**
 * \file manip.cpp
 * \brief Basic operations for manipulation of polynomials
 *       (exponents in reverse lexicographic order, number of coefficient in
 *       homogeneous polynomials...).
 *       Note that the algebraic operations on series (sum, product) are not done here
 *       but in ofts.h/ofts.tpp.
 * \author BLB using code by Angel Jorba.
 *
 *      Based on Jorba 1999 (http://www.maia.ub.es/~angel/soft.html).
 */

using namespace std;

int  Manip::m_max_deg=0;  //static maximum degree
int  Manip::m_nvar=0;     //static number of variables
long int **Manip::m_num_mon;    //contains m_num_mon(i,j) for i=2..m_nvar.


//----------------------------------------------------------------------------------------
//Class routines
//----------------------------------------------------------------------------------------
/**
 *   \brief Initializes the table m_num_mon(i,j), which contains the number of monomials
 *          of degree j with i variables.
 *
 *   parameters:
 *   t_max_deg: maximum degree we are going to work with. it can not be greater than 63.
 *   t_nvar: number of variables
 *
 *   returned value: number of kbytes allocated by the internal tables.
 *   Based on a routine by Angel Jorba, 1999.
 **/
int  Manip::init(int t_nvar, int t_max_deg)
{
    int i,j,l;
    unsigned long int mem; /* mem: to count the amount of memory used */

    if (m_max_deg != 0) /* this means that Manip is already initialized */
    {
        if (t_max_deg <= m_max_deg)
        {
            puts("**************************************************");
            printf("%s warning message:\n"," init");
            printf("%s is already initialized to degree %d\n"," init",m_max_deg);
            printf("now you want to initialize it again to the same degree.\n");
            printf("action taken: this call to %s is simply ignored.\n"," init");
            puts("**************************************************");
            return(0);
        }
        printf("%s error message:\n"," init");
        printf("%s is already initialized to degree %d\n"," init",m_max_deg);
        printf("and now you want to initialize it to degree %d.\n",t_max_deg);
        printf("you must call routine %s first.\n"," free");
        printf("action taken: program aborted\n");
        exit(1); /* this will stop the program */
    }
    m_max_deg=t_max_deg;
    m_nvar=t_nvar;
    mem=(t_nvar-1)*sizeof(int*); /* this is to count the number of kb. allocated */
    m_num_mon= (long  int**)malloc((t_nvar)*sizeof(long int*));
    if (m_num_mon == NULL)
    {
        puts("Manip error. no memory (1).");
        exit(1);
    }
    m_num_mon -= sizeof(long  int*);
    for (i=1; i<=t_nvar; i++)
    {
        m_num_mon[i]=(long int*)malloc((t_max_deg+1)*sizeof(long  int));
        if (m_num_mon[i] == NULL)
        {
            puts("Manip error. no memory (2).");
            exit(1);
        }
        mem += (t_max_deg+1)*sizeof(int);
    }
    for (j=0; j<=t_max_deg; j++) m_num_mon[1][j]=1;

    if(t_nvar > 1)
    {
        for (j=0; j<=t_max_deg; j++) m_num_mon[2][j]=j+1;
        for (i=3; i<=t_nvar; i++)
        {
            for (j=0; j<=t_max_deg; j++)
            {
                m_num_mon[i][j]=0;
                for (l=0; l<=j; l++) m_num_mon[i][j] += m_num_mon[i-1][l];
            }
        }
    }
    mem += (t_max_deg+1)*sizeof(int*);
    mem /= 1024;

    return(mem);
}

/**
 *  \brief Frees the space allocated by Manip::init.
 **/
void  Manip::free()
{
    int i;
    if (m_max_deg == 0)
    {
        puts("**************************************************");
        printf("%s warning message:\n"," free");
        printf("no memory to free\n");
        printf("action taken: this call is simply ignored.\n");
        puts("**************************************************");

        return;
    }
    cout << "m_nvar: " << m_nvar << endl;
    for (i=1; i<=m_nvar; i++)
    {
        cout << "i: " << i << endl;
        delete m_num_mon[i];
    }
    m_num_mon += sizeof(long  int*);
    delete m_num_mon;
    m_max_deg=0;
    m_nvar=0;
    return;
}

/**
 *   \brief Returns the number of monomials of degree t_max_deg with t_nvar variables,
 *          making use of the table Manip::m_num_mon.
 *
 *  parameters:
 *  t_nvar: number of variables
 *  t_max_deg: order we are interested in (input).
 *  returned value: number of monomials of order no.
 **/
long int Manip::nmon(int t_nvar, int t_max_deg)
{
    if (t_max_deg > m_max_deg)
    {
        printf("nmon: error, the requested degree %i is greater than m_max_deg=%i.\n", t_max_deg, m_max_deg );
        exit(1);
    }
    if (t_nvar > m_nvar)
    {
        puts("nmon: error, the requested numver of variable is greater than m_nvar.");
        exit(1);
    }

    if(t_nvar <= 1)
        return 1;
    else
    return(m_num_mon[t_nvar][t_max_deg]);

}

/**
 *  \brief given a multiindex k, this routine computes the next one
 *  according to the (reverse) lexicographic order.
 *
 *  parameters:
 *  t_array: array of t_nvar components containing the multiindex. It is overwritten on exit
 *  (input and output).
 **/
void  Manip::prxkt(int t_array[], int t_nvar)
{
    if(t_nvar == 0)
    {
        t_array[0]++;
        return;
    }
    else
    {
        int i;
        if (t_array[0] != 0)
        {
            t_array[0]--;
            t_array[1]++;
            return;
        }
        for(i=1; i<t_nvar-1; i++)
        {
            if (t_array[i] != 0)
            {
                t_array[0]=t_array[i]-1;
                t_array[i]=0;
                t_array[i+1]++;
                return;
            }
        }
    }
    puts("prxkt error 1.");
    exit(1);
}

/**
 *   \brief Returns the number of product operation in a taylor series product
 **/
long int pdk(int t_nvar, int t_order)
{
    long int result = 0;
    int k,i;
    for(k=0; k<=t_order; k++)
    {
        for(i=0; i<=k; i++) result+= Manip::nmon(t_nvar,i)* Manip::nmon(t_nvar,k-i); //result+=binomial(t_nvar+i-1, t_nvar-1)*binomial(t_nvar+k-i-1, t_nvar-1); //
    }
    return result;
}


//----------------------------------------------------------------------------------------
//Binomial coefficients. These routines do not make use of Manip::m_num_mon.
//----------------------------------------------------------------------------------------
/**
 * \brief Computes the binomial coefficient (x y).
 **/
unsigned long gcd_ui(unsigned long x, unsigned long y)
{
    unsigned long t;
    if (y < x)
    {
        t = x;
        x = y;
        y = t;
    }
    while (y > 0)
    {
        t = y;
        y = x % y;
        x = t;  /* y1 <- x0 % y0 ; x1 <- y0 */
    }
    return x;
}

/**
 * \brief Computes the binomial coefficient (n k).
 **/
unsigned long binomial(unsigned long n, unsigned long k)
{
    unsigned long d, g, r = 1;
    if (k == 0) return 1;
    if (k == 1) return n;
    if (k >= n) return (k == n);
    if (k > n/2) k = n-k;
    for (d = 1; d <= k; d++)
    {
        if (r >= ULONG_MAX/n)    /* Possible overflow */
        {
            unsigned long t_max_deg, dr;  /* reduced numerator / denominator */
            g = gcd_ui(n, d);
            t_max_deg = n/g;
            dr = d/g;
            g = gcd_ui(r, dr);
            r = r/g;
            dr = dr/g;
            if (r >= ULONG_MAX/t_max_deg) return 0;  /* Unavoidable overflow */
            r *= t_max_deg;
            r /= dr;
            n--;
        }
        else
        {
            r *= n--;
            r /= d;
        }
    }
    return r;
}


