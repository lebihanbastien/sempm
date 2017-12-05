/**
  @file     gnuplot_i.c
  @author   N. Devillard, Peter X, BLB
  @version  Custom
  @brief    C interface to gnuplot.

  gnuplot is a freely available, command-driven graphical display tool for
  Unix. It compiles and works quite well on a number of Unix flavours as
  well as other operating systems. The following module enables sending
  display requests to gnuplot through simple C calls.

*/


/*----------------------------------------------------------------------------------------
                                Includes
 ----------------------------------------------------------------------------------------*/
#include "gnuplot_i.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <assert.h>

#ifdef WIN32
#include <io.h>
#endif // #ifdef WIN32


//----------------------------------------------------------------------------------------
//                              Defines
//----------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------
//                          Prototype Functions
//----------------------------------------------------------------------------------------
/**
 * Creates a temporary file name for writing
 *
 * @author Peter (12/9/2011)
 *
 * @param handle
 *
 * @return char const * Pointer to file name string.
 */
char const * gnuplot_tmpfile(gnuplot_ctrl * handle);

/**
 * Plot a temporary file.
 *
 * @author Peter (12/9/2011)
 *
 * @param handle
 * @param tmp_filename
 * @param title
 */
void gnuplot_plot_atmpfile(gnuplot_ctrl * handle,
                           char const* tmp_filename,
                           char const* title,
                           char const* ls,
                           char const* lt,
                           char const* lw,
                           int lc);

/**
 * Plot3D a temporary file.
 *
 * @author Bastien (2015)
 *
 * @param handle
 * @param tmp_filename
 * @param title
 */
void gnuplot_plot_atmpfile_3d(gnuplot_ctrl * handle,
                              char const* tmp_filename,
                              char const* title,
                              char const* ls,
                              char const* lt,
                              char const* lw,
                              int lc);

//----------------------------------------------------------------------------------------
//                       Initialization and closing
//----------------------------------------------------------------------------------------
/**
  @brief    Opens up a gnuplot session, ready to receive commands.
  @return   Newly allocated gnuplot control structure.

  This opens up a new gnuplot session, ready for input. The struct
  controlling a gnuplot session should remain opaque and only be
  accessed through the provided functions.

  The session must be closed using gnuplot_close().
 */
gnuplot_ctrl * gnuplot_init(void)
{
    gnuplot_ctrl *  handle ;
    int i;

#ifndef WIN32
    if (getenv("DISPLAY") == NULL) {
        fprintf(stderr, "cannot find DISPLAY variable: is it set?\n") ;
    }
#endif // #ifndef WIN32


    /*
     * Structure initialization:
     */
    handle = (gnuplot_ctrl*)malloc(sizeof(gnuplot_ctrl)) ;
    handle->nplots = 0 ;
    gnuplot_setstyle(handle, "points") ;
    handle->ntmp = 0 ;
    gnuplot_setcolor(handle, 0);


    handle->gnucmd = popen("gnuplot", "w") ;
    if (handle->gnucmd == NULL) {
        fprintf(stderr, "error starting gnuplot, is gnuplot or gnuplot.exe in your path?\n") ;
        free(handle) ;
        return NULL ;
    }

    for (i=0;i<GP_MAX_TMP_FILES; i++)
    {
        handle->tmp_filename_tbl[i] = NULL;
    }
    return handle;
}

/**
  @brief    Closes a gnuplot session previously opened by gnuplot_init()
  @param    handle Gnuplot session control handle.
  @return   void

  Kills the child PID and deletes all opened temporary files.
  It is mandatory to call this function to close the handle, otherwise
  temporary files are not cleaned and child process might survive.

 */
void gnuplot_close(gnuplot_ctrl * handle)
{
    int     i ;

    if (pclose(handle->gnucmd) == -1) {
        fprintf(stderr, "problem closing communication to gnuplot\n") ;
        return ;
    }
    if (handle->ntmp) {
        for (i=0 ; i<handle->ntmp ; i++) {
            remove(handle->tmp_filename_tbl[i]) ;
            free(handle->tmp_filename_tbl[i]);
            handle->tmp_filename_tbl[i] = NULL;

        }
    }
    free(handle) ;
    return ;
}

/**
  @brief    Resets a gnuplot session (next plot will erase previous ones).
  @param    h Gnuplot session control handle.
  @return   void

  Resets a gnuplot session, i.e. the next plot will erase all previous
  ones.
 */
void gnuplot_resetplot(gnuplot_ctrl * h)
{
    int     i ;
    if (h->ntmp) {
        for (i=0 ; i<h->ntmp ; i++) {
            remove(h->tmp_filename_tbl[i]) ;
            free(h->tmp_filename_tbl[i]);
            h->tmp_filename_tbl[i] = NULL;

        }
    }
    h->ntmp = 0 ;
    h->nplots = 0 ;
    return ;
}

//----------------------------------------------------------------------------------------
//                       Gnuplot commands
//----------------------------------------------------------------------------------------
/**
  @brief    Sends a command to an active gnuplot session.
  @param    handle Gnuplot session control handle
  @param    cmd    Command to send, same as a printf statement.

  This sends a string to an active gnuplot session, to be executed.
  There is strictly no way to know if the command has been
  successfully executed or not.
  The command syntax is the same as printf.

  Examples:

  @code
  gnuplot_cmd(g, "plot %d*x", 23.0);
  gnuplot_cmd(g, "plot %.18e * cos(%.18e * x)", 32.0, -3.0);
  @endcode

  Since the communication to the gnuplot process is run through
  a standard Unix pipe, it is only unidirectional. This means that
  it is not possible for this interface to query an error status
  back from gnuplot.
 */
void gnuplot_cmd(gnuplot_ctrl *  handle, char const *  cmd, ...)
{
    va_list ap ;

    va_start(ap, cmd);
    vfprintf(handle->gnucmd, cmd, ap);
    va_end(ap);

    fputs("\n", handle->gnucmd) ;
    fflush(handle->gnucmd) ;
    return ;
}

/**
  @brief    Change the plotting style of a gnuplot session.
  @param    h Gnuplot session control handle
  @param    plot_style Plotting-style to use (character string)
  @return   void

  The provided plotting style is a character string. It must be one of
  the following:

  - lines
  - points
  - linespoints
  - impulses
  - dots
  - steps
  - errorbars
  - boxes
  - boxeserrorbars
 */
void gnuplot_setstyle(gnuplot_ctrl * h, char * plot_style)
{
    if (strcmp(plot_style, "lines") &&
        strcmp(plot_style, "points") &&
        strcmp(plot_style, "linespoints") &&
        strcmp(plot_style, "impulses") &&
        strcmp(plot_style, "dots") &&
        strcmp(plot_style, "steps") &&
        strcmp(plot_style, "errorbars") &&
        strcmp(plot_style, "boxes") &&
        strcmp(plot_style, "boxerrorbars")) {
        fprintf(stderr, "warning: unknown requested style: using points\n") ;
        strcpy(h->pstyle, "points") ;
    } else {
        strcpy(h->pstyle, plot_style) ;
    }
    return ;
}

/**
  @brief    Change the plotting color of a gnuplot session.
  @param    h Gnuplot session control handle
  @param    plot_color Plotting-color to use (int)
  @return   void
 */
void gnuplot_setcolor(gnuplot_ctrl * h, int plot_style)
{
    switch(plot_style)
    {
        case 1:
        strcpy(h->pcolor, "dark-violet") ;
        break;
        case 2:
        strcpy(h->pcolor, "#009e73") ;
        break;
        case 3:
        strcpy(h->pcolor, "#56b4e9") ;
        break;
        case 4:
        strcpy(h->pcolor, "#e69f00") ;
        break;
        case 5:
        strcpy(h->pcolor, "#f0e442") ;
        break;
        case 6:
        strcpy(h->pcolor, "#0072b2") ;
        break;
        case 7:
        strcpy(h->pcolor, "#e51e10") ;
        break;
        case 8:
        strcpy(h->pcolor, "black") ;
        break;
        case 9:
        strcpy(h->pcolor, "gray50") ;
        break;
        default:
        strcpy(h->pcolor, "default") ;
        break;
    }
    return ;
}

/**
  @brief    Sets the x label of a gnuplot session.
  @param    h Gnuplot session control handle.
  @param    label Character string to use for X label.
  @return   void

  Sets the x label for a gnuplot session.
 */
void gnuplot_set_xlabel(gnuplot_ctrl * h, char * label)
{
    gnuplot_cmd(h, "set xlabel \"%s\"", label) ;
}

/**
  @brief    Sets the y label of a gnuplot session.
  @param    h Gnuplot session control handle.
  @param    label Character string to use for Y label.
  @return   void

  Sets the y label for a gnuplot session.
 */
void gnuplot_set_ylabel(gnuplot_ctrl * h, char * label)
{
    gnuplot_cmd(h, "set ylabel \"%s\"", label) ;
}

/**
  @brief    Sets the z label of a gnuplot session.
  @param    h Gnuplot session control handle.
  @param    label Character string to use for Z label.
  @return   void

  Sets the z label for a gnuplot session.
 */
void gnuplot_set_zlabel(gnuplot_ctrl * h, char * label)
{
    gnuplot_cmd(h, "set zlabel \"%s\"", label) ;
}


//----------------------------------------------------------------------------------------
//                       Plotting
//----------------------------------------------------------------------------------------
/**
  @brief    Plots a 2d graph from a list of doubles.
  @param    handle  Gnuplot session control handle.
  @param    d       Array of doubles.
  @param    n       Number of values in the passed array.
  @param    title   Title of the plot.
  @return   void

  Plots out a 2d graph from a list of doubles. The x-coordinate is the
  index of the double in the list, the y coordinate is the double in
  the list.

  Example:

  @code
    gnuplot_ctrl    *h ;
    double          d[50] ;
    int             i ;

    h = gnuplot_init() ;
    for (i=0 ; i<50 ; i++) {
        d[i] = (double)(i*i) ;
    }
    gnuplot_plot_x(h, d, 50, "parabola") ;
    sleep(2) ;
    gnuplot_close(h) ;
  @endcode
 */
void gnuplot_plot_x(
    gnuplot_ctrl    *   handle,
    double          *   d,
    int                 n,
    char            *   title
)
{
    int     i ;
    FILE*   tmpfd ;
    char const * tmpfname;

    if (handle==NULL || d==NULL || (n<1)) return ;

    /* Open temporary file for output   */
    tmpfname = gnuplot_tmpfile(handle);
    tmpfd = fopen(tmpfname, "w");

    if (tmpfd == NULL) {
        fprintf(stderr,"cannot create temporary file: exiting plot") ;
        return ;
    }

    /* Write data to this file  */
    for (i=0 ; i<n ; i++) {
      fprintf(tmpfd, "%.18e\n", d[i]);
    }
    fclose(tmpfd) ;

    gnuplot_plot_atmpfile(handle,tmpfname,title, "lines", "1", "2", 1);
    return ;
}


/**
  @brief    Plot a 2d graph from a list of points.
  @param    handle      Gnuplot session control handle.
  @param    x           Pointer to a list of x coordinates.
  @param    y           Pointer to a list of y coordinates.
  @param    n           Number of doubles in x (assumed the same as in y).
  @param    title       Title of the plot.
  @return   void

  Plots out a 2d graph from a list of points. Provide points through a list
  of x and a list of y coordinates. Both provided arrays are assumed to
  contain the same number of values.

  @code
    gnuplot_ctrl    *h ;
    double          x[50] ;
    double          y[50] ;
    int             i ;

    h = gnuplot_init() ;
    for (i=0 ; i<50 ; i++) {
        x[i] = (double)(i)/10.0 ;
        y[i] = x[i] * x[i] ;
    }
    gnuplot_plot_xy(h, x, y, 50, "parabola") ;
    sleep(2) ;
    gnuplot_close(h) ;
  @endcode
 */
void gnuplot_plot_xy(
    gnuplot_ctrl    *   handle,
    double          *   x,
    double          *   y,
    int                 n,
    char const           *   title,
    char const      * ls,
    char const      * lt,
    char const      * lw,
    int lc)
{
    int     i ;
    FILE*   tmpfd ;
    char const * tmpfname;

    if (handle==NULL || x==NULL || y==NULL || (n<1)) return ;

    /* Open temporary file for output   */
    tmpfname = gnuplot_tmpfile(handle);
    tmpfd = fopen(tmpfname, "w");

    if (tmpfd == NULL) {
        fprintf(stderr,"cannot create temporary file: exiting plot") ;
        return ;
    }

    /* Write data to this file  */
    for (i=0 ; i<n; i++) {
        fprintf(tmpfd, "%.18e %.18e\n", x[i], y[i]) ;
    }
    fclose(tmpfd) ;

    gnuplot_plot_atmpfile(handle,tmpfname,title, ls, lt, lw, lc);
    return ;
}

/**
  @brief    Plot a 3d graph from a list of points.
  @param    handle      Gnuplot session control handle.
  @param    x           Pointer to a list of x coordinates.
  @param    y           Pointer to a list of y coordinates.
  @param    z           Pointer to a list of z coordinates.
  @param    n           Number of doubles in x (assumed the same as in y and z).
  @param    title       Title of the plot.
  @return   void

  Plots out a 3d graph from a list of points. Provide points through a list
  of x, y, and z coordinates. Provided arrays are assumed to
  contain the same number of values.
 */
void gnuplot_plot_xyz(
    gnuplot_ctrl    * handle,
    double          * x,
    double          * y,
    double          * z,
    int               n,
    char const      * title,
    char const      * ls,
    char const      * lt,
    char const      * lw,
    int lc)
{
    int     i ;
    FILE*   tmpfd ;
    char const * tmpfname;

    if (handle==NULL || x==NULL || y==NULL || (n<1)) return ;

    /* Open temporary file for output   */
    tmpfname = gnuplot_tmpfile(handle);
    tmpfd = fopen(tmpfname, "w");

    if (tmpfd == NULL) {
        fprintf(stderr,"cannot create temporary file: exiting plot") ;
        return ;
    }

    /* Write data to this file  */
    for (i=0 ; i<n; i++) {
        fprintf(tmpfd, "%.18e %.18e %.18e\n", x[i], y[i], z[i]) ;
    }
    fclose(tmpfd) ;

    gnuplot_plot_atmpfile_3d(handle,tmpfname,title, ls, lt, lw, lc);
    return ;
}

/**
  @brief    Store a 2d graph in a txt file from a list of points.
  @param    handle      Gnuplot session control handle.
  @param    x           Pointer to a list of x coordinates.
  @param    y           Pointer to a list of y coordinates.
  @param    n           Number of doubles in x (assumed the same as in y).
  @param    tmpfname    The name of the storing file
  @return   void

  Store a 2d graph in a txt file from a list of points. Provide points through a list
  of x and a list of y coordinates. Both provided arrays are assumed to
  contain the same number of values.
 */
void gnuplot_fplot_xy(
    double          * x,
    double          * y,
    int               n,
    char const      * tmpfname)
{
    int     i ;
    FILE*   tmpfd ;
    if (x==NULL || y==NULL || (n<1)) return ;

    /* Open file for output   */
    tmpfd = fopen(tmpfname, "w");

    if (tmpfd == NULL) {
        fprintf(stderr,"cannot create file: exiting plot") ;
        return ;
    }

    /* Write data to this file  */
    for (i=0 ; i<n; i++) {
        fprintf(tmpfd, "%.18e %.18e \n", x[i], y[i]) ;
    }
    fclose(tmpfd) ;

    return ;
}

/**
  @brief    Store a 3d graph in a txt file from a list of points.
  @param    handle      Gnuplot session control handle.
  @param    x           Pointer to a list of x coordinates.
  @param    y           Pointer to a list of y coordinates.
  @param    z           Pointer to a list of z coordinates.
  @param    n           Number of doubles in x (assumed the same as in y and z).
  @param    tmpfname    The name of the storing file
  @return   void

  Store a 3d graph in a txt file from a list of points. Provide points through a list
  of x, y, and z coordinates. Provided arrays are assumed to
  contain the same number of values.
 */
void gnuplot_fplot_xyz(
    double          * x,
    double          * y,
    double          * z,
    int               n,
    char const      * tmpfname)
{
    int     i ;
    FILE*   tmpfd ;
    if (x==NULL || y==NULL || (n<1)) return ;

    /* Open file for output   */
    tmpfd = fopen(tmpfname, "w");

    if (tmpfd == NULL) {
        fprintf(stderr,"cannot create file: exiting plot") ;
        return ;
    }

    /* Write data to this file  */
    for (i=0 ; i<n; i++) {
        fprintf(tmpfd, "%.18e, %.18e, %.18e\n", x[i], y[i], z[i]) ;
    }
    fclose(tmpfd) ;

    return ;
}

/**
  @brief    Store a 3d graph and the corresponding time grid in a txt file from a list of points.
  @param    handle      Gnuplot session control handle.
  @param    t           Pointer to a list of time coordinates.
  @param    x           Pointer to a list of x coordinates.
  @param    y           Pointer to a list of y coordinates.
  @param    z           Pointer to a list of z coordinates.
  @param    n           Number of doubles in x (assumed the same as in y and z).
  @param    tmpfname    The name of the storing file
  @return   void

  Store a 3d graph and the corresponding time grid in a txt file from a list of points.
  Provide points through a list of t, x, y, and z coordinates. Provided arrays are assumed to
  contain the same number of values.
 */
void gnuplot_fplot_txyz(
    double          * t,
    double          * x,
    double          * y,
    double          * z,
    int               n,
    char const      * tmpfname)
{
    int     i ;
    FILE*   tmpfd ;
    if (x==NULL || y==NULL || (n<1)) return ;

    /* Open file for output   */
    tmpfd = fopen(tmpfname, "w");

    if (tmpfd == NULL) {
        fprintf(stderr,"cannot create file: exiting plot") ;
        return ;
    }

    /* Write data to this file  */
    for (i=0 ; i<n; i++) {
        fprintf(tmpfd, "%.18e %.18e %.18e %.18e\n", t[i], x[i], y[i], z[i]) ;
    }
    fclose(tmpfd) ;

    return ;
}


/**
 *  \brief Generates a temporary txt file for gnuplot use.
 **/
char const * gnuplot_tmpfile(gnuplot_ctrl * handle)
{
    static char const * tmp_filename_template = "tmp/gnuplot_tmpdatafile_XXXXXX";
    char *              tmp_filename = NULL;
    int                 tmp_filelen = strlen(tmp_filename_template);

#ifndef WIN32
    int                 unx_fd;
#endif // #ifndef WIN32

    assert(handle->tmp_filename_tbl[handle->ntmp] == NULL);

    /* Open one more temporary file? */
    if (handle->ntmp == GP_MAX_TMP_FILES - 1) {
        fprintf(stderr,
                "maximum # of temporary files reached (%d): cannot open more",
                GP_MAX_TMP_FILES) ;
        return NULL;
    }

    tmp_filename = (char*) malloc(tmp_filelen+1);
    if (tmp_filename == NULL)
    {
        return NULL;
    }
    strcpy(tmp_filename, tmp_filename_template);

#ifdef WIN32
    if (_mktemp(tmp_filename) == NULL)
    {
        return NULL;
    }
#else // #ifdef WIN32
    unx_fd = mkstemp(tmp_filename);
    if (unx_fd == -1)
    {
        return NULL;
    }
    //close(unx_fd);

#endif // #ifdef WIN32

    handle->tmp_filename_tbl[handle->ntmp] = tmp_filename;
    handle->ntmp ++;
    return tmp_filename;
}


/**
 *   \brief Generates a gnuplot command in order to plot the 2D data stored in the
 *          temporary file tmp_filename.
 **/
void gnuplot_plot_atmpfile(gnuplot_ctrl * handle, char const* tmp_filename,
                           char const* title, char const* ls, char const* lt,
                           char const* lw, int lc)
{
    char const *    cmd    = (handle->nplots > 0) ? "replot" : "plot";
    title                  = (title == NULL)      ? "(none)" : title;

    switch(lc % 9)
    {
        case 1:
        strcpy(handle->pcolor, "dark-violet") ;
        break;
        case 2:
        strcpy(handle->pcolor, "#009e73") ;
        break;
        case 3:
        strcpy(handle->pcolor, "#56b4e9") ;
        break;
        case 4:
        strcpy(handle->pcolor, "#e69f00") ;
        break;
        case 5:
        strcpy(handle->pcolor, "#f0e442") ;
        break;
        case 6:
        strcpy(handle->pcolor, "#0072b2") ;
        break;
        case 7:
        strcpy(handle->pcolor, "#e51e10") ;
        break;
        case 8:
        strcpy(handle->pcolor, "black") ;
        break;
        case 9:
        strcpy(handle->pcolor, "gray50") ;
        break;
        default:
        strcpy(handle->pcolor, "dark-violet") ;
        break;
    }

    gnuplot_cmd(handle, "%s \"%s\" title \"%s\" with %s lt \"%s\" lw %s lc rgb \"%s\"", cmd, tmp_filename,title, ls, lt, lw, handle->pcolor);
    handle->nplots++ ;
    return ;
}

/**
 *   \brief Same as gnuplot_plot_atmpfile for 3D plots.
 **/
void gnuplot_plot_atmpfile_3d(gnuplot_ctrl * handle,
                              char const* tmp_filename,
                              char const* title,
                              char const* ls,
                              char const* lt,
                              char const* lw,
                              int lc)
{
    char const *    cmd    = (handle->nplots > 0) ? "replot" : "splot";
    title                  = (title == NULL)      ? "(none)" : title;

    switch(lc % 9)
    {
        case 1:
        strcpy(handle->pcolor, "dark-violet") ;
        break;
        case 2:
        strcpy(handle->pcolor, "#009e73") ;
        break;
        case 3:
        strcpy(handle->pcolor, "#56b4e9") ;
        break;
        case 4:
        strcpy(handle->pcolor, "#e69f00") ;
        break;
        case 5:
        strcpy(handle->pcolor, "#f0e442") ;
        break;
        case 6:
        strcpy(handle->pcolor, "#0072b2") ;
        break;
        case 7:
        strcpy(handle->pcolor, "#e51e10") ;
        break;
        case 8:
        strcpy(handle->pcolor, "black") ;
        break;
        case 9:
        strcpy(handle->pcolor, "gray50") ;
        break;
        default:
        strcpy(handle->pcolor, "dark-violet") ;
        break;
    }
    gnuplot_cmd(handle, "set view 60,75");
    gnuplot_cmd(handle, "%s \"%s\" title \"%s\" with %s lt \"%s\" lw %s lc rgb \"%s\"", cmd, tmp_filename, title, ls, lt, lw, handle->pcolor);
    handle->nplots++ ;
    return ;
}

