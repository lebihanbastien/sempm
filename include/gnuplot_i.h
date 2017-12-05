#ifndef _GNUPLOT_PIPES_H_
#define _GNUPLOT_PIPES_H_

/**
  @file     gnuplot_i.h
  @author   N. Devillard, Peter X, BLB
  @version  Custom
  @brief    C interface to gnuplot.

  gnuplot is a freely available, command-driven graphical display tool for
  Unix. It compiles and works quite well on a number of Unix flavours as
  well as other operating systems. The following module enables sending
  display requests to gnuplot through simple C calls.

*/

#include <stdio.h>

/** Maximal number of simultaneous temporary files */
#define GP_MAX_TMP_FILES    256


/**
  \typedef  gnuplot_ctrl
  \brief    gnuplot session handle (opaque type).

  This structure holds all necessary information to talk to a gnuplot
  session. It is built and returned by gnuplot_init() and later used
  by all functions in this module to communicate with the session, then
  meant to be closed by gnuplot_close().

  This structure is meant to remain opaque, you normally do not need
  to know what is contained in there.
 */
typedef struct gnuplot_ctrl {
    /** Pipe to gnuplot process */
    FILE    * gnucmd ;

    /** Number of currently active plots */
    int       nplots ;
    /** Current plotting style */
    char      pstyle[32] ;
    /** Current plotting color */
    char      pcolor[32] ;
    /** Pointer to table of names of temporary files */
    char*      tmp_filename_tbl[GP_MAX_TMP_FILES] ;
    /** Number of temporary files */
    int       ntmp ;
}gnuplot_ctrl;

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
gnuplot_ctrl * gnuplot_init(void);

/**
  @brief    Closes a gnuplot session previously opened by gnuplot_init()
  @param    handle Gnuplot session control handle.
  @return   void

  Kills the child PID and deletes all opened temporary files.
  It is mandatory to call this function to close the handle, otherwise
  temporary files are not cleaned and child process might survive.

 */
void gnuplot_close(gnuplot_ctrl * handle);

/**
  @brief    Resets a gnuplot session (next plot will erase previous ones).
  @param    h Gnuplot session control handle.
  @return   void

  Resets a gnuplot session, i.e. the next plot will erase all previous
  ones.
 */
void gnuplot_resetplot(gnuplot_ctrl * h);

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
void gnuplot_cmd(gnuplot_ctrl *  handle, char const *  cmd, ...);

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
void gnuplot_setstyle(gnuplot_ctrl * h, char * plot_style);

/**
  @brief    Change the plotting color of a gnuplot session.
  @param    h Gnuplot session control handle
  @param    plot_color Plotting-color to use (int)
  @return   void
 */
void gnuplot_setcolor(gnuplot_ctrl * h, int plot_style);

/**
  @brief    Sets the x label of a gnuplot session.
  @param    h Gnuplot session control handle.
  @param    label Character string to use for X label.
  @return   void

  Sets the x label for a gnuplot session.
 */
void gnuplot_set_xlabel(gnuplot_ctrl * h, char * label);

/**
  @brief    Sets the y label of a gnuplot session.
  @param    h Gnuplot session control handle.
  @param    label Character string to use for Y label.
  @return   void

  Sets the y label for a gnuplot session.
 */
void gnuplot_set_ylabel(gnuplot_ctrl * h, char * label);

/**
  @brief    Sets the z label of a gnuplot session.
  @param    h Gnuplot session control handle.
  @param    label Character string to use for Z label.
  @return   void

  Sets the z label for a gnuplot session.
 */
void gnuplot_set_zlabel(gnuplot_ctrl * h, char * label);

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
    char            *   title);

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
    int lc);

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
    int lc);

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
    char const      * tmpfname);

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
    char const      * tmpfname);

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
    char const      * tmpfname);

/**
 *  \brief Generates a temporary txt file for gnuplot use.
 **/
char const * gnuplot_tmpfile(gnuplot_ctrl * handle);


/**
 *   \brief Generates a gnuplot command in order to plot the 2D data stored in the
 *          temporary file tmp_filename.
 **/
void gnuplot_plot_atmpfile(gnuplot_ctrl * handle, char const* tmp_filename,
                           char const* title, char const* ls, char const* lt,
                           char const* lw, int lc);

/**
 *   \brief Same as gnuplot_plot_atmpfile for 3D plots.
 **/
void gnuplot_plot_atmpfile_3d(gnuplot_ctrl * handle,
                              char const* tmp_filename,
                              char const* title,
                              char const* ls,
                              char const* lt,
                              char const* lw,
                              int lc);
#endif
