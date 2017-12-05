#ifndef CONFIG_H
#define CONFIG_H

#include <iostream>
#include <iomanip>
#include <string>
#include <cstdio>

/**
 * \file config.h
 * \brief List of some parameters that are used throughout numerical compuration,
 *        such as precision in numerical integrators. Defines a Config Manager that allows
 *        to use common values accross the software.
 * \author BLB
 */


//------------------------------------------------------------------------------------
//   Global constants - such constants should be avoided and, in the long rung, all the
//   parameters declared hereafter should become private variables of the class Config.
//   However, it works well for now, and it certainly lightens the notations, code wise.
//   Indeed, a simple call to
//                      'OFTS_ORDER'
//    would become something like
//                      'Config::configManager().get_OFTS_ORDER()'
//------------------------------------------------------------------------------------
extern int OFTS_ORDER;  // Order of the Taylor series in Fourier-Taylor series
extern int OFS_ORDER;   // Order of the Fourier series in Fourier-Taylor series
extern int MODEL_TYPE;  // Type of model (chosen in the constants beginning by "M_")
extern int REDUCED_NV;  // Number of reduced variables (e.g. 4 for a center manifold, 5 for a center-stable...)


/**
 *  \brief Singleton class that initializes user parameters and allows
 *         for their common use accross the software.
 **/
class Config
{
    private:

        //--------------------------------------------------------------------------------
        // Parameters for ODE structure (see ode.h & cpp)
        //--------------------------------------------------------------------------------
        double PREC_HSTART;  //Initial step in numerical integrator
        double PREC_ABS;     //Absolute precision in numerical integrator
        double PREC_REL;     //Relative precision in numerical integrator
        double PREC_ROOT;    //Precision on root (zero) finding
        double PREC_LIB;     //Precision on the position of the libration point (maybe redundant with _ROOT)
        double DELTA_T;      //Very small delta of time necessary to avoid some errors in numerical procedure at t = 0.0.

        //--------------------------------------------------------------------------------
        // Parameters for differential correction procedures
        //--------------------------------------------------------------------------------
        int DC_ITERMAX;

        //--------------------------------------------------------------------------------
        // Parameters for aesthetics (cout, plotting)
        //--------------------------------------------------------------------------------
        int COUT_SMALL_PREC;    //Small cout precision
        int COUT_MEDIUM_PREC;   //Medium cout precision
        int COUT_LARGE_PREC;    //Large cout precision

    public:
        static Config& configManager()
        {
            static Config instance;
            return instance;
        }

        //--------------------------------------------------------------------------------
        // Getters
        //--------------------------------------------------------------------------------
        double G_PREC_ABS(){return PREC_ABS;}
        double G_PREC_REL(){return PREC_REL;}
        double G_PREC_ROOT(){return PREC_ROOT;}
        double G_PREC_LIB(){return PREC_LIB;}
        double G_PREC_HSTART(){return PREC_HSTART;}
        double G_DELTA_T(){return DELTA_T;}
        int    G_DC_ITERMAX(){return DC_ITERMAX;}

        //--------------------------------------------------------------------------------
        // Setters
        //--------------------------------------------------------------------------------
        /**
         *  \brief Hard precision in numerical integration
         **/
        void C_PREC_HARD();

        /**
         *  \brief Soft precision in numerical integration, to fasten correction loops
         *         It is the default case.
         **/
        void C_PREC_SOFT();

        //--------------------------------------------------------------------------------
        // Aesthetics
        //--------------------------------------------------------------------------------
        /**
         *  \brief Sets a big precision in cout.
         **/
        void coutlp();

        /**
         *  \brief Sets an average precision in cout.
         **/
        void coutmp();

        /**
         *  \brief Sets a small precision in cout.
         **/
        void coutsp();


    private:
        Config();
        Config(Config const&);
        ~Config();
        void operator=(Config const&);
};


/**
 *  \brief Prompt "Press Enter to go on"
 **/
void pressEnter(bool isFlag);

/**
 *  \brief Prompt msg
 **/
void pressEnter(bool isFlag, std::string msg);





#endif // CONFIG_H
