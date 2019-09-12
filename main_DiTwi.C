// $Id$
//==============================================================================
//!
//! \file main_DiTwi.C
//!
//! \date Sep 9 2019
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Main program for digital twin demo.
//!
//==============================================================================

#include "IFEM.h"
#include "NewmarkSIM.h"
#include "NewmarkDriver.h"
#include "SIMLinElModal.h"
#include "SIMargsBase.h"
#include "Utilities.h"
#include "VTF.h"
#include "Profiler.h"
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>


/*!
  \brief Driver for modal analysis of linear dynamic problems.
*/

class ModalDriver : public NewmarkDriver<NewmarkSIM>, public ControlCallback
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  //! \param sim Reference to the spline FE model
  explicit ModalDriver(SIMLinElModal<SIM3D>& sim) :
    NewmarkDriver<NewmarkSIM>(sim), model(sim) {}

  //! \brief Empty destructor.
  virtual ~ModalDriver() {}

  //! \brief Returns the number of solution vectors.
  virtual size_t numSolution() const
  {
    return dynamic_cast<SIMmodal*>(&model)->numExpSolution();
  }

  //! \brief Calculates the current real solution vectors.
  virtual const Vectors& realSolutions()
  {
    return dynamic_cast<SIMmodal*>(&model)->expandSolution(solution,true);
  }

  //! \brief Returns a const reference to the current real solution vector.
  virtual const Vector& realSolution(int i = 0) const
  {
    return dynamic_cast<SIMmodal*>(&model)->expandedSolution(i);
  }

  int solveProblem()
  {
    // Invoke the time-step loop
    int status = 0;
    while (status == 0 && this->advanceStep(params))
    {
      auto backupSol = this->solution;
      while (!model.m_advance) {
        // Solve the dynamic FE problem at this time step
        if (model.m_solve) {
          if (this->solveStep(params,SIM::DYNAMIC,1e-6,6) != SIM::CONVERGED)
          {
            status = 5;
            break;
          }
          // Print solution components at the user-defined points
          this->dumpResults(params.time.t,IFEM::cout,16);
          IFEM::cout << "--- end of iteration ---" << std::endl;
        }

        model.m_solve = false;

        IFEM::pollControllerFifo();

        if (!model.m_advance) {
          this->solution = backupSol;
          params.time.it = 0;
        }
      }
      model.m_advance = false;
    }

    return 0;
  }

  void OnControl(const TiXmlElement* elem) override
  {
    const TiXmlElement* child = elem->FirstChildElement();
    const char* value = nullptr;
    for (; child; child = child->NextSiblingElement())
      if (!strcasecmp(child->Value(),"new_load")) {
        value = utl::getValue(child, "new_load");
        model.setLoad(1000000, value);
        model.m_solve = true;
      } else if (!strcasecmp(child->Value(), "step_ok"))
        model.m_advance = true;
  }

  std::string GetContext() const override { return "ditwi"; }

  //! \brief Dumps solution variables at user-defined points.
  //! \param[in] time Current time
  //! \param[in] os The output stream to write the solution to
  //! \param[in] precision Number of digits after the decimal point
  //! \param[in] formatted If \e false, write all result points on a single line
  void dumpResults(double time, utl::LogStream& os,
                   std::streamsize precision = 3,
                   bool formatted = true) const override
  {
    model.dumpResults(this->realSolution(),time,os,formatted,precision);
  }

  SIMLinElModal<SIM3D>& model;
};


/*!
  \brief Main program for the NURBS-based isogeometric linear elasticity solver.

  The input to the program is specified through the following
  command-line arguments. The arguments may be given in arbitrary order.

  \arg \a input-file : Input file with model definition
  \arg -dense :   Use the dense LAPACK matrix equation solver
  \arg -spr :     Use the SPR direct equation solver
  \arg -superlu : Use the sparse SuperLU equation solver
  \arg -samg :    Use the sparse algebraic multi-grid equation solver
  \arg -petsc :   Use equation solver from PETSc library
  \arg -lag : Use Lagrangian basis functions instead of splines/NURBS
  \arg -spec : Use Spectral basis functions instead of splines/NURBS
  \arg -LR : Use LR-spline basis functions instead of tensorial splines/NURBS
  \arg -nGauss \a n : Number of Gauss points over a knot-span in each direction
  \arg -vtf \a format : VTF-file format (-1=NONE, 0=ASCII, 1=BINARY)
  \arg -nviz \a nviz : Number of visualization points over each knot-span
  \arg -nu \a nu : Number of visualization points per knot-span in u-direction
  \arg -nv \a nv : Number of visualization points per knot-span in v-direction
  \arg -nw \a nw : Number of visualization points per knot-span in w-direction
  \arg -hdf5 : Write primary and projected secondary solution to HDF5 file
  \arg -ztol \a eps : Zero tolerance for printing of solution components
  \arg -ignore \a p1, \a p2, ... : Ignore these patches in the analysis
  \arg -eig \a iop : Eigenproblem solver to use (1...6)
  \arg -nev \a nev : Number of eigenvalues to compute
  \arg -ncv \a ncv : Number of Arnoldi vectors to use in the eigenvalue analysis
  \arg -shift \a shf : Shift value to use in the eigenproblem solver
  \arg -free : Ignore all boundary conditions (use in free vibration analysis)
  \arg -dynamic : Solve the dynamic problem using modal transformation
  \arg -check : Data check only, read model and output to VTF (no solution)
  \arg -checkRHS : Check that the patches are modelled in a right-hand system
  \arg -fixDup : Resolve co-located nodes by merging them into a single node
  \arg -DGL2 : Estimate error using discrete global L2 projection
  \arg -CGL2 : Estimate error using continuous global L2 projection
  \arg -SCR : Estimate error using Superconvergent recovery at Greville points
  \arg -VDSA: Estimate error using Variational Diminishing Spline Approximations
  \arg -LSQ : Estimate error using through Least Square projections
  \arg -QUASI : Estimate error using Quasi-interpolation projections
*/


int main (int argc, char** argv)
{
  Profiler prof(argv[0]);
  utl::profiler->start("Initialization");

  std::vector<std::string> topSets;
  std::vector<int> ignoredPatches;
  int  i;
  char* infile = nullptr;
  Elasticity::wantPrincipalStress = true;
  SIMargsBase args("elasticity");

  IFEM::Init(argc,argv,"Linear Elasticity solver");

  for (i = 1; i < argc; i++)
    if (argv[i] == infile || args.parseArg(argv[i]))
      ; // ignore the input file on the second pass
    else if (SIMoptions::ignoreOldOptions(argc,argv,i))
      ; // ignore the obsolete option
    else if (!infile)
    {
      infile = argv[i];
      if (strcasestr(infile,".xinp"))
      {
        if (!args.readXML(infile,false))
          return 1;
        i = 0; // start over and let command-line options override input file
      }
    }
    else
      std::cerr <<"  ** Unknown option ignored: "<< argv[i] << std::endl;

  if (!infile)
  {
    std::cout <<"usage: "<< argv[0] <<" <inputfile>\n";
    return 0;
  }

  IFEM::cout <<"\nInput file: "<< infile;
  IFEM::getOptions().print(IFEM::cout);
  utl::profiler->stop("Initialization");
  utl::profiler->start("Model input");

  // Create the simulation model
  std::vector<Mode> modes;

  SIMLinElModal<SIM3D> model(modes, false);
  ModalDriver solver(model);

  // Read in model definitions
  if (!solver.read(infile))
    return 1;

  model.opt.eig = 4;

  model.opt.print(IFEM::cout,true) << std::endl;

  utl::profiler->stop("Model input");

  if (!model.preprocess())
    return 7;

  model.setMode(SIM::VIBRATION);
  model.setQuadratureRule(model.opt.nGauss[0],true,true);
  model.initSystem(model.opt.solver,2,0);
  if (!model.assembleSystem())
    return 8;

  // Solve the generalized eigenvalue problem
  if (!model.systemModes(modes))
    return 9;

  /*for (const Mode& mode : modes) {
    size_t patch = 1;
    for (const ASMbase* pch : model.getFEModel()) {
      Vector locVec;
      model.extractPatchSolution(mode.eigVec, locVec, pch, 3);
      std::stringstream str;
      str << "mode_" << mode.eigNo <<"_patch" << patch << ".asc";
      std::ofstream of(str.str());
      for (size_t i = 0; i < locVec.size() / 3; ++i)
        of << locVec[i*3] << " " << locVec[i*3+1] << " " << locVec[i*3+2] << std::endl;
      ++patch;
    }
  }
*/

  // Initialize the modal time-domain simulation
  solver.initPrm();
  solver.initSolution(modes.size(),3);
  solver.printProblem();

  IFEM::registerCallback(solver);

  // Initialize the linear equation system.
  // Actually, we don't need any system matrices here
  // since we are only integrating the external load vector in space.
  if (!model.initSystem(LinAlg::DENSE,0,1))
    return 3;

  Elasticity::wantStrain = true;
  return solver.solveProblem();
}
