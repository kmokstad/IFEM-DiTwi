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
#include "ModalDriver.h"
#include "SIMLinElModal.h"
#include "SIMargsBase.h"
#include "Utilities.h"
#include "Profiler.h"
#include <fstream>
#include <stdlib.h>
#include <string.h>


/*!
  \brief Simulator driver for modal linear elastic digital twin models.
*/

class SIMDigitalTwin : public SIMLinElModal<SIM3D>
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  SIMDigitalTwin(std::vector<Mode>& modes) : SIMLinElModal<SIM3D>(modes) {}
  //! \brief Empty destructor.
  ~SIMDigitalTwin() {}

  //! \brief Redefines the (constant) traction load.
  void setLoad(size_t propInd, const char* value)
  {
    delete SIM3D::myTracs[propInd];
    SIM3D::myTracs[propInd] = utl::parseTracFunc(value,"constant",1);
  }

  bool m_solve   = false;
  bool m_advance = false;
};


/*!
  \brief Driver for digital twin simulators.
*/

class DigTwinDriver : public ModalDriver, public ControlCallback
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  explicit DigTwinDriver(SIMDigitalTwin& sim) : ModalDriver(sim), myModel(sim) {}
  //! \brief Empty destructor.
  virtual ~DigTwinDriver() {}

  //! \brief Invokes the main time stepping simulation loop.
  int solveProblem()
  {
    // Invoke the time-step loop
    int status = 0;
    while (status == 0 && this->advanceStep(params))
    {
      auto backupSol = this->solution;
      while (!myModel.m_advance) {
        // Solve the dynamic FE problem at this time step
        if (myModel.m_solve) {
          if (this->solveStep(params,SIM::DYNAMIC,1e-6,6) != SIM::CONVERGED)
          {
            status = 5;
            break;
          }
          // Print solution components at the user-defined points
          this->dumpResults(params.time.t,IFEM::cout,16,true);
          IFEM::cout << "--- end of iteration ---" << std::endl;
        }

        myModel.m_solve = false;

        IFEM::pollControllerFifo();

        if (!myModel.m_advance) {
          this->solution = backupSol;
          params.time.it = 0;
        }
      }
      myModel.m_advance = false;
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
        myModel.setLoad(1000000, value);
        myModel.m_solve = true;
      } else if (!strcasecmp(child->Value(), "step_ok"))
        myModel.m_advance = true;
  }

  std::string GetContext() const override { return "ditwi"; }

private:
  SIMDigitalTwin& myModel;
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
  \arg -eig \a iop : Eigenproblem solver to use (1...6)
  \arg -nev \a nev : Number of eigenvalues to compute
  \arg -ncv \a ncv : Number of Arnoldi vectors to use in the eigenvalue analysis
  \arg -shift \a shf : Shift value to use in the eigenproblem solver
*/


int main (int argc, char** argv)
{
  Profiler prof(argv[0]);
  utl::profiler->start("Initialization");

  bool dumpModes = false;
  char* infile = nullptr;
  SIMargsBase args("elasticity");

  IFEM::Init(argc,argv,"Linear Elasticity solver");

  for (int i = 1; i < argc; i++)
    if (argv[i] == infile || args.parseArg(argv[i]))
      ; // ignore the input file on the second pass
    else if (SIMoptions::ignoreOldOptions(argc,argv,i))
      ; // ignore the obsolete option
    else if (!strcmp(argv[i],"-principal"))
      Elasticity::wantPrincipalStress = true;
    else if (!strcmp(argv[i],"-dumpModes"))
      dumpModes = true;
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
  SIMDigitalTwin    model(modes);
  DigTwinDriver     solver(model);

  // Read in model definitions
  if (!solver.read(infile))
    return 1;

  model.opt.eig = 4;

  model.opt.print(IFEM::cout,true) << std::endl;

  utl::profiler->stop("Model input");

  if (!model.preprocess())
    return 7;

  // Free vibration: Assemble [Km] and [M]
  model.setMode(SIM::VIBRATION);
  model.setQuadratureRule(model.opt.nGauss[0],true,true);
  model.initSystem(model.opt.solver,2,0);
  if (!model.assembleSystem())
    return 8;

  // Solve the generalized eigenvalue problem
  if (!model.systemModes(modes))
    return 9;

  // Print out control point stresses for the eigenmodes
  if (dumpModes)
    solver.dumpModes(IFEM::cout,10);

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


// Dummy implementations (no analytical solutions here),
template<> bool SIMLinEl3D::parseDimSpecific (char*, std::istream&)
{
  return false;
}

template<> bool SIMLinEl3D::parseDimSpecific (const TiXmlElement*)
{
  return false;
}
