//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file reconnection.cpp
//  \brief Problem generator for magnetic reconnection problem.  Works in Cartesian,
//         coordinate.
//
// REFERENCE:

// C++ headers
#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../utils/utils.hpp" // ran2()

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

#undef RELATIVISTIC_DYNAMICS

// initial pressure
static Real pres_init(const Real bx, const Real by, const Real bz, const int pres_balance);

// functions to compute vector potential to initialize the solution
static Real A1(const Real x, const Real y, const Real z);
static Real A2(const Real x, const Real y, const Real z);
static Real A3(const Real x, const Real y, const Real z);
static Real A3_pert(const Real x, const Real y, const Real z);

/* // Vorticity */
/* Real ***wx, ***wy, ***wz; */

/* // Apply a density floor - useful for large |z| regions */
/* static Real D_FLOOR = 1.e-3; */

// Parameters in initial configurations
static int forcefree, num_cs, pres_balance, uniform_rho;
static Real cs_width, beta0;
static Real b0 = 1.0, Bguide; // Normalized magnetic field
static int pert_B, pert_V, random_vpert; // Perturbation
static Real phi_pert, vin_pert;
static Real xmin, xmax, lx; // Grid dimensions
static Real ymin, ymax, ly;
static Real xleft, xright, xcenter;
static Real gamma_adi, gamma_adi_red, gm1, iso_cs;

#ifdef RELATIVISTIC_DYNAMICS
static Real rho, pgas;
#endif

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief magnetic reconnection problem generator
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  forcefree = 1;    // In default, forcefree current sheet is used
  num_cs = 2;       // In default, two current sheets are used
  pres_balance = 1; // In default, the initial pressure is balanced
  uniform_rho = 0;  // In default, the initial density is not uniform
  pert_B = 1;       // In default, magnetic field is perturbed
  pert_V = 0;       // In default, velocity field is not perturbed
  random_vpert = 1; // In default, random velocity perturbation is used

  xmin  = pin->GetReal("mesh", "x1min");
  xmax  = pin->GetReal("mesh", "x1max");
  ymin  = pin->GetReal("mesh", "x2min");
  ymax  = pin->GetReal("mesh", "x2max");
  beta0 = pin->GetReal("problem", "beta0");
  vin_pert = pin->GetReal("problem", "vin_pert");
  random_vpert = pin->GetInteger("problem", "random_vpert");
  forcefree = pin->GetInteger("problem", "forcefree");
  cs_width = pin->GetReal("problem", "cs_width");
  Bguide   = pin->GetReal("problem", "Bguide");
  num_cs   = pin->GetInteger("problem", "num_cs");
  phi_pert = pin->GetReal("problem", "phi_pert");
  pres_balance = pin->GetInteger("problem", "pres_balance");
  uniform_rho = pin->GetInteger("problem", "uniform_rho");
  pert_B = pin->GetInteger("problem", "pert_B");
  pert_V = pin->GetInteger("problem", "pert_V");

  // initialize global variables
  if (NON_BAROTROPIC_EOS) {
    gamma_adi = peos->GetGamma();
    gamma_adi_red = gamma_adi / (gamma_adi - 1.0);
    gm1 = (gamma_adi - 1.0);
  } else {
    iso_cs = pin->GetReal("hydro", "iso_sound_speed");
  }

#ifdef RELATIVISTIC_DYNAMICS
  b0   = pin->GetReal("problem", "b0");
  rho  = pin->GetReal("problem", "rho");
  pgas = pin->GetReal("problem", "pgas");
  Real sigma = b0 * b0 / (rho + gamma_adi_red * pgas);
  if (Globals::my_rank == 0) {
    std::cout << "Magnetization (sigma) is: " << sigma << std::endl;
  }
#endif

  lx  = xmax - xmin;
  ly  = ymax - ymin;

  AthenaArray<Real> a1, a2, a3;
  AthenaArray<Real> b1i_pert, b2i_pert, b1c_pert, b2c_pert;
  int nx1 = (ie-is)+1 + 2*(NGHOST);
  int nx2 = (je-js)+1 + 2*(NGHOST);
  int nx3 = (ke-ks)+1 + 2*(NGHOST);
  a1.NewAthenaArray(nx3,nx2,nx1);
  a2.NewAthenaArray(nx3,nx2,nx1);
  a3.NewAthenaArray(nx3,nx2,nx1);

  // Positions of the current sheets
  xleft   = xmin + 0.25 * lx;
  xright  = xmin + 0.75 * lx;
  xcenter = (xmax + xmin) * 0.5;

  int level=loc.level;
  // Initialize components of the vector potential
  if (block_size.nx3 > 1) {
    for (int k=ks; k<=ke+1; k++) {
      for (int j=js; j<=je+1; j++) {
        for (int i=is; i<=ie+1; i++) {
          if ((pbval->nblevel[1][0][1]>level && j==js)
           || (pbval->nblevel[1][2][1]>level && j==je+1)
           || (pbval->nblevel[0][1][1]>level && k==ks)
           || (pbval->nblevel[2][1][1]>level && k==ke+1)
           || (pbval->nblevel[0][0][1]>level && j==js   && k==ks)
           || (pbval->nblevel[0][2][1]>level && j==je+1 && k==ks)
           || (pbval->nblevel[2][0][1]>level && j==js   && k==ke+1)
           || (pbval->nblevel[2][2][1]>level && j==je+1 && k==ke+1)) {
            Real x1l = pcoord->x1f(i)+0.25*pcoord->dx1f(i);
            Real x1r = pcoord->x1f(i)+0.75*pcoord->dx1f(i);
            a1(k,j,i) = 0.5*(A1(x1l, pcoord->x2f(j), pcoord->x3f(k)) +
                             A1(x1r, pcoord->x2f(j), pcoord->x3f(k)));
          } else {
            a1(k,j,i) = A1(pcoord->x1v(i), pcoord->x2f(j), pcoord->x3f(k));
          }

          if ((pbval->nblevel[1][1][0]>level && i==is)
           || (pbval->nblevel[1][1][2]>level && i==ie+1)
           || (pbval->nblevel[0][1][1]>level && k==ks)
           || (pbval->nblevel[2][1][1]>level && k==ke+1)
           || (pbval->nblevel[0][1][0]>level && i==is   && k==ks)
           || (pbval->nblevel[0][1][2]>level && i==ie+1 && k==ks)
           || (pbval->nblevel[2][1][0]>level && i==is   && k==ke+1)
           || (pbval->nblevel[2][1][2]>level && i==ie+1 && k==ke+1)) {
            Real x2l = pcoord->x2f(j)+0.25*pcoord->dx2f(j);
            Real x2r = pcoord->x2f(j)+0.75*pcoord->dx2f(j);
            a2(k,j,i) = 0.5*(A2(pcoord->x1f(i), x2l, pcoord->x3f(k)) +
                             A2(pcoord->x1f(i), x2r, pcoord->x3f(k)));
          } else {
            a2(k,j,i) = A2(pcoord->x1f(i), pcoord->x2v(j), pcoord->x3f(k));
          }

          if ((pbval->nblevel[1][1][0]>level && i==is)
           || (pbval->nblevel[1][1][2]>level && i==ie+1)
           || (pbval->nblevel[1][0][1]>level && j==js)
           || (pbval->nblevel[1][2][1]>level && j==je+1)
           || (pbval->nblevel[1][0][0]>level && i==is   && j==js)
           || (pbval->nblevel[1][0][2]>level && i==ie+1 && j==js)
           || (pbval->nblevel[1][2][0]>level && i==is   && j==je+1)
           || (pbval->nblevel[1][2][2]>level && i==ie+1 && j==je+1)) {
            Real x3l = pcoord->x3f(k)+0.25*pcoord->dx3f(k);
            Real x3r = pcoord->x3f(k)+0.75*pcoord->dx3f(k);
            a3(k,j,i) = 0.5*(A3(pcoord->x1f(i), pcoord->x2f(j), x3l) +
                             A3(pcoord->x1f(i), pcoord->x2f(j), x3r));
            if (pert_B == 1)
              a3(k,j,i) += 0.5*(A3_pert(pcoord->x1f(i), pcoord->x2f(j), x3l) +
                                A3_pert(pcoord->x1f(i), pcoord->x2f(j), x3r));
          } else {
            a3(k,j,i) = A3(pcoord->x1f(i), pcoord->x2f(j), pcoord->x3v(k));
            a3(k,j,i) += A3_pert(pcoord->x1f(i), pcoord->x2f(j), pcoord->x3v(k));
          }
        }
      }
    }
  } else {
    for (int k=ks; k<=ke+1; k++) {
      for (int j=js; j<=je+1; j++) {
        for (int i=is; i<=ie+1; i++) {
          if (i != ie+1)
            a1(k,j,i) = A1(pcoord->x1v(i), pcoord->x2f(j), pcoord->x3f(k));
          if (j != je+1)
            a2(k,j,i) = A2(pcoord->x1f(i), pcoord->x2v(j), pcoord->x3f(k));
          if (k != ke+1) {
            a3(k,j,i) = A3(pcoord->x1f(i), pcoord->x2f(j), pcoord->x3v(k));
            a3(k,j,i) += A3_pert(pcoord->x1f(i), pcoord->x2f(j), pcoord->x3v(k));
          }
        }
      }
    }
  }

  // Initialize interface fields
  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie+1; i++) {
      pfield->b.x1f(k,j,i) = (a3(k  ,j+1,i) - a3(k,j,i))/pcoord->dx2f(j) -
                             (a2(k+1,j  ,i) - a2(k,j,i))/pcoord->dx3f(k);
    }
  }}

  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je+1; j++) {
    for (int i=is; i<=ie; i++) {
      pfield->b.x2f(k,j,i) = (a1(k+1,j,i  ) - a1(k,j,i))/pcoord->dx3f(k) -
                             (a3(k  ,j,i+1) - a3(k,j,i))/pcoord->dx1f(i);
    }
  }}

  for (int k=ks; k<=ke+1; k++) {
  for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie; i++) {
      pfield->b.x3f(k,j,i) = (a2(k,j  ,i+1) - a2(k,j,i))/pcoord->dx1f(i) -
                             (a1(k,j+1,i  ) - a1(k,j,i))/pcoord->dx2f(j);
    }
  }}

  // Compute cell-centered fields
  pfield->CalculateCellCenteredField(pfield->b, pfield->bcc, pcoord, is, ie, js, je, ks, ke);

  // Initialize hydro variables
#ifdef RELATIVISTIC_DYNAMICS
  AthenaArray<Real> bb;
  bb.NewAthenaArray(3, ke+1, je+1, ie+1);
  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie; i++) {
      // Set primitives
      phydro->w(IDN,k,j,i) = phydro->w1(IDN,k,j,i) = rho;
      phydro->w(IPR,k,j,i) = phydro->w1(IPR,k,j,i) = pgas;
      phydro->w(IVX,k,j,i) = phydro->w1(IM1,k,j,i) = 0.0;
      phydro->w(IVY,k,j,i) = phydro->w1(IM2,k,j,i) = 0.0;
      phydro->w(IVZ,k,j,i) = phydro->w1(IM3,k,j,i) = 0.0;
      // Set magnetic fields
      bb(IB1,k,j,i) = pfield->b.x1f(k,j,i);
      bb(IB2,k,j,i) = pfield->b.x2f(k,j,i);
      bb(IB3,k,j,i) = pfield->b.x3f(k,j,i);
    }
  }}
  peos->PrimitiveToConserved(phydro->w, bb, phydro->u, pcoord, is, ie, js, je, ks, ke);
  bb.DeleteAthenaArray();
#else
  Real pgas_nr;
  int64_t iseed = -1 - gid; // Ensure a different initial random seed for each meshblock.
  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie; i++) {
      pgas_nr = pres_init(pfield->bcc(IB1,k,j,i), pfield->bcc(IB2,k,j,i),
          pfield->bcc(IB3,k,j,i), pres_balance);
      // density
      if (uniform_rho)
        phydro->u(IDN,k,j,i) = b0 * b0;
      else
        phydro->u(IDN,k,j,i) = pgas_nr / (0.5 * beta0);
      // momentum
      if (pert_V) {
        if (random_vpert) {
          phydro->u(IM1,k,j,i) = vin_pert * (2.0*ran2(&iseed) - 1.0);
          phydro->u(IM2,k,j,i) = vin_pert * (2.0*ran2(&iseed) - 1.0);
        } else {
          Real x1 = pcoord->x1v(i);
          Real x2 = pcoord->x2v(j);
          if (num_cs == 1) {
            phydro->u(IM1,k,j,i) = -vin_pert * 2 * PI * cos(PI*(x1-xcenter)/lx) * sin(2*PI*x2/ly) / ly;
            phydro->u(IM2,k,j,i) = vin_pert * PI * sin(PI*(x1-xcenter)/lx) * cos(2*PI*x2/ly) / lx;
          } else if (num_cs == 2) {
            phydro->u(IM1,k,j,i) = -vin_pert * 2 * PI * (cos(2.0*PI*(x1-xleft)/lx) -
                cos(2.0*PI*(x1-xright)/lx)) * sin(2*PI*x2/ly) / ly;
            phydro->u(IM2,k,j,i) = vin_pert * 2 * PI * (sin(2.0*PI*(x1-xleft)/lx) -
                sin(2.0*PI*(x1-xright)/lx)) * cos(2*PI*x2/ly) / lx;
          } else {
            std::stringstream msg;
            msg << "### FATAL ERROR in reconnection.cpp ProblemGenerator" << std::endl
                << "More than 2 current sheets is not supported now" << std::endl;
            throw std::runtime_error(msg.str().c_str());
          }
        }
      } else {
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
      }
      phydro->u(IM1,k,j,i) *= phydro->u(IDN,k,j,i);
      phydro->u(IM2,k,j,i) *= phydro->u(IDN,k,j,i);
      phydro->u(IM3,k,j,i) = 0.0;

      if (NON_BAROTROPIC_EOS) {
        phydro->u(IEN,k,j,i) = pgas_nr/gm1;
        phydro->u(IEN,k,j,i) += 0.5*(SQR(phydro->u(IM1,k,j,i)) +
                                     SQR(phydro->u(IM2,k,j,i)) +
                                     SQR(phydro->u(IM3,k,j,i))) / phydro->u(IDN,k,j,i);
        phydro->u(IEN,k,j,i) += 0.5*(SQR(pfield->bcc(IB1,k,j,i)) +
                                     SQR(pfield->bcc(IB2,k,j,i)) +
                                     SQR(pfield->bcc(IB3,k,j,i)));
      }
    }
  }}
#endif

  a1.DeleteAthenaArray();
  a2.DeleteAthenaArray();
  a3.DeleteAthenaArray();

  return;
}

//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//  \brief
//========================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
  return;
}

int RefinementCondition(MeshBlock *pmb);

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  if(adaptive==true)
    EnrollUserRefinementCondition(RefinementCondition);

  return;
}

// refinement condition: density and pressure curvature
int RefinementCondition(MeshBlock *pmb)
{
  AthenaArray<Real> &w = pmb->phydro->w;
  Real maxeps=0.0;
  for(int k=pmb->ks; k<=pmb->ke; k++) {
    for(int j=pmb->js; j<=pmb->je; j++) {
      for(int i=pmb->is; i<=pmb->ie; i++) {
        Real epsr= (std::abs(w(IDN,k,j,i+1)-2.0*w(IDN,k,j,i)+w(IDN,k,j,i-1))
                   +std::abs(w(IDN,k,j+1,i)-2.0*w(IDN,k,j,i)+w(IDN,k,j-1,i))
                   +std::abs(w(IDN,k+1,j,i)-2.0*w(IDN,k,j,i)+w(IDN,k-1,j,i)))/w(IDN,k,j,i);
        Real epsp= (std::abs(w(IPR,k,j,i+1)-2.0*w(IPR,k,j,i)+w(IPR,k,j,i-1))
                   +std::abs(w(IPR,k,j+1,i)-2.0*w(IPR,k,j,i)+w(IPR,k,j-1,i))
                   +std::abs(w(IPR,k+1,j,i)-2.0*w(IPR,k,j,i)+w(IPR,k-1,j,i)))/w(IPR,k,j,i);
        Real eps = std::max(epsr, epsp);
        maxeps = std::max(maxeps, eps);
      }
    }
  }
  /* // refine : curvature > 0.01 */
  /* if(maxeps > 0.01) return 1; */
  /* // derefinement: curvature < 0.005 */
  /* if(maxeps < 0.005) return -1; */
  /* // otherwise, stay */
  // refine : curvature > 10.0
  if(maxeps > 100.0) return 1;
  // derefinement: curvature < 5.0
  if(maxeps < 10.0) return -1;
  // otherwise, stay
  return 0;
}

/* // refinement condition: check the pressure gradient */
/* int RefinementCondition(MeshBlock *pmb) */
/* { */
/*   AthenaArray<Real> &w = pmb->phydro->w; */
/*   Real maxeps=0.0; */
/*   Real dx = pmb->pcoord->dx1f(0); */
/*   for(int k=pmb->ks; k<=pmb->ke; k++) { */
/*     for(int j=pmb->js; j<=pmb->je; j++) { */
/*       for(int i=pmb->is; i<=pmb->ie; i++) { */
/*         Real eps= sqrt(SQR(0.5*(w(IPR,k,j,i+1)-w(IPR,k,j,i-1))) */
/*                       +SQR(0.5*(w(IPR,k,j+1,i)-w(IPR,k,j-1,i))) */
/*                       +SQR(0.5*(w(IPR,k+1,j,i)-w(IPR,k-1,j,i))))/w(IPR,k,j,i); */
/*         maxeps = std::max(maxeps, eps); */
/*       } */
/*     } */
/*   } */
/*   if(maxeps > 10.0) return 1; */
/*   if(maxeps < 2.0) return -1; */
/*   return 0; */
/* } */

//========================================================================================
//! \fn static Real pres_init(const Real bx, const Real by, const Real bz, const int pres_balance)
//  \brief initial gas pressure
//========================================================================================
static Real pres_init(const Real bx, const Real by, const Real bz, const int pres_balance)
{
  Real p0, pB, pmag_max;
  if (pres_balance && pert_B) {
    if (num_cs == 1) {
      pmag_max = 0.5 * b0 * b0 * (1.0+phi_pert*PI/lx)*(1.0+phi_pert*PI/lx) + 0.5 * Bguide * Bguide;
    } else if (num_cs == 2) {
      pmag_max = 0.5 * b0 * b0 * (1.0+4.0*phi_pert*PI/lx)*(1.0+4.0*phi_pert*PI/lx) + 0.5 * Bguide * Bguide;
    } else {
      std::stringstream msg;
      msg << "### FATAL ERROR in reconnection.cpp pres_init" << std::endl
          << "More than 2 current sheets is not supported now" << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
    pB = 0.5 * (bx*bx + by*by + bz*bz);
    p0 = (beta0 * pmag_max + pmag_max) - pB;
  } else {
    p0 = beta0 * (0.5 * b0 * b0);
  }
  return p0;
}

//========================================================================================
//! \fn static Real A1(const Real x, const Real y, const Real z)
//! \fn static Real A2(const Real x, const Real y, const Real z)
//! \fn static Real A3(const Real x, const Real y, const Real z)
//  \brief three components of the vector potential
//========================================================================================

static Real A1(const Real x, const Real y, const Real z)
{
  return 0.0;
}

static Real A2(const Real x, const Real y, const Real z)
{
  Real pix, Ay, a, a2, s2, b2;

  a = Bguide / b0;
  a2 = a * a;
  s2 = sqrt(2);

  if (num_cs == 1) {
    if (forcefree == 0) {
      Ay = Bguide * x;
    } else {
      pix = (x - xcenter) / cs_width;
      b2 = a2*cosh(2*pix) + a2 + 2;
      Ay = b0 * cs_width * (1.0/b2) * s2 * cosh(pix) *
          sqrt(a2 + 1.0/(cosh(pix)*cosh(pix))) *
          (sqrt(a2*(a2 + 1)) * asinh(sqrt(a2/(a2+1)) * sinh(pix)) *
           sqrt(b2 / (a2 + 1)) + sqrt(b2) * atan(s2*sinh(pix) / sqrt(b2)));
    }
  } else if (num_cs == 2) {
    if (forcefree == 0) {
      Ay = Bguide * x;
    } else {
      pix = (x - xleft) / cs_width;
      b2 = a2*cosh(2*pix) + a2 + 2;
      Ay = b0 * cs_width * (1.0/b2) * s2 * cosh(pix) *
          sqrt(a2 + 1.0/(cosh(pix)*cosh(pix))) *
          (sqrt(a2*(a2 + 1)) * asinh(sqrt(a2/(a2+1)) * sinh(pix)) *
           sqrt(b2 / (a2 + 1)) + sqrt(b2) * atan(s2*sinh(pix) / sqrt(b2)));
      pix = (x - xright) / cs_width;
      b2 = a2*cosh(2*pix) + a2 + 2;
      Ay += b0 * cs_width * (1.0/b2) * s2 * cosh(pix) *
          sqrt(a2 + 1.0/(cosh(pix)*cosh(pix))) *
          (sqrt(a2*(a2 + 1)) * asinh(sqrt(a2/(a2+1)) * sinh(pix)) *
           sqrt(b2 / (a2 + 1)) + sqrt(b2) * atan(s2*sinh(pix) / sqrt(b2)));
      Ay -= Bguide * x;
    }
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in reconnection.cpp A2" << std::endl
        << "More than 2 current sheets is not supported now" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  return Ay;
}

static Real A3(const Real x, const Real y, const Real z)
{
  Real Az;

  if (num_cs == 1) {
    Az = -b0 * cs_width * log(cosh((x-xcenter)/cs_width));
  } else if (num_cs == 2) {
    Az = -b0 * cs_width * (log(cosh((x-xleft)/cs_width)) -
                           log(cosh((x-xright)/cs_width))) + b0 * x;
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in reconnection.cpp A3" << std::endl
        << "More than 2 current sheets is not supported now" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  return Az;
}

//========================================================================================
//! \fn static Real A3_pert(const Real x, const Real y, const Real z)
//  \brief A3_pert: 3-component of perturbed vector potential
//========================================================================================
static Real A3_pert(const Real x, const Real y, const Real z)
{
  Real dAz;

  if (num_cs == 1) {
    dAz = phi_pert * cos(PI*(x-xcenter)/lx) * cos(2*PI*y/ly);
  } else if (num_cs == 2) {
    dAz = phi_pert * (cos(2.0*PI*(x-xleft)/lx) - cos(2.0*PI*(x-xright)/lx)) * cos(2*PI*y/ly);
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in reconnection.cpp A3_pert" << std::endl
        << "More than 2 current sheets is not supported now" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  return dAz;
}
