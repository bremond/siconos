/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

/*!\file BallNewtonEulerOnMovingPlane.cpp
  \brief
  A Ball bouncing on the moving ground.
  Direct description of the model.
  Simulation with a Time-Stepping scheme.
*/


#include "SiconosKernel.hpp"

//#define WITH_PROJ
#define WITH_FC3D
using namespace std;
#ifdef WITH_FC3D
#define R_CLASS NewtonEulerFrom3DLocalFrameR
#else
#define R_CLASS NewtonEulerFrom1DLocalFrameR
#endif

class my_NewtonEulerR : public R_CLASS
{

  double _sBallRadius ;

public:

  my_NewtonEulerR(double radius): R_CLASS(), _sBallRadius(radius) { };

  virtual void computeOutput(double t, Interaction& inter, unsigned int derivativeNumber)
  {
    VectorOfBlockVectors& DSlink = inter.linkToDSVariables();
    if (derivativeNumber == 0)
    {
      computeh(t, *DSlink[NewtonEulerR::q0], *inter.y(0));
    }
    else
    {
      R_CLASS::computeOutput(t, inter, derivativeNumber);
    }

  }
  void computeh(double time, BlockVector& q0, SiconosVector& y)
  {

    //std::cout <<"my_NewtonEulerR:: computeh" << std:: endl;
    //std::cout <<"q0.size() = " << q0.size() << std:: endl;
    double height = q0.getValue(0) - _sBallRadius - q0.getValue(7);

    y.setValue(0, height);
    _Nc->setValue(0, 1);
    _Nc->setValue(1, 0);
    _Nc->setValue(2, 0);
    _Pc1->setValue(0, q0.getValue(0) - _sBallRadius);
    _Pc1->setValue(1, q0.getValue(1));
    _Pc1->setValue(2, q0.getValue(2));

    _Pc2->setValue(0,q0.getValue(7));
    _Pc2->setValue(1,q0.getValue(8));
    _Pc2->setValue(2,q0.getValue(9));
    //printf("my_NewtonEulerR N, Pc\n");
    //_Nc->display();
    //_Pc1->display();
    //_Pc2->display();
    //std::cout <<"my_NewtonEulerR:: computeh ends" << std:: endl;
  }
  //ACCEPT_VISITORS();
};
TYPEDEF_SPTR(my_NewtonEulerR);





int main(int argc, char* argv[])
{
  try
  {


    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int nDof = 3;           // degrees of freedom for the ball
    unsigned int qDim = 7;           // degrees of freedom for the ball
    unsigned int nDim = 6;           // degrees of freedom for the ball
    double t0 = 0;                   // initial computation time
    double T = 10.0;                  // final computation time
    double h = 0.005;                // time step
    double position_init = 1.0;      // initial position for lowest bead.
    double velocity_init = 2.0;      // initial velocity for lowest bead.
    double omega_initx = 0.0;
    double omega_initz = 0.0;// initial velocity for lowest bead.
    double theta = 0.5;              // theta for MoreauJeanOSI integrator
    double m = 1; // Ball mass
    double g = 9.81; // Gravity
    double radius = 0.0;
    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    cout << "====> Model loading ..." << endl << endl;

    // -- Initial positions and velocities --
    SP::SiconosVector q0(new SiconosVector(qDim));
    SP::SiconosVector v0(new SiconosVector(nDim));
    SP::SimpleMatrix I(new SimpleMatrix(3, 3));
    v0->zero();
    q0->zero();

    I->eye();
    (*q0)(0) = position_init;
    /*initial quaternion equal to (1,0,0,0)*/
    (*q0)(3) = 1.0;

    (*v0)(0) = velocity_init;
    (*v0)(3) = omega_initx;
    (*v0)(5) = omega_initz;
    // -- The dynamical system --
    SP::NewtonEulerDS ball(new NewtonEulerDS(q0, v0, m, I));

    // -- Set external forces (weight) --
    SP::SiconosVector weight(new SiconosVector(nDof));
    (*weight)(0) = -m * g;
    ball->setFExtPtr(weight);

    // -- Moving Plane --


    // -- Initial positions and velocities --
    SP::SiconosVector q02(new SiconosVector(qDim));
    SP::SiconosVector v02(new SiconosVector(nDim));
    v02->zero();
    q02->zero();
    (*q02)(3) = 1.0;
    // -- The dynamical system --
    SP::NewtonEulerDS movingplane(new NewtonEulerDS(q02, v02, m, I ));

    // // -- Set external forces (weight) --
    movingplane->setFExtPtr(weight);

    SP::IndexInt bdindex(new IndexInt(1));
    (*bdindex)[0] = 0;

    // SP::SiconosVector bdPrescribedVelocity(new SiconosVector(1));
    // bdPrescribedVelocity->setValue(0,0.5);
    // SP::BoundaryCondition bd (new BoundaryCondition(bdindex,bdPrescribedVelocity));


    SP::BoundaryCondition bd(new BoundaryCondition(bdindex));
    bd->setComputePrescribedVelocityFunction("BallOnMovingPlanePlugin", "prescribedvelocity");

    movingplane->setBoundaryConditions(bd);




    // --------------------
    // --- Interactions ---
    // --------------------

    // -- nslaw --
    double e = 0.9;

    // Interaction ball-floor
    //

    //     vector<SP::SiconosMatrix> vecMatrix1;
    //     vecMatrix1.push_back(H);
    //     SP::BlockMatrix H_block(new BlockMatrix(vecMatrix1,1,1));

    //     SP::SiconosMatrix HT(new SimpleMatrix(1,nDim));
    //     vector<SP::SiconosMatrix> vecMatrix2;
    //     vecMatrix2.push_back(HT);
    //     SP::BlockMatrix HT_block(new BlockMatrix(vecMatrix2,1,1));

#ifdef WITH_FC3D
    int nslawsize = 3;
    SP::NonSmoothLaw nslaw0(new NewtonImpactFrictionNSL(e, e, 0.6, 3));
#else
    int nslawsize = 1;
    SP::NonSmoothLaw nslaw0(new NewtonImpactNSL(e));
#endif

    // Version with my_NewtonEulerR()
    SP::NewtonEulerR relation0(new my_NewtonEulerR(radius));
    SP::Interaction inter(new Interaction(nslaw0, relation0));

    // -------------
    // --- Model ---
    // -------------
    SP::NonSmoothDynamicalSystem bouncingBall(new NonSmoothDynamicalSystem(t0, T));
    // add the dynamical system in the non smooth dynamical system
    bouncingBall->insertDynamicalSystem(ball);
    bouncingBall->insertDynamicalSystem(movingplane);

    // link the interaction and the dynamical system
    //r
    bouncingBall->link(inter, ball, movingplane);

    // ------------------
    // --- Simulation ---
    // ------------------

    // -- (1) OneStepIntegrators --
#ifdef WITH_PROJ
    SP::MoreauJeanDirectProjectionOSI OSI(new MoreauJeanDirectProjectionOSI(theta));
#else
    SP::MoreauJeanOSI OSI(new MoreauJeanOSI(theta));
#endif
    // -- (2) Time discretisation --
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));

    // -- (3) one step non smooth problem
    SP::OneStepNSProblem osnspb(new GenericMechanical());
#ifdef WITH_PROJ
    SP::OneStepNSProblem osnspb_pos(new MLCPProjectOnConstraints(SICONOS_MLCP_ENUM, 1.0));
#endif
    // -- (4) Simulation setup with (1) (2) (3)
#ifdef WITH_PROJ
    SP::TimeSteppingDirectProjection s(new TimeSteppingDirectProjection(bouncingBall, t, OSI, osnspb, osnspb_pos));
    s->setProjectionMaxIteration(20);
    s->setConstraintTolUnilateral(1e-08);
    s->setConstraintTol(1e-08);
#else
    SP::TimeStepping s(new TimeStepping(bouncingBall, t, OSI, osnspb));
#endif
    s->setNewtonTolerance(1e-10);
    s->setNewtonMaxIteration(10);
    
    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    int N = ceil((T - t0) / h); // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 19;
    SimpleMatrix dataPlot(N + 1, outputSize);

    SP::SiconosVector q = ball->q();
    SP::SiconosVector v = ball->velocity();
    SP::SiconosVector p = ball->p(1);
    SP::SiconosVector qplane = movingplane->q();
    SP::SiconosVector vplane = movingplane->velocity();
    SP::SiconosVector pplane = movingplane->p(1);

    SP::SiconosVector lambda = inter->lambda(1);
    SP::SiconosVector y = inter->y(0);

    SP::SiconosVector reaction = movingplane->reactionToBoundaryConditions();


    dataPlot(0, 0) = bouncingBall->t0();
    dataPlot(0, 1) = (*q)(0);
    dataPlot(0, 2) = (*v)(0);
    dataPlot(0, 3) = (*p)(0);
    dataPlot(0, 4) = (*lambda)(0);
    dataPlot(0, 5) = acos((*q)(3));
    dataPlot(0, 6) = relation0->contactForce()->norm2();

    dataPlot(0, 7) = (*q)(0);
    dataPlot(0, 8) = (*q)(1);
    dataPlot(0, 9) = (*q)(2);
    dataPlot(0, 10) = (*q)(3);
    dataPlot(0, 11) = (*q)(4);
    dataPlot(0, 12) = (*q)(5);
    dataPlot(0, 13) = (*q)(6);

    dataPlot(0, 14) = (*v)(1);
    dataPlot(0, 15) = (*v)(2);

    dataPlot(0, 16) = (*qplane)(2);
    dataPlot(0, 17) = (*vplane)(2);
    dataPlot(0, 18) = (*reaction)(0);

    // --- Time loop ---
    cout << "====> Start computation ... " << endl << endl;
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;
    boost::progress_display show_progress(N);

    boost::timer time;
    time.restart();
    dataPlot(k, 6) = relation0->contactForce()->norm2();
    while (s->hasNextEvent() && k <5000000)
    {
      //      s->computeOneStep();
      s->advanceToEvent();
      // --- Get values to be plotted ---
      dataPlot(k, 0) =  s->nextTime();
      dataPlot(k, 1) = (*q)(0);
      dataPlot(k, 2) = (*v)(0);
      dataPlot(k, 3) = (*p)(0);
      // dataPlot(k, 4) = (*lambda)(0);
      dataPlot(k, 5) = acos((*q)(3));
      dataPlot(k, 6) = relation0->contactForce()->norm2();
      dataPlot(k, 7) = (*q)(0);
      dataPlot(k, 8) = (*q)(1);
      dataPlot(k, 9) = (*q)(2);
      dataPlot(k, 10) = (*q)(3);
      dataPlot(k, 11) = (*q)(4);
      dataPlot(k, 12) = (*q)(5);
      dataPlot(k, 13) = (*q)(6);
      dataPlot(k, 14) = (*v)(1);
      dataPlot(k, 15) = (*v)(2);

      dataPlot(k, 16) = (*qplane)(0);
      dataPlot(k, 17) = (*vplane)(0);

      dataPlot(k, 18) = (*reaction)(0);


      s->nextStep();
      ++show_progress;
      k++;
    }
    cout << endl << "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation Time " << time.elapsed()  << endl;

    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    dataPlot.resize(k, outputSize);
    ioMatrix::write("BallNewtonEulerOnMovingPlane.dat", "ascii", dataPlot, "noDim");

    // Comparison with a reference file
    cout << "====> Comparison with a reference file ..." << endl;
    SimpleMatrix dataPlotRef(dataPlot);
    dataPlotRef.zero();
#ifdef WITH_PROJ
    ioMatrix::read("resultNETS-WITHPROJ.ref", "ascii", dataPlotRef);
#else
    ioMatrix::read("BallNewtonEulerOnMovingPlane.ref", "ascii", dataPlotRef);
#endif
    std::cout << "Error w.r.t reference file = " << (dataPlot - dataPlotRef).normInf() << std::endl;

    if ((dataPlot - dataPlotRef).normInf() > 1e-10)
    {
      std::cout << "Warning. The results is rather different from the reference file. err = " << (dataPlot - dataPlotRef).normInf() << std::endl;
      //(dataPlot-dataPlotRef).display();

      double maxerror = -1e+24;
      unsigned int imax = -1;
      unsigned int jmax = -1;
      double error;
      for (unsigned int ii = 0;  ii < dataPlot.size(0); ii++)
      {
        for (unsigned int jj = 0;  jj < dataPlot.size(1); jj++)
        {
          error = std::abs(dataPlot.getValue(ii, jj) - dataPlotRef.getValue(ii, jj)) ;
          if (error > 1e-12)
          {

            std::cout << "error = " << error << std::endl;
            std::cout << "ii  = " << ii << "  jj  = " << jmax << std::endl;
            std::cout << "dataPlot.getValue(ii,jj) = " <<  dataPlot.getValue(ii, jj) << std::endl;
            std::cout << "dataPlotRef.getValue(ii,jj) =" <<  dataPlotRef.getValue(ii, jj) << std::endl;
          }

          if (error > maxerror)
          {
            maxerror = error ;
            imax = ii;
            jmax = jj;
          }
        }
      }
      std::cout << "max error = " << maxerror << std::endl;
      std::cout << "imax  = " << imax << "  jmax  = " << jmax << std::endl;
      std::cout << "dataPlot.getValue(imax,jmax) = " <<  dataPlot.getValue(imax, jmax) << std::endl;
      std::cout << "dataPlotRef.getValue(imax,jmax) =" <<  dataPlotRef.getValue(imax, jmax) << std::endl;


      return 1;
    }
  }

  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in BouncingBallNETS.cpp" << endl;
  }

}
