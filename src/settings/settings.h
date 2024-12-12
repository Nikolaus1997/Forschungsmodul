#pragma once

#include <iostream>
#include <array>

/** All settings that parametrize a simulation run.
 */
struct Settings
{
  std::string OutputName = "linearFlux";
  std::array<int,1> nCells;             //< number of cells in x and y direction
  std::array<double,2> physicalSize;    //< physical size of the domain
  double endTime = 10.0;                //< end time of the simulation
  double dt = 0.1;                      //< maximum time step width
  double CFL = 0.9;                     //< CFL number
  int PP_N = 5;                         //< Polynomial degree of the basis functions

  std::string fluxFunction = "linear";        //< Flux function 
  std::string initialCondition = "unitStep";  //< Initial condition
  std::string RiemannSolver = "upwind";  //< Initial condition
  double initCondA = 0.0;                     //< Starting point of the initial Condition
  double initCondB = 0.0;                     //< End point of the initial Condition

  int BarenblattM = 0;
  double BarenblattTime= 1.0;

  std::array<double,2> dirichletBcBottom;  //< prescribed values of u,v at bottom of domain
  std::array<double,2> dirichletBcTop;     //< prescribed values of u,v at top of domain
  std::array<double,2> dirichletBcLeft;    //< prescribed values of u,v at left of domain
  std::array<double,2> dirichletBcRight;   //< prescribed values of u,v at right of domain

  int maximumNumberOfIterations = 1e5;    //< maximum number of iterations in the solver

  //! parse a text file with settings, each line contains "<parameterName> = <value>"
  void loadFromFile(std::string filename);

  //! output all settings to console
  void printSettings();
};