#include "settings/settings.h"
#include <fstream>   // for file operations
#include <iomanip>


void Settings::loadFromFile(std::string filename)
{
  std::string parameterName;
  std::string parameterValue;

  // if(parameterName.find_first_of("/") !=0){
  //  filename.erase(0,filename.find_first_of('/')+1);
  // }

  // open file
  std::ifstream file(filename.c_str(), std::ios::in);


  
  // check if file is open
  if (!file.is_open())
  {
    std::cout << "Could not open parameter file \"" << filename << "\"." << std::endl;
    return;
  }

  // loop over lines of file
  for (int lineNo = 0;; lineNo++)
  {
    // read line
    std::string line;
    getline(file, line);

    // at the end of the file break for loop
    if (file.eof()){
      file.close();
      break;
    }
    //skip lines that start with #
    if (line[0]=='#')
    {
      continue;
    }

    //skip lines with no = in them
    if (line.find('=')==std::string::npos)
    {
      continue;
    }
    
    //erase blank spaces at the beginning
    while (line.find_first_of("\t")==0||line.find_first_of(" ")==0)
    {
      line.erase(0,1);
    }
    
    parameterName = line.substr(0,line.find_first_of('='));
    //remove spaces
    if (parameterName.find_first_of(" \t") != std::string::npos)
    {
      parameterName.erase(parameterName.find_first_of(" \t"));
    }
    
    
    parameterValue = line.substr(line.find_first_of('=')+1);
    //remove spaces
    while(parameterValue.find(' ')!= std::string::npos)
    {
      parameterValue.erase(parameterValue.find_first_of(' '),1);
    }

    //remove comments at the end of the line
    if(parameterValue.find('#')!= std::string::npos)
    {
      parameterValue.erase(parameterValue.find('#'));
    }

    //find corresponding parameter name in settings struct and set the correct value
    if (parameterName=="endTime")
    {
      endTime = atof(parameterValue.c_str());
    }

    if (parameterName=="physicalSizeStart")
    {
      physicalSize[0] = atof(parameterValue.c_str());
    }

    if (parameterName=="physicalSizeEnd")
    {
      physicalSize[1] = atof(parameterValue.c_str());
    }

    if (parameterName=="dirichletBottomX")
    {
     dirichletBcBottom[0] = atof(parameterValue.c_str());  
    }

    if (parameterName=="dirichletBottomY")
    {
      dirichletBcBottom[1] = atof(parameterValue.c_str());  
    }    

    if (parameterName=="dirichletTopX")
    {
        dirichletBcTop[0] = atof(parameterValue.c_str());  
    }

    if (parameterName=="dirichletTopY")
    {
        dirichletBcTop[1] = atof(parameterValue.c_str());      
    }

    if (parameterName=="dirichletLeftX")
    {
      dirichletBcLeft[0] = atof(parameterValue.c_str());
    }

    if (parameterName=="dirichletLeftY")
    {
        dirichletBcLeft[1] = atof(parameterValue.c_str());
    }

    if (parameterName=="dirichletRightX")
    {
       dirichletBcRight[0] = atof(parameterValue.c_str()); 
    }

    if (parameterName=="dirichletRightY")
    {
       dirichletBcRight[1] = atof(parameterValue.c_str()); 
    }

    if (parameterName=="nCells")
    {
         nCells[0] = atof(parameterValue.c_str());
    }

    if (parameterName=="maximumDt")
    {
        dt = atof(parameterValue.c_str());
    }

    if (parameterName=="CFL")
    {
        CFL = atof(parameterValue.c_str());
    }

    if (parameterName=="maximumNumberOfIterations")
    {
       maximumNumberOfIterations = atof(parameterValue.c_str());
    }
    if (parameterName=="PP_N")
    {
       PP_N = atof(parameterValue.c_str());
    }
    if (parameterName=="fluxFunction"){
      fluxFunction=parameterValue;
    }
    if (parameterName=="initialCondition"){
      initialCondition=parameterValue;
    }

    if (parameterName=="OutputName"){
      OutputName=parameterValue;
    }

    if (parameterName=="RiemannSolver"){
      RiemannSolver=parameterValue;
    }

    if(parameterName=="initCondA"){
      initCondA = atof(parameterValue.c_str());
    }
    if(parameterName=="initCondB"){
      initCondB = atof(parameterValue.c_str());
    }
    if(parameterName=="BarenblattM"){
      BarenblattM = atoi(parameterValue.c_str());
    }
    if(parameterName=="BarenblattTime"){
      BarenblattTime = atof(parameterValue.c_str());
    }
  }
}

//output all settings to console
void Settings::printSettings()
{
  std::cout << "Settings: " << std::endl
    <<"PP_N: "<<PP_N<<std::endl
    << "  Interval: [" << physicalSize[0] << " , " << physicalSize[1] << "] , nCells: " << nCells[0]  << std::endl
    << "  endTime: " << endTime <<  ", maximum dt: " << dt << ", maximumNumberOfIterations: " << maximumNumberOfIterations << std::endl;
}