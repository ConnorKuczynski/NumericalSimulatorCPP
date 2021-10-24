#include <iostream>
#include <Eigen/Dense>

const int DIM = 2; // simulator dimension
const int X = 0; //enumerate accessing values from pos/vel/acc arrays 
const int Y = 1;
const int Z = 2;

/*
 * stores all the physical parameters of the rocket
 * list of param:
 * dry mass, fuel mass, oxidizer mass, length, width, height, (idea is to model each rocket part as a simple shape where we know the moment of inertia
 *       For example, we seperate the fuel tanks from the dry body, and model them as a sphere (or the best approximate shape), while body is a rectangle)
 * 
 * Summing the center of mass of each individual shape divided by total mass is equal to the center of mass of entire rocket 
 * 
 * we will model the initial rocket as a rectangular prism with a moment of inertia equal to 
 * 
 */
class Tank {
public:
  Tank() {
  }

  /*
   *  we only need to calc z because h_fuel changes as fuel decreases but other dimensions do not change
   */
  void updateCoP(int dltaH)
  {
    com[Z] -= dltaH; 
  }

private:
  Shape s; //simple shape with known area, volume, other parameters, etc
  int mass; 
  int[] cop; // x and y component of cof is constant if shape is symmetric/ z decreases as fuel in tank decreases, relative to P
  
}


class Rocket {
public:
  Rocket() {
    
  }
  /*
   * updatePARAM vs getPARAM (sometimes you just want to get the PARAM vs getting and updating it) this is done to avoid recalculating values that we already know
   *
   */
  double* updateCoM()
  {
    com[X] = dry_mass*length;
    com[Y] = dry_mass*width;
    com[Z] = dry_mass*height;

    for (int i = 0; i < tanks.length; i++)
    {
      for (int j = 0; j < DIM; j++)
      {
        com[j] += tanks[i].mass*tanks[i].cof[j];
      }
    }
    return com;
  }
  /*
   * return Coeffient of drag based on the geometry of the surface of each direction of the rocket
   */
  float* getCD() {
    return cd;
  }
  double* getCOM() {
    return com;
  }
  float* getSA() {
    return sa;
  }
  int tot_mass(){
    return tot_mass;
  }
  
private:
  //constants 

  const int height; // meters
  const int radius; // meters

  const float[] cd; 
  //Shape[] s, eventually we will allow for arbitrary simple shapes but use rectangle for now
  const int P; //gimbal base, all rotations are taken with respect to this point
  const int dry_mass; 
  const float[] sa;

  // variables



  int tot_mass; // Newtons , tot_mass = dry_mass + m_tank1 + m_tank2 + ... + m_tankN
  double[] com = double[DIM];
  
  Tank[] tanks; //arbitary number of fuel tanks (lander should only have two but future designs may have different numbers)

};



float RK4() {


}

class Dynamics {
public:
  Dynamics () {

  }
  /*
   * @param double[] air, model air resistance as a vector 
   * vel is the velocity of the rocket, need to find relative wind vector
   * cd is coffiecient of drag of the rocket in each vector direction / it is proportional to surface area of that part of the rocket
   * TODO account for air resistance that flow around the rocket
   */
  void updateFdrag(double[] air)
  {
    double P = getP(pos[Z]); 
    for (int i = 0; i < DIM; i++)
    {
      fDrag[i] = .5*Math.pow((-vel[i] + air[i]),2)*P*r.getCD()[i]*r.getSA()[i];
    } 
  }
  void updateFg()
  {
    fG[Z] = getG(pos[Z])*r.getMass(); 
  }
  /*
   * given by the control block
   */
  void updateThrust(double[] thrust)
  {
    for (int i = 0; i < DIM; i++)
    {
      fThrust[i] = thrust[i];
    }
  }
  
private:
  double[] ang; //global frame variables
  double[] pos;
  double[] vel;
  double[] acc;
  double[] fDrag = double[DIM]; // just update forces so we don't have to make new memory have time the function is called
  double[] fThrust = double[DIM];
  double[] fG; 
  Rocket* r = new Rocket();
  
};

class Environment {
public:

  const double GM = 3.986004418*Math.pow(10,14); 
  const double R = 5927399.424; //radius of Earth in San Diego (meters)
  double getG(double h)
  {
    return GM/Math.pow((R + h),2);
  }
  /*
   * density (kg/m) of atmosphere as a function of height (meters)
   * https://www.grc.nasa.gov/www/k-12/rocket/atmosmet.html
   */
  double getP(int h)
  {
    if (h < 11000) {
      return 101.29*Math.pow(((15.04-.00649*h + 273.1)/288.08),5.256);
    }
    else if (h > 25000) {
      return 2.488*Math.pow(((-131.21+.00299*h + 273.1)/216.6),-11.388);    }
    else {
      return 22.65*Math.pow(Math.e,(1.73-.000157*h));
    }
  }
  
  Environment() {

  }
};

int main() {
  std::cout << "Hello World!\n";
  Eigen::MatrixXd m;
}