
#include <iostream>
#include <cmath>
#include "main.hpp"

//#include "raylib.h"
using namespace std;




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
 * Momentum of inertia will be stored relative to the P of the entire rocket
 */
class Environment {
public:
  const double GM = 3.986004418*pow(10,14); 
  const double R = 6371927; //radius of Earth in San Diego (meters)
  double* air;
  Environment(double* a) {
    this->air = a;
  }
  double* getAir()
  {
    return air;
  }
  double getG(double h)
  {
    return GM/pow((R + h),2);
  }
  /*
   * density (kg/m) of atmosphere as a function of height (meters) above sea level
   * https://www.grc.nasa.gov/www/k-12/rocket/atmosmet.html
   * WARNING, this pressure does not scale for h < 0, obviously because that means it crashed into the ocean
   * this leads to the rocket not reaching terminal velocity due to the factors of gravity increasing and pressure increasing
   * when h < 0
   */
  double getP(int h)
  {
    if (h < 11000) {
      return (101.29*pow(((15.04-.00649*h + 273.1)/288.08),5.256))/(.2869*((15.04-.00649*h + 273.1)));
    }
    else if (h > 25000) {
      return 2.488*pow(((-131.21+.00299*h + 273.1)/216.6),-11.388)/(.2869*((-131.21+.00299*h + 273.1)));   
    }
    else {
      return 22.65*pow(E,(1.73-.000157*h))/(.2869*((-56.46 + 273.1)));
    }
  }
};
//all coordinate by default are measured relative to P of rocket
//to specify local coordinate system use suffix Local

//Shape 

float* Shape::getSA() {
  return sa; 
}
double* Shape::getI() {
  return I;
}
float* Shape::getOrigin() {
  return P;
}
double* Shape::updateCOVLocal() {
  for (int i = 0; i < DIM; i++) {
    covLocal[i] = cov[i] - P[i];
  } 
  return covLocal;
}
double* Shape::getCOVLocal() {
  return covLocal;
}
double* Shape::getCOV() {
  return cov;
}
double* Shape::updateI(float) { throw std::invalid_argument("Shape updateI() was run instead of a subclass shape"); }
void Shape::updateCOM(float*, float, float, float) { throw std::invalid_argument("Shape updateCOM() was run instead of a subclass shape"); }
void Shape::updateFilledVolume(Tank* t, double dltaM) { throw std::invalid_argument("Shape updateFilledVolume() was run instead of a subclass shape"); };
double* Shape::updateIFuel(Tank* t) {  throw std::invalid_argument("Shape updateIFuel() was run instead of a subclass shape"); }
  //  virtual double* Shape::updateI(float) { throw std::invalid_argument("Shape updateI() was run instead of a subclass shape");};
   // virtual void Shape::updateCOM(float*, float, float, float) { throw std::invalid_argument("Shape updateCOM() was run instead of a subclass shape"); };
   // virtual void Shape::updateFilledVolume(Tank* t, double dltaM);
   // virtual double* Shape::updateIFuel(Tank* t);
    

/*
 * length = X, width = Y, height = Z
 */

//END Shape

//RectangularPrism

int RectagularPrism::getAreaLW()
{
  return length*width;
}
int RectagularPrism::getLength()
{
  return length;
}
int RectagularPrism::getWidth()
{
  return width;
}
int RectagularPrism::getHeight()
{
  return height;
}
void RectagularPrism::updateFilledVolume(Tank* t, double dltaM) {
  t->updateFilledVolumeHelper(this, dltaM);
}
double* RectagularPrism::updateIFuel(Tank* t) {
  return t->updateIFuelHelper(this);
}
double* RectagularPrism::updateI(float mass)
{
  double* r = cov;

  I[X] = mass*((1/12.)*(pow(width,2) + pow(height,2)) + pow(r[Y],2) + pow(r[Z],2));
  I[Y] = mass*((1/12.)*(pow(length,2) + pow(height,2)) + pow(r[X],2) + pow(r[Z],2));
  I[Z] = mass*((1/12.)*(pow(length,2) + pow(width,2)) + pow(r[X],2) + pow(r[Y],2));

  return I;
}
void RectagularPrism::updateCOM(float* com, float fuelH, float fuelM, float dry_mass) 
{
  for (int i = 0; i < DIM; i++) {
    com[i] = 0;
    if (i == Z) { 
      com[Z] += ((fuelH/2 + P[Z])*fuelM)/(fuelM + dry_mass); 
      com[Z] += ((height/2. + P[Z])*dry_mass)/(fuelM + dry_mass);
    }
    else {
      com[i] += ((P[i])*fuelM)/(fuelM + dry_mass);
      com[i] += ((P[i])*dry_mass)/(fuelM + dry_mass);
    }
    
  }
  
}

//End RectangularPrism

//Sphere
float Sphere::getRadius() {
  return radius;
}
void Sphere::updateFilledVolume(Tank* t, double dltaM) {
  t->updateFilledVolumeHelper(this, dltaM);
}
float Sphere::getVolume() {
  return (4/3.)*PI*pow(radius,3);
}
double* Sphere::updateIFuel(Tank* t) {
  return t->updateIFuelHelper(this);
}
//hollow sphere
double* Sphere::updateI(float mass) {
  double* r = cov;
  I[X] = mass*((2/3.)*pow(radius,2) + pow(r[Y],2) + pow(r[Z],2));
  I[Y] = mass*((2/3.)*pow(radius,2) + pow(r[X],2) + pow(r[Z],2));
  I[Z] = mass*((2/3.)*pow(radius,2) + pow(r[Y],2) + pow(r[X],2));

  return I;
}
void Sphere::updateCOM(float* com, float fuelH, float fuelM, float dry_mass) 
{
  for (int i = 0; i < DIM; i++) {
    com[i] = 0;
    if (i == Z) { 
      float ang = 2*acos(1-fuelH/radius);
      com[Z] += ((radius - 4*radius*pow(sin(ang/2),3)/(3*(ang-sin(ang))) + P[Z])*fuelM)/(fuelM+dry_mass);
      com[Z] += ((radius/2. + P[Z])*dry_mass)/(fuelM + dry_mass);
    }
    else {
      com[i] += ((P[i] + radius)*dry_mass)/(fuelM + dry_mass);
      com[i] += ((P[i])*fuelH)/(fuelM + dry_mass);
    }
    
  }
  
}

//End Sphere

//Tank
/*
* currently we are assuming that the fuel stays perpendicular to the z axis of the rocket which obviously will need to changed in the future.
*/
//calculate momentum of inertia caused by fuel, calculate base on container of fuel
double* Tank::updateIFuel(Shape* s) {
  return (s->updateIFuel(this));
}
double* Tank::updateIFuelHelper(RectagularPrism* s)
{
  float Z_coor = com[Z];
  double* r = s->getCOV();
  int w = s->getWidth();
  int l = s->getLength();
  fuelI[X] = fuelM*((1/12.)*(pow(w,2) + pow(fuelH,2)) + pow(r[Y],2) + pow(Z_coor,2));
  fuelI[Y] = fuelM*((1/12.)*(pow(l,2) + pow(fuelH,2)) + pow(r[X],2) + pow(Z_coor,2));
  fuelI[Z] = fuelM*((1/12.)*(pow(w,2) + pow(l,2)) + pow(r[X],2) + pow(r[Y],2));

  return fuelI;
}
//http://blitiri.blogspot.com/2014/05/mass-moment-of-inertia-of-hemisphere.html
double* Tank::updateIFuelHelper(Sphere* s)
{
  int radius = s->getRadius();
  float* origin = s->getOrigin();
  float* r = com;
  float vol = s->getVolume();
  //calculate X and Y is same as finding filled sphere but find ratio of filled to full and then use parallel axis to account for center of mass not being center of volume
  //next shift to pivot point of the rocket. For Z direction solve integral for finding sphere but instead of -R to R calculate 0 to fuelH.
  //Account for difference of center of volume of sphere and center of mass and shift to axis of rocket pivot
  fuelI[X] = fuelM*((volumeFilled/vol)*(2/5.)*pow(radius,2) + pow(r[Y]-origin[Y],2) + pow(r[Z]-origin[Z],2) + (pow(r[Y],2) + pow(r[Z],2)));
  fuelI[Y] = fuelM*((volumeFilled/vol)*(2/5.)*pow(radius,2) + pow(r[X]-origin[X],2) + pow(r[Z]-origin[Z],2) + pow(r[X],2) + pow(r[Z],2));
  fuelI[Z] = fuelM*(PI/volumeFilled*(pow(fuelH,5)/5 - (2/3.)*pow(fuelH,3)*pow(radius,2) + fuelH*pow(radius,4)) + pow(r[Y]-origin[Y],2) + pow(r[X]-origin[X],2) + pow(r[Y],2) + pow(r[X],2));

  return fuelI;
}
//change in mass is in kg 
int Tank::updateMass(double dt)
{
  if (fuelM < 0) { return dry_mass; }
  double dltaM = this->mass_vel_out*dt;
  updateFilledVolume(s, dltaM);
  s->updateCOM(com, fuelH, fuelM, dry_mass);
  updateIFuel(s);
  fuelM -= dltaM;

  return (dry_mass+fuelM);
}
float Tank::getDryMass()
{
  return dry_mass;
}
float Tank::getFuelMass()
{
  return fuelM;
}
float Tank::getMass()
{
  return dry_mass + fuelM;
}
float* Tank::getCOM()
{
  s->updateCOM(com, fuelH, fuelM, dry_mass);
  return com;
}
Shape* Tank::getShape() {
  return s;
} 
void Tank::updateFilledVolume(Shape* s, double dltaM) {
  s->updateFilledVolume(this, dltaM);
}
void Tank::updateFilledVolumeHelper(RectagularPrism* s, double dltaM) {
  float dltaV = dltaM*densityFuel;
  float dltaH = dltaV/s->getAreaLW();
  volumeFilled -= dltaV;
  fuelH -= dltaH;
}
double* Tank::getIFuel() {
  return fuelI;
}
float Tank::getFuelHeight(){
  return fuelH;
}
void Tank::updateFilledVolumeHelper(Sphere* s, double dltaM) {
  float dltaV = dltaM*densityFuel;
  float dltaH = dltaV/(PI*fuelH*(2*s->getRadius()-fuelH));
  volumeFilled -= dltaV;
  fuelH -= dltaH;
}

    
    

    
    
//End Tank



double* Control::updateThrust(Dynamics::State* s, Rocket* r, float fuelMass, double dt) {
  //cout << "TIME ELAPSED: " << s->timeElapsed << " " << time_itr << endl;
  
  updateProfile(dt); 
  
  if (fuelMass < 0) {
    for (int i = 0; i < DIM; i++) {
      thrust[i] = 0;
    }
    return thrust;
  }
  else {
    for (int i = 0; i < DIM; i++) {
      thrust[i] = flight_profile[time_itr][i];
    }
  }
  //cout << time_itr << "/" << (int)(TIME_FINAL/SECS_PER_ITR) << endl;
  /*
  if (flight_profile[time_itr] == HOVER)
  {
    thrust[Z] += -s->acc[Z]*s->mass;
  }
  else if (flight_profile[time_itr] == MAX_THRUST) {
    thrust[Z] = r->getMaxThrust();
  }
  else {
    thrust[Z] = 0;
    cout << "BAD NEWS";
  }
  */
  return thrust;
}
double** Control::initProfile() {

  flight_profile = new double*[(int)(TIME_FINAL/SECS_PER_ITR)];
  for (int i = 0; i < (int)(TIME_FINAL/SECS_PER_ITR); i++) {
    flight_profile[i] = new double[DIM];
  }
  
  float x_ang = (PI/180)*5; //radians
  float y_ang = (PI/180)*5; 
  float x = sin(x_ang); 
  float y = sin(y_ang); 
  float z = sqrt(1 - pow(x,2) - pow(y,2)); //x^2 + y^2 + z^2 = 1
  for (int i = 0; i < 100; i++) {
    flight_profile[i][X] = MAX_THRUST*x;
    flight_profile[i][Y] = MAX_THRUST*y;
    flight_profile[i][Z] = MAX_THRUST*z;
  }
  
  for (int i = 100; i < TIME_FINAL/SECS_PER_ITR; i++) {
      flight_profile[i][X] = 0;
      flight_profile[i][Y] = 0;
      flight_profile[i][Z] = MAX_THRUST;
  }
  return flight_profile;
}
void Control::updateProfile(double dt) {
  if (secs >= SECS_PER_ITR) 
  { 
    time_itr++;
    secs = dt; 
  }
  else { secs += dt; }
}
double* Control::getThrust() {
  return thrust;
}


//dry mass does not include dry mass of tanks
    /*
    * updatePARAM vs getPARAM (sometimes you just want to get the PARAM vs getting and updating it) this is done to avoid recalculating values that we already know
    *
    */
void Rocket::updateComDry()
{
  for (int i = 0; i < DIM; i++) {

    comDry[i] = cov[i]*dry_mass;
  }
}
Tank** Rocket::getTanks()
{
  return tanks;
}
double* Rocket::updateIDry()
{
  for (int j = 0; j < numShapes; j++)
  {
    //cout << "Shapes (X,Y,Z)" << endl;
    double* IShapes = s[j]->updateI(massParts[j]); 
    //cout << IShapes[X] << " " << IShapes[Y] << " " << IShapes[Z] << endl;
    //cout << "Dry mass tanks" << endl;
    double* ITanksDry = tanks[j]->getShape()->updateI(tanks[j]->getDryMass());
    //cout << ITanksDry[X] << " " << ITanksDry[Y] << " " << ITanksDry[Z] << endl;
    //cout << "Fuel height" << tanks[j]->getFuelHeight() << endl; 
    for (int i = 0; i < DIM; i++)
    {
      IDry[i] = IShapes[i];
      IDry[i] += ITanksDry[i];
    }
  }
  return IDry;
}
double* Rocket::updateI() {     
  for (int j = 0; j < numTanks; j++)
  {
    
    double* ITanks = tanks[j]->getIFuel();
    for (int i = 0; i < DIM; i++)
    {
      I[i] = IDry[i];
      I[i] += ITanks[i];  
    }
    //cout << I[X] << " " << I[Y] <<  " " << I[Z] << endl;
  }
  
  
  return I;
}
double* Rocket::updateCOM()
{
  
  for (int j = 0; j < DIM; j++)
  {
    com[j] = comDry[j];
  }
  
  for (int i = 0; i < numTanks; i++)
  {
    for (int j = 0; j < DIM; j++)
    {
      com[j] += (tanks[i]->getMass()*(tanks[i]->getCOM()[j]));
      com[j] /= tot_mass; 
    }    
  }
  return com;
}
double* Rocket::updateCOP()
{
  return cop;
}
double* Rocket::getI()
{
  return I;
}
int Rocket::updateMass(double dt)
{
  tot_mass = dry_mass;
  for(int i = 0; i < numTanks; i++)
  {
    Tank* tank = tanks[i];
    tot_mass += tank->updateMass(dt);
  }
  updateI();
  return tot_mass;
}
/*
* return Coeffient of drag based on the geometry of the surface of each direction of the rocket
*/
float* Rocket::getCD() {
  return cd;
}
double* Rocket::getCOM() {
  return com;
}
double* Rocket::getCOP() {
  return cop;
}
float* Rocket::getSA() {
  return sa;
}
int Rocket::getMass(){
  return tot_mass;
}
const int* Rocket::getP() {
  return P;
}
const int Rocket::getNumParts() {
  return numShapes;
}
const int Rocket::getNumTanks() {
  return numTanks;
}
const double Rocket::getMaxThrust() {
  return max_thrust;
}


 /* 
void RK4(double[] initialVec, double[] vec, double (*diffEq)(float, double, double), float dt) {
    int[] n = int[DIM];
 
    float[] k1, k2, k3, k4, k5;
 
    for (int i = 0; i < DIM; i++)
    {
        // Apply Runge Kutta Formulas to find
        // next value of y
        k1[i] = //dt*diffEq(0, initialVec);
        k2[i] = dt*diffEq(.5*dt, initialVec + 0.5*k1[i], );
        k3[i] = dt*diffEq(.5*dt, initialVec + 0.5*k2[i]);
        k4[i] = dt*diffEq(dt, initialVec + k3[i]);
 
        vec[i] += (1.0/6.0)*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
    }
} */


/*
 * ORDER of calculating values
 * ALWAYS update F before M 
 * 
 *  updateFdrag
 *  updateFg
 *  updateFThrust
 * 
 *  updateNetF
 * 
 *  updateMdrag
 *  updateMg
 *  updateMThrust
 * 
 *  updateNetM
 *  updateAcc
 * 
 *  updateNetM
 *  updateAngAcc
 *  
 * 
 *  clearNet
 */


//TODO updatePrev copies current state pointer to prev state which is just the same memory so need to fix that eventually
void Dynamics::State::updatePrev(State* curr)
{
  ang = curr->ang; //global frame variables 
  pos = curr->pos;
  vel = curr->vel;
  timeElapsed = curr->timeElapsed;
  ang_vel = curr->ang_vel;
  ang_acc = curr->ang_acc;
  acc = curr->acc;
  mass = curr->mass;
}
string Dynamics::State::to_string() 
{
  string output = "(X,Y,Z) \n";
  string posStr = "POSITION: (";
  string velStr = "VELOCITY: (";
  string accStr = "ACCERLATION: (";
  string angStr = "ANGLE: (";
  string ang_velStr = "ANGULAR VELOCITY: (";
  string ang_accStr = "ANGULAR ACCERLATION: (";
  for (int i = 0; i < DIM; i++) {
    if (i < DIM - 1) { 
      posStr += (std::to_string(pos[i]) + ","); 
      velStr += (std::to_string(vel[i]) + ","); 
      accStr += (std::to_string(acc[i]) + ","); 
      angStr += (std::to_string((int)(ang[i]*(180/PI)) % 360) + ","); 
      ang_velStr += (std::to_string((int)(ang_vel[i]*(180/PI))) + ","); 
      ang_accStr += (std::to_string((int)(ang_acc[i]*(180/PI))) + ","); 
    }
    else { 
      posStr += (std::to_string(pos[i]) + ")\n"); 
      velStr += (std::to_string(vel[i]) + ")\n"); 
      accStr += (std::to_string(acc[i]) + ")\n"); 
      angStr += (std::to_string((int)(ang[i]*(180/PI)) % 360) + ")\n"); 
      ang_velStr += (std::to_string((int)(ang_vel[i]*(180/PI))) + ")\n"); 
      ang_accStr += (std::to_string((int)(ang_acc[i]*(180/PI))) + ")\n"); 
    }
  } 
  output = posStr + velStr + accStr + angStr + ang_velStr + ang_accStr;
  return output;
}
  string Dynamics::to_string()
  {
    string output = s->to_string();
    string Istr = "MOMENT OF INERTIA: (";
    string COMstr = "CENTER OF MASS: (";
    string massStr = "Mass: (" + std::to_string(r->getMass()) + ")\n";
    for (int i = 0; i < DIM; i++) {
      if (i < DIM - 1) { 
        Istr += (std::to_string(r->getI()[i]) + ","); 
        COMstr += (std::to_string(r->getCOM()[i]) + ","); 
      }
      else { 
        Istr += (std::to_string(r->getI()[i]) + ")\n"); 
        COMstr += (std::to_string(r->getCOM()[i]) + ")\n"); 
      }
    } 
    output += (Istr + COMstr + massStr);
    return output;
  }
  void Dynamics::RK4(double dt)
  {
    int n[DIM];
  
    double k1Pos[DIM], k2Pos[DIM], k3Pos[DIM], k4Pos[DIM];
    double k1Ang[DIM], k2Ang[DIM], k3Ang[DIM], k4Ang[DIM];  

//initial force calculations
    Dynamics* state0 = new Dynamics(*this); //deep clone current state
    double* forces[3] = {mDrag,mThrust,mG};

    state0->updateInerMom();
    state0->updateFdrag();
    state0->updateFg();
   // state0->updateFThrust();
    
    state0->updateNetF(forces);
    state0->updateAcc();

    state0->updateMdrag();
    state0->updateMg();
    state0->updateMThrust();
    state0->updateNetM(forces);
    state0->updateAngAcc();

    for (int i = 0; i < DIM; i++) {
      k1Pos[i] = dt*dPosdt(0*dt, state0->getVel()[i], state0->getAcc()[i]); //dltaPos
      k1Ang[i] = dt*dPosdt(0*dt, state0->getAngVel()[i], state0->getAngAcc()[i]); //dltaAng
    }
    

    //this order is important we update position based off .5*dltaPos
    //then we update vel because its already applied from dPos/dt
    state0->updatePos(.5, k1Pos);
    state0->updateVel(.5*dt);
    state0->updateAng(.5, k1Ang);
    state0->updateAngVel(.5*dt);

    state0->r->updateMass(.5*dt);
    state0->r->updateCOM();
    state0->updateInerMom();
    state0->updateFdrag();
    state0->updateFg();
    state0->updateNetF(forces);
    state0->updateAcc();

    state0->updateMdrag();
    state0->updateMg();
    state0->updateMThrust();
    state0->updateNetM(forces);
    state0->updateAngAcc();
    state0->clearNet();

    for (int i = 0; i < DIM; i++) {
      k2Pos[i] = dt*dPosdt(.5*dt, getVel()[i], getAcc()[i]);
      k2Ang[i] = dt*dPosdt(.5*dt, getAngVel()[i], getAngAcc()[i]);
    }

    state0->updatePos(.5, k2Pos);
    state0->updateVel(.5*dt);
    state0->updateAng(.5, k2Ang);
    state0->updateAngVel(.5*dt);

    state0->r->updateMass(.5*dt);
    state0->r->updateCOM();
    state0->updateInerMom();
    state0->updateFdrag();
    state0->updateFg();
    state0->updateNetF(forces);
    state0->updateAcc();
    
    state0->updateMdrag();
    state0->updateMg();
    state0->updateMThrust();
    state0->updateNetM(forces);
    state0->updateAngAcc();
    state0->clearNet();

    for (int i = 0; i < DIM; i++) {
      k3Pos[i] = dt*dPosdt(.5*dt, getVel()[i], getAcc()[i]);
      k3Ang[i] = dt*dPosdt(.5*dt, getAngVel()[i], getAngAcc()[i]);
    }

    state0->updatePos(dt, k3Pos);
    state0->updateVel(dt);
    state0->updateAng(dt, k3Ang);
    state0->updateAngVel(dt);

    state0->r->updateMass(dt);
    state0->r->updateCOM();
    state0->updateInerMom();
    state0->updateFdrag();
    state0->updateFg();
    state0->updateNetF(forces);
    state0->updateAcc();

    state0->updateMdrag();
    state0->updateMg();
    state0->updateMThrust();
    state0->updateNetM(forces);
    state0->updateAngAcc();
    state0->clearNet();

    for (int i = 0; i < DIM; i++) {
      k4Pos[i] = dt*dPosdt(dt, getVel()[i], getAcc()[i]);
      k4Ang[i] = dt*dPosdt(dt, getAngVel()[i], getAngAcc()[i]);

      this->s->pos[i] += (1.0/6.0)*(k1Pos[i] + 2*k2Pos[i] + 2*k3Pos[i] + k4Pos[i]);
      this->s->ang[i] += (1.0/6.0)*(k1Ang[i] + 2*k2Ang[i] + 2*k3Ang[i] + k4Ang[i]);
    }

    
  }
  void Dynamics::updateState(double dt) {
    double* forces[3] = {fDrag,fThrust,fG};
    double* moments[3] = {mDrag,mThrust,mG};
    

    this->updateFdrag();
    if (abs(fDrag[X]) > 1000000) {
      cout << "STOP" << endl;
    } 
    this->updateFg();
    if (abs(fG[X]) > 100000) {
      cout << "STOP" << endl;
    } 
    this->updateFThrust(s->prev, dt);
    if (abs(fThrust[X]) > 100000) {
      cout << "STOP" << endl;
    } 
    this->updateNetF(forces);
    this->updateAcc();

    this->updateMdrag();
    this->updateMg();
    this->updateMThrust();
    this->updateInerMom();
    this->updateNetM(moments);
    this->updateAngAcc();
    s->timeElapsed += dt;
  }
  void Dynamics::euler(double dt) {
    this->updateState(dt);

    for (int i = 0; i < DIM; i++)
    {
      s->pos[i] += dPosdt(dt, s->vel[i], s->acc[i]);
      s->vel[i] += s->acc[i]*dt;
    }
    for (int i = 0; i < DIM; i++)
    {
      s->ang[i] += dPosdt(dt, s->ang_vel[i], s->ang_acc[i]);
      s->ang_vel[i] += s->ang_acc[i]*dt;
    }
    s->mass = this->r->updateMass(dt);
    this->r->updateCOM();

    s->prev->updatePrev(s);
    clearNet();
  }

  double Dynamics::dPosdt(double dt, double vel, double acc)
  {
    return .5*acc*pow(dt,2) + vel*dt;
  }
  /*
   * this may get replaced/moved in the future but this succesfully calculates the cross_p, assumes DIM = 3
   */ 
  void Dynamics::crossProduct(double dist[], double force[], double moment[])
  {
    moment[0] = dist[1] * force[2] - dist[2] * force[1];
    moment[1] = dist[2] * force[0] - dist[0] * force[2];
    moment[2] = dist[0] * force[1] - dist[1] * force[0];
  }

  /*
   * @param double[] air, model air resistance as a vector 
   * vel is the velocity of the rocket, need to find relative wind vector
   * cd is coffiecient of drag of the rocket in each vector direction / it is proportional to surface area of that part of the rocket
   * TODO account for air resistance that flow around the rocket
   */
  void Dynamics::updateFdrag()
  {
    double altitude = s->pos[Z];
    double P = e->getP(altitude);
    double vel_relativeTo_Air_Glo[DIM];
    double vel_relativeTo_Air_Roc[DIM];
    double fDrag_Roc[DIM];
    //double fDrag_Ang_Roc[DIM];
    for (int i = 0; i < DIM; i++) {
      vel_relativeTo_Air_Glo[i] = (-s->vel[i] + e->air[i]); 
    } 

    globalToRocketFrame(vel_relativeTo_Air_Roc, vel_relativeTo_Air_Glo);
  
    for (int i = 0; i < DIM; i++)
    {
      if (vel_relativeTo_Air_Roc[i] > 0) {
        fDrag_Roc[i] = .5*pow(vel_relativeTo_Air_Roc[i],2)*P*r->getCD()[i]*r->getSA()[i];
      }
      else {
        fDrag_Roc[i] = -.5*pow(vel_relativeTo_Air_Roc[i],2)*P*r->getCD()[i]*r->getSA()[i];
      }
      //angular drag
    } 
    /*
    if (s->ang_vel[X] > 0) {
        fDrag_Ang_Roc[X] = .5*pow(s->ang_vel[X],2)*P*((r->getCD()[Y]*r->getSA()[Y])+(r->getCD()[Z]*r->getSA()[Z]));
    }
    else {
        fDrag_Ang_Roc[X] = -.5*pow(s->ang_vel[X],2)*P*((r->getCD()[Y]*r->getSA()[Y])+(r->getCD()[Z]*r->getSA()[Z]));
    }
    if (s->ang_vel[Y] > 0) {
        fDrag_Ang_Roc[Y] = .5*pow(s->ang_vel[Y],2)*P*((r->getCD()[X]*r->getSA()[X])+(r->getCD()[Z]*r->getSA()[Z]));
    }
    else {
        fDrag_Ang_Roc[Y] = -.5*pow(s->ang_vel[Y],2)*P*((r->getCD()[X]*r->getSA()[X])+(r->getCD()[Z]*r->getSA()[Z]));
    }
    if (s->ang_vel[Z] > 0) {
        fDrag_Ang_Roc[Z] = .5*pow(s->ang_vel[Z],2)*P*((r->getCD()[X]*r->getSA()[X])+(r->getCD()[Y]*r->getSA()[Y]));
    }
    else {
        fDrag_Ang_Roc[Z] = -.5*pow(s->ang_vel[Z],2)*P*((r->getCD()[X]*r->getSA()[X])+(r->getCD()[Y]*r->getSA()[Y]));
    }
    */
    rocketToGlobalFrame(fDrag_Roc, fDrag);
    //rocketToGlobalFrame(fDrag_Ang_Roc, fDrag_Ang);
  }
  /*
   * the moment created from the force of drag
   */
  void Dynamics::updateMdrag()
  {
    double dist[DIM];
    for (int i = 0; i < DIM; i++) {
      dist[i] = r->getCOP()[i];
    }
    /*
    double fDrag_Tot[DIM];
    for (int i = 0; i < DIM; i++) {
      fDrag_Tot[i] = fDrag[i]; //+ fDrag_Ang[i];
    }
    */
    double fDrag_Roc[DIM];
    globalToRocketFrame(fDrag_Roc, fDrag);
    crossProduct(dist, fDrag_Roc, mDrag);
  }


//set getG parameters to 0 for constant fG
  void Dynamics::updateFg()
  {
    //global plane angles and forces
    double altitude = s->pos[Z];

    double fG_Z = -e->getG(altitude)*r->getMass();

    fG[X] = 0;
    fG[Y] = 0;
    fG[Z] = fG_Z; 

  }
  void Dynamics::updateMg()
  {
    //convert from global to rocket frame
    double a = s->ang[X];
    double b = s->ang[Y];
    double y = s->ang[Z];

    double fG_roc[DIM];

    //IMPORTANT: if we do decide to change transformation matrix remember to change this one
    //also because this assumes fG[X] and fG[Y] are 0 to simpify calculations compared to globalToRocFrame function

    fG_roc[X] = (cos(a)*sin(b)*cos(y)+sin(a)*sin(y))*fG[Z];
    fG_roc[Y] = (sin(a)*sin(b)*cos(y)-cos(a)*sin(y))*fG[Z];
    fG_roc[Z] = cos(b)*cos(y)*fG[Z]; 

    crossProduct(r->getCOM(), fG_roc, mG);
  }

   /*
   * given by the control block
   */
  void Dynamics::updateFThrust(Dynamics::State* s, double dt)
  {
    float fuelMass = 0;
    for (int i = 0; i < r->getNumTanks(); i++)
    {
      fuelMass += r->getTanks()[i]->getFuelMass();
    }
    
    double* thrust_roc = c->updateThrust(s, r, fuelMass, dt);

    rocketToGlobalFrame(thrust_roc, fThrust);

  }
  
  void Dynamics::updateMThrust()
  {
    //c->getThrust() is thrust from rocket plane
    crossProduct(r->getCOM(), c->getThrust(), mThrust);
  }

  //moment caused by movement
  void Dynamics::updateInerMom()
  {
    double netF_roc[DIM];

    globalToRocketFrame(netF_roc, netF);

    crossProduct(r->getCOM(), netF_roc, inerM);
  }
  void Dynamics::rocketToGlobalFrame(double* roc, double* glo) {
    double a = s->ang[X];
    double b = s->ang[Y];
    double y = s->ang[Z];

    //use transpose matrix to go from rocket to global plane

    glo[X] = (cos(a)*cos(b)*roc[X]) + (sin(a)*cos(b))*roc[Y] + (-sin(b))*roc[Z];
    glo[Y] = (cos(a)*sin(b)*sin(y)-sin(a)*cos(y))*roc[X] + (sin(a)*sin(b)*sin(y)+cos(a)*cos(y))*roc[Y] + (cos(b)*sin(y))*roc[Z];
    glo[Z] = (cos(a)*sin(b)*cos(y)+sin(a)*sin(y))*roc[X] + (sin(a)*sin(b)*cos(y)-cos(a)*sin(y))*roc[Y] + (cos(b)*cos(y))*roc[Z];
  }
  void Dynamics::globalToRocketFrame(double* roc, double* glo) {
    double a = s->ang[X];
    double b = s->ang[Y];
    double y = s->ang[Z];
    
    roc[X] = cos(a)*cos(b)*glo[X] + (cos(a)*sin(b)*sin(y) - sin(a)*cos(y))*glo[Y] + (cos(a)*sin(b)*cos(y) + sin(a)*sin(y))*glo[Z];
    roc[Y] = sin(a)*cos(b)*glo[X] + (sin(a)*sin(b)*sin(y) + cos(a)*cos(y))*glo[Y] + (sin(a)*sin(b)*cos(y) - cos(a)*sin(y))*glo[Z];
    roc[Z] = -sin(b)*glo[X] + cos(b)*sin(y)*glo[Y] + cos(b)*cos(y)*glo[Z]; 
  }
  void Dynamics::updateNetF(double** forces)
  {
    for (int i = 0; i < NUM_FORCES; i++) {
      for (int j = 0; j < DIM; j++) {
        netF[j] += forces[i][j]; 
      }
    }

  }
  void Dynamics::updateNetM(double** moments)
  {
    for (int i = 0; i < NUM_FORCES; i++) {
      for (int j = 0; j < DIM; j++) {
        netM[j] += moments[i][j]; 
      }
    }
  } 

  void Dynamics::clearNet()
  {
    for (int i = 0; i < DIM; i++){
      netF[i] = 0;
      netM[i] = 0;
    }
  }
  void Dynamics::updateAcc()
  {
    for (int i = 0; i < DIM; i++){
      s->acc[i] = netF[i]/r->getMass();
    }
  }
  void Dynamics::updateVel(double dt) {
    for (int i = 0; i < DIM; i++)
    {
      s->vel[i] += dt*s->acc[i];
    }
  }
  double* Dynamics::getVel() {
    return s->vel;
  }
  double* Dynamics::getAcc() {
    return s->acc;
  }
  double* Dynamics::getAngVel() {
    return s->ang_vel;
  }
  double* Dynamics::getAngAcc() {
    return s->ang_acc;
  }
  void Dynamics::updatePos(double dt) {
    for (int i = 0; i < DIM; i++)
    {
      s->pos[i] += dt*s->vel[i];
    }
  }

  void Dynamics::updateVel(double dt, double acc[]) {
    for (int i = 0; i < DIM; i++)
    {
      s->vel[i] += dt*acc[i];
    }
  }
  void Dynamics::updatePos(double dt, double vel[]) {
    for (int i = 0; i < DIM; i++)
    {
      s->pos[i] += dt*vel[i];
    }
  }

  void Dynamics::updateAngAcc()
  {
    double I_Glo[DIM];
    rocketToGlobalFrame(r->getI(), I_Glo);
    for (int i = 0; i < DIM; i++){
      s->ang_acc[i] = (netM[i]- inerM[i])/(I_Glo[i]);
    }
  }
  void Dynamics::updateAngVel(double dt) {
    for (int i = 0; i < DIM; i++)
    {
      s->ang_vel[i] += dt*s->ang_acc[i];
    }
  }
  void Dynamics::updateAng(double dt) {
    for (int i = 0; i < DIM; i++)
    {
      s->ang[i] += dt*s->ang_vel[i];
    }
  }
  void Dynamics::updateAng(double dt, double ang_vel[]) {
    for (int i = 0; i < DIM; i++)
    {
      s->ang[i] += dt*ang_vel[i];
    }
  }

Dynamics::State* initState(Rocket* r)
{
  
  double* ang = new double[DIM] {0,0,0}; //global frame
  double* pos = new double[DIM] {0,0,0}; 
  double* vel = new double[DIM] {0,0,0};
  double* ang_vel = new double[DIM] {0,0,0};
  double* acc = new double[DIM] {0,0,0};
  double* ang_acc = new double[DIM] {0,0,0};

  Dynamics::State* s = new Dynamics::State(ang,pos,vel,ang_vel,acc,ang_acc,r->getMass());
  return s;
}
Shape** initStructure()
{
  int length = 10; //meters X
  int height = 10; //meters Y
  int width = 10; //meters Z
  float* P = new float[DIM] {0, 0, 0}; //origin of shape is same as origin of rocket
  double* COV = new double[DIM] {0., 0, height/2.};
  float* surfaceArea = new float[DIM] {(float)(height*width), (float)(length*width), (float)(length*height)};

  Shape** shapes = new Shape*[1] {new RectagularPrism(surfaceArea,COV,P,length*width*height,length,width,height)};
  return shapes;
}
Tank** initTanks() 
{
  int length = 10; //meters X
  int height = 10; //meters Y
  int width = 10; //meters Z
  float* P = new float[DIM] {0, 0, 0}; //origin of tank is same as origin of rocket
  double* COV = new double[DIM] {0, 0, height/2.};
  float* surfaceArea = new float[DIM] {(float)(height*width), (float)(length*width), (float)(length*height)};

  Shape* tankShape = new RectagularPrism(surfaceArea,COV,P,length*width*height,length,width,height);
  float dry_mass = 1000; //kg
  float massF = 10; //kg
  double mass_vel_out = .1; //kg/s
  double fuelD = 0.657; //kg/m^3
  float volFill = 1; //m^3
  float fuelHeight = .5; //meters
  
  Tank** tanks = new Tank*[1] {new Tank(tankShape, dry_mass, massF, mass_vel_out, fuelD, volFill, fuelHeight)}; 
  return tanks;
}
Rocket* initRocket()
{
  Shape** shapes = initStructure();
  int* massParts = new int[1] {1000};
  double* cov = (shapes[0])->getCOV();
  float* coeffDrag = new float[DIM] {1.28, 1.28, 1.28}; //https://www.grc.nasa.gov/www/k-12/airplane/shaped.html
  int d_m = 1000; //kg
  float* surfaceArea = shapes[0]->getSA();
  float m_vel_out = .1; //m/s set as a parameters before, not actually used for anything currently, mass changes are handled in the tanks
  int numTanks = 1;
  int numShape = 1;
  float max_t = MAX_THRUST; //N
  Tank** t = initTanks();

  Rocket* r = new Rocket(shapes, cov, coeffDrag, d_m, surfaceArea, m_vel_out, t, numTanks, numShape, massParts, max_t); 
  return r;
}


Environment* initEnvironment()
{
  return (new Environment(new double[DIM] {0,0,0}));
}
Control* initControl()
{
  return (new Control(new double[DIM] {0,0,0}));
}
void sampleData(Dynamics* sys, double timeElapsed)
{
  if (timeElapsed > 60*80) {
    cout << "CHECK VALUES" << "\n";
  }
  cout << sys->to_string() << "\n";
}
int main() {
  double timeElapsed = 0;
  
  double secsPerSample = dt*loopPerSample;
  double secs = secsPerSample;
  Rocket* r = initRocket();
  Dynamics::State* s = initState(r);
  Dynamics* sys = new Dynamics(r, initEnvironment(), initControl(), s);

  while (timeElapsed < TIME_FINAL) {
    sys->euler(dt);

    timeElapsed += dt;
    if (secs >= secsPerSample) 
    {
      sampleData(sys, timeElapsed); 
      secs = dt; 
    }
    else { secs += dt; }
  }
  cout << "DONE";

  //RENDERING POSITIONAL DATA
/*
// Initialization
    //--------------------------------------------------------------------------------------
    const int screenWidth = 800;
    const int screenHeight = 450;

    InitWindow(screenWidth, screenHeight, "raylib [core] example - 3d camera free");

    // Define the camera to look into our 3d world
    Camera3D camera = { 0 };
    camera.position = (Vector3){ 10.0f, 10.0f, 10.0f }; // Camera position
    camera.target = (Vector3){ 0.0f, 0.0f, 0.0f };      // Camera looking at point
    camera.up = (Vector3){ 0.0f, 1.0f, 0.0f };          // Camera up vector (rotation towards target)
    camera.fovy = 45.0f;                                // Camera field-of-view Y
    camera.projection = CAMERA_PERSPECTIVE;                   // Camera mode type

    Vector3 cubePosition = { 0.0f, 0.0f, 0.0f };

    SetCameraMode(camera, CAMERA_FREE); // Set a free camera mode

    SetTargetFPS(60);                   // Set our game to run at 60 frames-per-second
    //--------------------------------------------------------------------------------------

    // Main game loop
    while (!WindowShouldClose())        // Detect window close button or ESC key
    {
        // Update
        //----------------------------------------------------------------------------------
        UpdateCamera(&camera);          // Update camera

        if (IsKeyDown('Z')) camera.target = (Vector3){ 0.0f, 0.0f, 0.0f };
        //----------------------------------------------------------------------------------

        // Draw
        //----------------------------------------------------------------------------------
        BeginDrawing();

            ClearBackground(RAYWHITE);

            BeginMode3D(camera);

                DrawCube(cubePosition, 2.0f, 2.0f, 2.0f, RED);
                DrawCubeWires(cubePosition, 2.0f, 2.0f, 2.0f, MAROON);

                DrawGrid(10, 1.0f);

            EndMode3D();

            DrawRectangle( 10, 10, 320, 133, Fade(SKYBLUE, 0.5f));
            DrawRectangleLines( 10, 10, 320, 133, BLUE);

            DrawText("Free camera default controls:", 20, 20, 10, BLACK);
            DrawText("- Mouse Wheel to Zoom in-out", 40, 40, 10, DARKGRAY);
            DrawText("- Mouse Wheel Pressed to Pan", 40, 60, 10, DARKGRAY);
            DrawText("- Alt + Mouse Wheel Pressed to Rotate", 40, 80, 10, DARKGRAY);
            DrawText("- Alt + Ctrl + Mouse Wheel Pressed for Smooth Zoom", 40, 100, 10, DARKGRAY);
            DrawText("- Z to zoom to (0, 0, 0)", 40, 120, 10, DARKGRAY);

        EndDrawing();
        //----------------------------------------------------------------------------------
    }

    // De-Initialization
    //--------------------------------------------------------------------------------------
    CloseWindow();        // Close window and OpenGL context
    //--------------------------------------------------------------------------------------

    return 0;
*/
}
